/**
 * @file sweep.c
 * @brief The threshold trade-off under unmodelled friction and parameter error.
 *
 * Everything before this ran with the observer holding the same model as the
 * plant, so the residual was exactly the external torque and any threshold
 * worked. This program removes that assumption, which is the only interesting
 * version of the question.
 *
 * Two independent sources of error, each a single severity number:
 *
 *   - **friction** present in the plant and absent from the observer;
 *   - **inertial parameter error**, the observer's masses and inertias
 *     scaled away from the plant's.
 *
 * For each severity it measures the two quantities a threshold trades
 * between:
 *
 *   - **false alarms**: triggers over a run with no contact at all;
 *   - **minimum detectable torque**: the smallest contact that fires,
 *     by bisection.
 *
 * The useful output is the pair. A threshold low enough to catch a light
 * contact is also low enough to fire on friction, and the width of the gap
 * between them is what says whether the method is usable on a given arm.
 *
 * The arm is driven along a trajectory rather than held still, because
 * friction and inertial error only appear when the joints move. A detector
 * evaluated at rest would report a trade-off that does not exist.
 */
#define _POSIX_C_SOURCE 199309L

#include "qcollide/detector.h"
#include "qcollide/friction.h"
#include "qcollide/sensors.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#define N 2
#define DT 1e-4

static const double PI = 3.14159265358979323846;

typedef struct {
    double friction_scale;   /* 0 = none */
    double param_error;      /* 0 = perfect model */
    double gain;
    double threshold;
    unsigned confirm;
    int    sensors;          /* 0 = ideal measurements */
    double torque_noise;     /* Nm, standard deviation */
    qcollide_velocity_mode vel_mode;
} config;

/**
 * Commanded torque for a sinusoidal joint trajectory, from inverse dynamics
 * on the *plant* model.
 *
 * Using the plant model here is deliberate: it isolates the experiment. The
 * residual then carries only friction and the observer's parameter error,
 * not a tracking error from an imperfect controller as well. Adding a
 * controller would make the measurement depend on how it was tuned.
 */
static void commanded_torque(const qcollide_model *plant, double t,
                             const double *q, const double *qd, double *tau)
{
    double qdd[N];
    int i;
    const double w = 2.0 * PI * 0.4;      /* 0.4 Hz */
    const double amp = 0.6;               /* rad */

    for (i = 0; i < N; ++i) {
        const double phase = w * t + (double)i * 0.7;
        qdd[i] = -amp * w * w * sin(phase);
    }
    (void)qcollide_inverse_dynamics(plant, q, qd, qdd, tau);
}


/**
 * External torque seen by the plant: friction, which depends on the joint
 * velocity, plus the contact, which does not.
 *
 * Passed to qcollide_step_ext() so friction is re-evaluated at every
 * Runge-Kutta stage. Evaluating it once per step and holding it constant, as
 * this file used to, silently drops the integrator to first order in the
 * friction term, which is present in most of the cases swept below.
 */
typedef struct {
    const qcollide_friction *fric;
    const double *contact;   /* constant across the step */
} plant_ext;

static void plant_ext_torque(const double *q, const double *qd, double *tau_ext,
                             size_t n, void *user_data)
{
    const plant_ext *c = (const plant_ext *)user_data;
    size_t i;
    (void)q;
    qcollide_friction_vector(c->fric, n, qd, tau_ext);
    for (i = 0u; i < n; ++i) {
        tau_ext[i] += c->contact[i];
    }
}

/**
 * Run the chain once.
 *
 * @param contact  may be NULL for a run with no contact.
 * @return 1 if the detector fired at any point.
 */
static int run(const config *cfg, const qcollide_contact *contact,
               double duration, unsigned *events, double *peak_residual)
{
    qcollide_model plant, belief;
    qcollide_friction fric;
    qcollide_observer obs;
    qcollide_detector det;
    qcollide_sensors sens;
    qcollide_sensor_config scfg;
    double q[N] = { 0.4, -0.3 };
    double qd[N] = { 0.0, 0.0 };
    double tau[N], ext[N];
    double qm[N], qdm[N], taum[N];
    const int steps = (int)(duration / DT + 0.5);
    int k, fired = 0;
    double worst = 0.0;

    (void)qcollide_model_uniform(&plant, N, 1.0, 0.6, 9.81);
    (void)qcollide_model_perturb(&belief, &plant, cfg->param_error);
    qcollide_friction_typical(&fric, cfg->friction_scale);

    (void)qcollide_observer_init_uniform(&obs, &belief, cfg->gain);
    (void)qcollide_detector_init(&det, cfg->threshold,
                                 0.5 * cfg->threshold, cfg->confirm);

    qcollide_sensors_typical(&scfg, 1.0);
    scfg.torque_noise = cfg->torque_noise;
    scfg.velocity_mode = cfg->vel_mode;
    if (!cfg->sensors) {
        scfg.encoder_counts = 0u;
        scfg.torque_noise = 0.0;
        scfg.velocity_mode = QCOLLIDE_VEL_EXACT;
    }
    /* Fixed seed: a sweep whose numbers move between runs cannot be
       compared against itself. */
    (void)qcollide_sensors_init(&sens, &scfg, 20260801u);

    (void)qcollide_observer_reset(&obs, q, qd);

    for (k = 0; k < steps; ++k) {
        const double t = (double)k * DT;

        commanded_torque(&plant, t, q, qd, tau);
        qcollide_contact_torque(contact, N, t, ext);

        /* Friction acts on the plant only. The observer is given the
           commanded torque, exactly as a real controller would report it,
           and never learns that friction removed some of it. It is applied
           through the callback form so that it is re-evaluated at each
           Runge-Kutta stage rather than frozen at the start of the step. */
        {
            plant_ext ctx;
            ctx.fric = &fric;
            ctx.contact = ext;
            (void)qcollide_step_ext(&plant, q, qd, tau,
                                    plant_ext_torque, &ctx, DT);
        }

        /* The observer sees measurements, not the plant state. This is the
           only place the two are distinguished, and keeping it to one place
           is what makes the ideal-sensor baseline exactly comparable.

           Note the one sample offset: the state measured here is the state
           after the step, while tau is the torque that was applied across it.
           That is the pairing a real controller reports, the command it just
           issued alongside the encoder read that followed, so it is kept
           deliberately rather than corrected. */
        (void)qcollide_sensors_measure(&sens, N, q, qd, tau,
                                       qm, qdm, taum, DT);
        (void)qcollide_observer_update(&obs, qm, qdm, taum, DT);

        {
            const double mag = qcollide_observer_magnitude(&obs);
            if (t > 0.3 && mag > worst) worst = mag;   /* skip the start-up */
            if (t > 0.3 && qcollide_detector_update(&det, mag)) fired = 1;
        }
    }

    if (events != NULL) *events = qcollide_detector_events(&det);
    if (peak_residual != NULL) *peak_residual = worst;
    return fired;
}

/**
 * Smallest contact torque that fires, by bisection, as a worst case over
 * when the contact happens.
 *
 * Measuring at a single instant gives a number that depends on where the
 * arm was in its trajectory. Under friction that dependence is strong: a
 * contact arriving while the friction residual happens to be near its peak
 * needs almost nothing to push the total over the threshold, while the same
 * contact half a cycle later needs far more. The first version of this
 * function used a fixed start time and reported the lucky case.
 *
 * The contact is therefore tried at several phases of the trajectory and the
 * hardest one is reported, since that is the contact the detector has to be
 * able to catch.
 */
static double minimum_detectable(const config *cfg)
{
    /* Spread across one 2.5 s period of the 0.4 Hz trajectory. */
    static const double starts[5] = { 1.00, 1.50, 2.00, 2.50, 3.00 };
    double worst = -1.0;
    int phase;

    for (phase = 0; phase < 5; ++phase) {
        qcollide_contact c;
        double lo = 0.0, hi = 40.0;
        int it, i;

        c.start = starts[phase];
        c.duration = 1.0;
        for (i = 0; i < N; ++i) c.torque[i] = 0.0;

        c.torque[0] = hi;
        if (!run(cfg, &c, c.start + 1.5, NULL, NULL)) {
            return -1.0;   /* unbounded at some phase: report it as such */
        }

        for (it = 0; it < 18; ++it) {
            const double mid = 0.5 * (lo + hi);
            c.torque[0] = mid;
            if (run(cfg, &c, c.start + 1.5, NULL, NULL)) hi = mid; else lo = mid;
        }
        if (hi > worst) worst = hi;
    }
    return worst;
}

/* ------------------------------------------------------------------ */
/* Collected results                                                   */
/*                                                                     */
/* Everything printed below is also stored, so that --json can write it */
/* out. tools/make_figures.py and tools/check_readme.py read that file  */
/* rather than carrying transcribed copies of these numbers, which is   */
/* what they used to do and what let the figures drift away from the    */
/* measurements they claimed to show.                                  */
/* ------------------------------------------------------------------ */

#define MAX_THRESHOLDS 8
#define MAX_CASES 8
#define MAX_FLOORS 8

typedef struct {
    double threshold;
    unsigned events;
    double min_detectable;   /* negative when nothing was detected below the ceiling */
    int usable;              /* 0 when the detector fires with no contact present */
} threshold_row;

typedef struct {
    const char *title;
    double friction;
    double param_error;
    int sensors;
    threshold_row rows[MAX_THRESHOLDS];
    int n_rows;
} sweep_case;

typedef struct {
    double severity;
    double floor;
} floor_row;

typedef struct {
    const char *name;
    double floor;
} named_floor_row;

static sweep_case g_cases[MAX_CASES];
static int g_n_cases = 0;
static floor_row g_friction_floors[MAX_FLOORS];
static int g_n_friction = 0;
static floor_row g_param_floors[MAX_FLOORS];
static int g_n_param = 0;
static named_floor_row g_sensor_floors[MAX_FLOORS];
static int g_n_sensor = 0;

static void sweep_full(const char *title, double friction, double param_error,
                       int sensors, double torque_noise,
                       qcollide_velocity_mode vel_mode)
{
    const double thresholds[7] = { 0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0 };
    config cfg;
    sweep_case *rec = &g_cases[g_n_cases++];
    int i;

    rec->title = title;
    rec->friction = friction;
    rec->param_error = param_error;
    rec->sensors = sensors;
    rec->n_rows = 0;

    cfg.friction_scale = friction;
    cfg.param_error = param_error;
    cfg.gain = 50.0;
    cfg.confirm = 3;
    cfg.sensors = sensors;
    cfg.torque_noise = torque_noise;
    cfg.vel_mode = vel_mode;

    printf("\n## %s\n\n", title);
    printf("| threshold | false alarms in 3 s | min. detectable torque |\n");
    printf("|---|---|---|\n");

    for (i = 0; i < 7; ++i) {
        unsigned events = 0u;
        double mdt = -1.0;
        threshold_row *row = &rec->rows[rec->n_rows++];

        cfg.threshold = thresholds[i];
        (void)run(&cfg, NULL, 4.5, &events, NULL);
        if (events == 0u) {
            mdt = minimum_detectable(&cfg);
        }

        row->threshold = thresholds[i];
        row->events = events;
        row->min_detectable = mdt;
        row->usable = (events == 0u);

        printf("| %.2f Nm | %u | ", thresholds[i], events);
        if (events > 0u) {
            /* With the detector already firing on model error alone, the
               bisection converges to zero and would print 0.00 Nm, which
               reads as perfect sensitivity when it means the opposite: the
               threshold is unusable. */
            printf("unusable (fires without contact) |\n");
        } else if (mdt < 0.0) {
            printf("not detected below 40 Nm |\n");
        } else {
            printf("%.2f Nm |\n", mdt);
        }
    }
}

/* ------------------------------------------------------------------ */
/* Machine-readable output                                             */
/* ------------------------------------------------------------------ */

static void write_json(const char *path)
{
    FILE *f = fopen(path, "w");
    int c, i;

    if (f == NULL) {
        fprintf(stderr, "could not open %s for writing\n", path);
        return;
    }

    fprintf(f, "{\n");
    fprintf(f, "  \"gain\": 50.0,\n");
    fprintf(f, "  \"confirm_samples\": 3,\n");
    fprintf(f, "  \"dt\": %g,\n", DT);

    fprintf(f, "  \"friction_floors\": [");
    for (i = 0; i < g_n_friction; ++i) {
        fprintf(f, "%s\n    {\"scale\": %.6g, \"floor\": %.6f}",
                i ? "," : "", g_friction_floors[i].severity,
                g_friction_floors[i].floor);
    }
    fprintf(f, "\n  ],\n");

    fprintf(f, "  \"parameter_floors\": [");
    for (i = 0; i < g_n_param; ++i) {
        fprintf(f, "%s\n    {\"error\": %.6g, \"floor\": %.6f}",
                i ? "," : "", g_param_floors[i].severity,
                g_param_floors[i].floor);
    }
    fprintf(f, "\n  ],\n");

    fprintf(f, "  \"sensor_floors\": [");
    for (i = 0; i < g_n_sensor; ++i) {
        fprintf(f, "%s\n    {\"case\": \"%s\", \"floor\": %.6f}",
                i ? "," : "", g_sensor_floors[i].name,
                g_sensor_floors[i].floor);
    }
    fprintf(f, "\n  ],\n");

    fprintf(f, "  \"cases\": [");
    for (c = 0; c < g_n_cases; ++c) {
        const sweep_case *sc = &g_cases[c];
        fprintf(f, "%s\n    {\n      \"title\": \"%s\",\n", c ? "," : "", sc->title);
        fprintf(f, "      \"friction_scale\": %.6g,\n", sc->friction);
        fprintf(f, "      \"parameter_error\": %.6g,\n", sc->param_error);
        fprintf(f, "      \"sensors\": %s,\n", sc->sensors ? "true" : "false");
        fprintf(f, "      \"thresholds\": [");
        for (i = 0; i < sc->n_rows; ++i) {
            const threshold_row *r = &sc->rows[i];
            fprintf(f, "%s\n        {\"threshold\": %.4f, \"false_alarms\": %u, "
                       "\"usable\": %s, \"min_detectable\": ",
                    i ? "," : "", r->threshold, r->events,
                    r->usable ? "true" : "false");
            if (r->usable && r->min_detectable >= 0.0) {
                fprintf(f, "%.4f}", r->min_detectable);
            } else {
                fprintf(f, "null}");
            }
        }
        fprintf(f, "\n      ]\n    }");
    }
    fprintf(f, "\n  ]\n}\n");

    fclose(f);
    fprintf(stderr, "wrote %s\n", path);
}

int main(int argc, char **argv)
{
    config probe;
    double peak;
    const char *json_path = NULL;
    int i;

    for (i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--json") == 0 && i + 1 < argc) {
            json_path = argv[++i];
        } else {
            fprintf(stderr, "usage: %s [--json PATH]\n", argv[0]);
            return 2;
        }
    }

    printf("# qcollide - threshold trade-off under model error\n\n");
    printf("Two-link arm, sinusoidal trajectory at 0.4 Hz, observer gain "
           "K = 50 (time constant 20 ms), 3 consecutive samples to confirm.\n");

    /* Peak residual with no contact: the floor any threshold must clear. */
    printf("\n## Residual floor with no contact\n\n");
    printf("| model error | peak residual, no contact |\n|---|---|\n");

    probe.gain = 50.0;
    probe.threshold = 1e6;      /* never fires: measuring the residual only */
    probe.confirm = 3;
    probe.sensors = 0;
    probe.torque_noise = 0.0;
    probe.vel_mode = QCOLLIDE_VEL_EXACT;

    {
        const double frictions[4] = { 0.0, 0.25, 0.5, 1.0 };
        for (i = 0; i < 4; ++i) {
            probe.friction_scale = frictions[i];
            probe.param_error = 0.0;
            (void)run(&probe, NULL, 3.0, NULL, &peak);
            g_friction_floors[g_n_friction].severity = frictions[i];
            g_friction_floors[g_n_friction].floor = peak;
            g_n_friction++;
            printf("| friction x%.2f | %.3f Nm |\n", frictions[i], peak);
        }
    }
    {
        const double errors[4] = { 0.0, 0.02, 0.05, 0.10 };
        for (i = 0; i < 4; ++i) {
            probe.friction_scale = 0.0;
            probe.param_error = errors[i];
            (void)run(&probe, NULL, 3.0, NULL, &peak);
            g_param_floors[g_n_param].severity = errors[i];
            g_param_floors[g_n_param].floor = peak;
            g_n_param++;
            printf("| inertia +%.0f%% | %.3f Nm |\n", errors[i] * 100.0, peak);
        }
    }

    /* Sensor imperfection, measured on its own before being combined. */
    printf("\n## Residual floor from sensors alone (perfect model)\n\n");
    printf("| measurement | peak residual, no contact |\n|---|---|\n");
    probe.sensors = 1;
    {
        struct { const char *name; double noise; qcollide_velocity_mode m; }
        cases[5] = {
            { "exact velocity, no torque noise",   0.00, QCOLLIDE_VEL_EXACT },
            { "exact velocity, 0.05 Nm noise",     0.05, QCOLLIDE_VEL_EXACT },
            { "encoder difference, no filter",     0.00, QCOLLIDE_VEL_DIFF },
            { "encoder difference, 50 Hz filter",  0.00, QCOLLIDE_VEL_FILTERED },
            { "filtered velocity + 0.05 Nm noise", 0.05, QCOLLIDE_VEL_FILTERED },
        };
        for (i = 0; i < 5; ++i) {
            probe.friction_scale = 0.0;
            probe.param_error = 0.0;
            probe.torque_noise = cases[i].noise;
            probe.vel_mode = cases[i].m;
            (void)run(&probe, NULL, 3.0, NULL, &peak);
            g_sensor_floors[g_n_sensor].name = cases[i].name;
            g_sensor_floors[g_n_sensor].floor = peak;
            g_n_sensor++;
            printf("| %s | %.3f Nm |\n", cases[i].name, peak);
        }
    }

    sweep_full("Perfect model, ideal sensors (baseline)",
               0.0, 0.0, 0, 0.0, QCOLLIDE_VEL_EXACT);
    sweep_full("Friction only, ideal sensors",
               1.0, 0.0, 0, 0.0, QCOLLIDE_VEL_EXACT);
    sweep_full("Inertial parameters 5% heavy, ideal sensors",
               0.0, 0.05, 0, 0.0, QCOLLIDE_VEL_EXACT);
    sweep_full("Real sensors, perfect model",
               0.0, 0.0, 1, 0.05, QCOLLIDE_VEL_FILTERED);
    sweep_full("Everything: friction, 5% parameter error, real sensors",
               1.0, 0.05, 1, 0.05, QCOLLIDE_VEL_FILTERED);

    printf("\nThe pair of columns is the result. A threshold below the "
           "residual floor fires without any contact at all; one above it "
           "misses every contact smaller than itself. The lowest usable "
           "threshold, and the smallest contact it can still catch, are what "
           "decide whether the method is usable on a given arm.\n");
    printf("\nThe minimum detectable torque is a worst case over five "
           "contact phases within the trajectory. Measured at a single "
           "instant it is considerably lower and considerably less "
           "meaningful: under friction it depends on whether the contact "
           "happens to arrive while the model error is already pushing the "
           "residual upwards.\n");

    if (json_path != NULL) {
        write_json(json_path);
    }
    return 0;
}
