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
    double tau[N], ext[N], fr[N], applied[N];
    double qm[N], qdm[N], taum[N];
    const int steps = (int)(duration / DT + 0.5);
    int k, i, fired = 0;
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

    commanded_torque(&plant, 0.0, q, qd, tau);
    (void)qcollide_observer_reset(&obs, q, qd);

    for (k = 0; k < steps; ++k) {
        const double t = (double)k * DT;

        commanded_torque(&plant, t, q, qd, tau);
        qcollide_contact_torque(contact, N, t, ext);

        /* Friction acts on the plant only. The observer is given the
           commanded torque, exactly as a real controller would report it,
           and never learns that friction removed some of it. */
        qcollide_friction_vector(&fric, N, qd, fr);
        for (i = 0; i < N; ++i) {
            applied[i] = ext[i] + fr[i];
        }

        (void)qcollide_step(&plant, q, qd, tau, applied, DT);

        /* The observer sees measurements, not the plant state. This is the
           only place the two are distinguished, and keeping it to one place
           is what makes the ideal-sensor baseline exactly comparable. */
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

static void sweep_full(const char *title, double friction, double param_error,
                       int sensors, double torque_noise,
                       qcollide_velocity_mode vel_mode)
{
    const double thresholds[7] = { 0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0 };
    config cfg;
    int i;

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

        cfg.threshold = thresholds[i];
        (void)run(&cfg, NULL, 4.5, &events, NULL);
        if (events == 0u) {
            mdt = minimum_detectable(&cfg);
        }

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

int main(void)
{
    config probe;
    double peak;
    int i;

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
            printf("| friction x%.2f | %.3f Nm |\n", frictions[i], peak);
        }
    }
    {
        const double errors[4] = { 0.0, 0.02, 0.05, 0.10 };
        for (i = 0; i < 4; ++i) {
            probe.friction_scale = 0.0;
            probe.param_error = errors[i];
            (void)run(&probe, NULL, 3.0, NULL, &peak);
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
    return 0;
}
