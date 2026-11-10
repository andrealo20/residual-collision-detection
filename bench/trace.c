/**
 * @file trace.c
 * @brief Dump a residual trace to CSV, for the figures in the README.
 *
 * One run of the full chain with a contact in the middle, under a chosen
 * combination of model error and sensor imperfection. The figures are
 * generated from this rather than drawn, so they can be regenerated when the
 * numbers change and cannot drift away from them.
 *
 * Usage: qcollide_trace [friction_scale] [param_error] [sensors 0|1]
 */
#include "qcollide/detector.h"
#include "qcollide/friction.h"
#include "qcollide/sensors.h"

#include "arm.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define N ARM_N
#define DT 1e-4


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

int main(int argc, char **argv)
{
    const double friction_scale = (argc > 1) ? atof(argv[1]) : 0.0;
    const double param_error    = (argc > 2) ? atof(argv[2]) : 0.0;
    const int    use_sensors    = (argc > 3) ? atoi(argv[3]) : 0;

    qcollide_model plant, belief;
    qcollide_friction fric;
    qcollide_observer obs;
    qcollide_sensors sens;
    qcollide_sensor_config scfg;
    qcollide_contact contact;
    double q[N], qd[N], qdd[N];
    double tau[N], ext[N];
    double qm[N], qdm[N], taum[N];
    int k, i;

    /* Start on the trajectory; see arm.h for why that matters. */
    arm_desired(0.0, q, qd, qdd);
    for (i = 0; i < N; ++i) {
        qm[i] = q[i];
        qdm[i] = qd[i];
    }

    (void)qcollide_model_uniform(&plant, N, 1.0, 0.6, 9.81);
    (void)qcollide_model_perturb(&belief, &plant, param_error);
    qcollide_friction_typical(&fric, friction_scale);
    (void)qcollide_observer_init_uniform(&obs, &belief, 50.0);

    qcollide_sensors_typical(&scfg, 1.0);
    if (!use_sensors) {
        scfg.encoder_counts = 0u;
        scfg.torque_noise = 0.0;
        scfg.velocity_mode = QCOLLIDE_VEL_EXACT;
    }
    (void)qcollide_sensors_init(&sens, &scfg, 20260801u);

    contact.start = 1.5;
    contact.duration = 0.4;
    contact.torque[0] = 5.0;
    contact.torque[1] = 0.0;

    printf("t,residual0,residual1,magnitude,tau_ext0\n");

    for (k = 0; k < 30000; ++k) {
        const double t = (double)k * DT;

        arm_command(&plant, t, qm, qdm, tau);
        qcollide_contact_torque(&contact, N, t, ext);
        {
            plant_ext ctx;
            ctx.fric = &fric;
            ctx.contact = ext;
            (void)qcollide_step_ext(&plant, q, qd, tau,
                                    plant_ext_torque, &ctx, DT);
        }
        (void)qcollide_sensors_measure(&sens, N, q, qd, tau, qm, qdm, taum, DT);

        if (k == 0) {
            (void)qcollide_observer_reset(&obs, qm, qdm);
            continue;
        }
        (void)qcollide_observer_update(&obs, qm, qdm, taum, DT);

        /* Decimated to keep the CSV small; the residual has no content above
           a few hundred hertz once the observer has filtered it. */
        if (k % 10 == 0) {
            const double *r = qcollide_observer_residual(&obs);
            printf("%.4f,%.6f,%.6f,%.6f,%.3f\n",
                   t, r[0], r[1], qcollide_observer_magnitude(&obs), ext[0]);
        }
    }
    return 0;
}
