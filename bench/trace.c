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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define N 2
#define DT 1e-4
static const double PI = 3.14159265358979323846;

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
    double q[N] = { 0.4, -0.3 }, qd[N] = { 0.0, 0.0 };
    double tau[N], qdd[N], ext[N], fr[N], applied[N];
    double qm[N], qdm[N], taum[N];
    const double w = 2.0 * PI * 0.4;
    int k, i;

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

        for (i = 0; i < N; ++i) {
            qdd[i] = -0.6 * w * w * sin(w * t + (double)i * 0.7);
        }
        (void)qcollide_inverse_dynamics(&plant, q, qd, qdd, tau);
        qcollide_contact_torque(&contact, N, t, ext);
        qcollide_friction_vector(&fric, N, qd, fr);
        for (i = 0; i < N; ++i) applied[i] = ext[i] + fr[i];

        (void)qcollide_step(&plant, q, qd, tau, applied, DT);
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
