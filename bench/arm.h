/**
 * @file arm.h
 * @brief The trajectory and controller shared by the sweep and the trace.
 *
 * Both programs used to command the arm with feedforward alone, the inverse
 * dynamics of a sinusoidal acceleration, starting from rest. The sinusoid's
 * velocity is not zero at t = 0, so the arm began off its trajectory with a
 * constant velocity error that nothing corrected: without friction both
 * joints wound up through several radians, and with friction the second
 * joint stuck and never moved. Every number measured on top of that was
 * measured on a motion nobody described.
 *
 * The arm now starts on the trajectory and tracks it with computed-torque
 * control: feedforward from the plant model plus a PD correction, so friction
 * and contacts are rejected instead of accumulating. The controller reads the
 * measured state, as a real one would. The observer is untouched by this
 * choice: it is given whatever torque was commanded.
 */
#ifndef QCOLLIDE_BENCH_ARM_H
#define QCOLLIDE_BENCH_ARM_H

#include "qcollide/dynamics.h"

#include <math.h>

#define ARM_N 2

/* 0.4 Hz, 0.6 rad about a fixed centre, the joints 0.7 rad apart in phase. */
#define ARM_FREQ_HZ 0.4
#define ARM_AMP     0.6
#define ARM_LAG     0.7

/* Tracking bandwidth 40 rad/s, critically damped. The residual does not
   depend on the gains, since the observer is given the torque actually
   commanded; they only decide how closely the arm follows the trajectory. */
#define ARM_KP 1600.0
#define ARM_KD 80.0

static const double arm_centre[ARM_N] = { 0.4, -0.3 };

/** Desired position, velocity and acceleration at time t. */
static inline void arm_desired(double t, double *q, double *qd, double *qdd)
{
    const double w = 2.0 * 3.14159265358979323846 * ARM_FREQ_HZ;
    int i;
    for (i = 0; i < ARM_N; ++i) {
        const double ph = w * t + (double)i * ARM_LAG;
        q[i] = arm_centre[i] + ARM_AMP * sin(ph);
        qd[i] = ARM_AMP * w * cos(ph);
        qdd[i] = -ARM_AMP * w * w * sin(ph);
    }
}

/** Computed-torque command from the measured state qm, qdm. */
static inline void arm_command(const qcollide_model *plant, double t,
                               const double *qm, const double *qdm,
                               double *tau)
{
    double q[ARM_N], qd[ARM_N], a[ARM_N];
    int i;
    arm_desired(t, q, qd, a);
    for (i = 0; i < ARM_N; ++i) {
        a[i] += ARM_KP * (q[i] - qm[i]) + ARM_KD * (qd[i] - qdm[i]);
    }
    (void)qcollide_inverse_dynamics(plant, qm, qdm, a, tau);
}

#endif /* QCOLLIDE_BENCH_ARM_H */
