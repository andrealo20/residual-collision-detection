/**
 * @file sensors.h
 * @brief Encoder quantisation, velocity estimation and torque noise.
 *
 * M3 gave the observer a wrong model. This gives it wrong *measurements*,
 * which is the other half of what separates a simulation from a robot.
 *
 * ## The velocity problem
 *
 * The momentum observer was chosen precisely to avoid joint acceleration,
 * because obtaining it means differentiating position twice and each
 * differentiation multiplies quantisation noise by the sample rate. But the
 * observer still needs *velocity*, and on a real arm velocity is not
 * measured either: it comes from differencing encoder counts.
 *
 * That single differentiation is not free. With a 17-bit encoder at 10 kHz,
 * one count of quantisation becomes
 *
 *     2*pi / 2^17 * 10000  =  0.48 rad/s
 *
 * of velocity noise from a backward difference -- larger than the velocities
 * of interest in slow motion. The usual answer is to filter, which trades
 * noise for lag, and lag in the velocity estimate appears in the residual as
 * a torque that is not there.
 *
 * So this module offers three velocity estimators, and the point of M4 is to
 * measure what each costs in detection sensitivity rather than to assert
 * that one is best.
 *
 * ## Determinism
 *
 * The noise generator is a small explicit PRNG with a caller-supplied seed
 * rather than rand(), so a result is reproducible across platforms and runs.
 * A sweep whose numbers move between machines cannot be compared against
 * itself.
 */
#ifndef QCOLLIDE_SENSORS_H
#define QCOLLIDE_SENSORS_H

#include <stddef.h>
#include <stdint.h>

#include "qcollide/dynamics.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Deterministic Gaussian source. xorshift plus Box-Muller. */
typedef struct {
    uint64_t state;
    double   spare;
    int      have_spare;
} qcollide_noise;

void   qcollide_noise_seed(qcollide_noise *n, uint64_t seed);
double qcollide_noise_gaussian(qcollide_noise *n);

/** Velocity estimation strategy. */
typedef enum {
    QCOLLIDE_VEL_EXACT = 0,   /**< the plant's true velocity: the baseline */
    QCOLLIDE_VEL_DIFF,        /**< backward difference of quantised position */
    QCOLLIDE_VEL_FILTERED     /**< backward difference, first-order low pass */
} qcollide_velocity_mode;

typedef struct {
    /** Encoder resolution in counts per revolution; 0 disables quantisation. */
    uint32_t encoder_counts;

    /** Standard deviation of torque measurement noise, Nm. */
    double torque_noise;

    qcollide_velocity_mode velocity_mode;

    /** Cut-off of the velocity filter, Hz. Used only in FILTERED mode. */
    double velocity_cutoff;
} qcollide_sensor_config;

typedef struct {
    qcollide_sensor_config cfg;
    qcollide_noise         noise;

    double last_q[QCOLLIDE_MAX_DOF];   /**< previous quantised position */
    double vel[QCOLLIDE_MAX_DOF];      /**< current velocity estimate */
    int    have_last;
} qcollide_sensors;

/** Defaults representative of an industrial joint: 17-bit encoder, 1 % torque noise. */
void qcollide_sensors_typical(qcollide_sensor_config *cfg, double torque_scale);

qcollide_status qcollide_sensors_init(qcollide_sensors *s,
                                      const qcollide_sensor_config *cfg,
                                      uint64_t seed);

/** Clear the difference history, so the next sample starts a new estimate. */
void qcollide_sensors_reset(qcollide_sensors *s);

/**
 * Produce measured position, velocity and torque from the true values.
 *
 * @param q_true, qd_true, tau_true  plant state and commanded torque.
 * @param q_meas, qd_meas, tau_meas  outputs; any may alias its input.
 * @param dt  sample period, needed for the velocity difference.
 *
 * In EXACT velocity mode qd_meas is copied from qd_true, which is the
 * baseline the other modes are measured against rather than a claim that
 * velocity is ever known exactly.
 */
qcollide_status qcollide_sensors_measure(qcollide_sensors *s, size_t n,
                                         const double *q_true,
                                         const double *qd_true,
                                         const double *tau_true,
                                         double *q_meas, double *qd_meas,
                                         double *tau_meas, double dt);

#ifdef __cplusplus
}
#endif

#endif /* QCOLLIDE_SENSORS_H */
