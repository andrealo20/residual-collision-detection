/**
 * @file observer.h
 * @brief Momentum-based residual observer for sensorless collision detection.
 *
 * The residual is an estimate of the external joint torque acting on the arm,
 * built from measured quantities only. When nothing is touching the robot it
 * settles at zero; when something is, it settles at the external torque.
 *
 * ## Why momentum
 *
 * Reading the external torque straight off the equation of motion,
 *
 *     tau_ext = tau - M(q) qdd - C(q,qd) qd - g(q)
 *
 * requires joint acceleration, which is not measured. Obtaining it means
 * differentiating position twice, and each differentiation multiplies
 * high-frequency noise by the sample rate. At 1 kHz the resulting noise
 * swamps the contact one is trying to detect.
 *
 * The generalized momentum p = M(q) qd has the derivative
 *
 *     pdot = tau + C(q,qd)' qd - g(q) + tau_ext
 *
 * with no acceleration term. That identity is a consequence of Mdot - 2C
 * being skew-symmetric, which is why dynamics.h is careful about which
 * Coriolis matrix it produces; the property is verified numerically in the
 * dynamics tests rather than assumed.
 *
 * Defining the residual as
 *
 *     r = K ( p - integral( tau + C' qd - g + r ) dt - p(0) )
 *
 * gives, for each joint independently,
 *
 *     rdot = K ( tau_ext - r )
 *
 * a first-order low-pass response with time constant 1/K. So r tracks
 * tau_ext with a settling time set by one gain per joint, and nothing else.
 *
 * ## What the gain trades
 *
 * A large K detects faster but passes more of whatever the model gets wrong:
 * friction, inertial parameter error, torque measurement noise. That
 * trade-off is the subject of M3 and is measured there rather than argued
 * about here.
 *
 * ## Integration
 *
 * The integral is evaluated with the trapezoidal rule rather than forward
 * Euler. That is not a detail: the integrand contains pdot, which grows with
 * joint velocity, so a first-order rule leaves an error proportional to the
 * timestep times the speed of the arm. Measured on a two-link arm driven to
 * 19 rad/s, forward Euler produced a spurious residual of 0.33 Nm on a joint
 * with no external torque at all -- large enough to be mistaken for a light
 * contact. See docs/design.md.
 *
 * ## What this is not
 *
 * The residual estimates external torque at the *joints*. Recovering the
 * contact force and its location on the link needs more information than a
 * joint-torque vector carries: an n-joint arm gives n numbers, and a contact
 * has six degrees of freedom plus an unknown position along the structure.
 * Detection is well posed; localisation is not, without extra assumptions.
 */
#ifndef QCOLLIDE_OBSERVER_H
#define QCOLLIDE_OBSERVER_H

#include <stddef.h>

#include "qcollide/dynamics.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    const qcollide_model *model;

    double gain[QCOLLIDE_MAX_DOF];      /**< K per joint, 1/s */
    double integral[QCOLLIDE_MAX_DOF];  /**< running integral term */
    double residual[QCOLLIDE_MAX_DOF];  /**< current estimate of tau_ext */
    double p0[QCOLLIDE_MAX_DOF];        /**< momentum at initialisation */
    double prev[QCOLLIDE_MAX_DOF];      /**< previous integrand, for the
                                             trapezoidal rule */
    int    started;                     /**< p0 captured */
    int    have_prev;                   /**< prev is populated */
} qcollide_observer;

/**
 * Initialise the observer against a model and a per-joint gain.
 *
 * The model given here is the observer's *belief* about the arm. Passing the
 * same model the plant uses gives a perfect-model observer, which is what M1
 * tests; passing a perturbed one is how M3 measures the effect of parameter
 * error. Keeping them separate objects is what makes that experiment
 * possible at all.
 *
 * @param gain n elements, each strictly positive. A gain of K gives a
 *             residual time constant of 1/K seconds.
 * @return QCOLLIDE_OK, or QCOLLIDE_ERR_ARG on a null pointer or a
 *         non-positive gain.
 */
qcollide_status qcollide_observer_init(qcollide_observer *o,
                                       const qcollide_model *model,
                                       const double *gain);

/** Convenience: the same gain on every joint. */
qcollide_status qcollide_observer_init_uniform(qcollide_observer *o,
                                               const qcollide_model *model,
                                               double gain);

/**
 * Reset the observer to the current state.
 *
 * Called at start-up and after any discontinuity the observer should not
 * interpret as contact. The residual is zeroed and the momentum reference is
 * recaptured, so the estimate restarts from "nothing is touching the arm".
 */
qcollide_status qcollide_observer_reset(qcollide_observer *o,
                                        const double *q, const double *qd);

/**
 * Advance the observer by one sample.
 *
 * @param q, qd  measured joint position and velocity.
 * @param tau    commanded or measured joint torque.
 * @param dt     sample period, seconds.
 *
 * On the first call after init the momentum reference is captured and the
 * residual stays zero: there is no history to integrate yet, and returning a
 * spurious spike would be indistinguishable from a contact at start-up.
 */
qcollide_status qcollide_observer_update(qcollide_observer *o,
                                         const double *q, const double *qd,
                                         const double *tau, double dt);

/** Current residual, n elements. Estimates the external joint torque. */
const double *qcollide_observer_residual(const qcollide_observer *o);

/** Largest absolute residual across joints, the usual detection statistic. */
double qcollide_observer_magnitude(const qcollide_observer *o);

#ifdef __cplusplus
}
#endif

#endif /* QCOLLIDE_OBSERVER_H */
