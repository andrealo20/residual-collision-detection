/**
 * @file friction.h
 * @brief Joint friction, applied to the plant and hidden from the observer.
 *
 * Up to this point the plant and the observer share the same model, so the
 * residual is exactly the external torque and any threshold works. That is
 * not the situation on a real arm. Friction is present, it is imperfectly
 * known, and whatever the model gets wrong appears in the residual as if it
 * were a contact.
 *
 * This module exists to be given to the plant and withheld from the
 * observer. That asymmetry is the experiment: it is what turns the threshold
 * from a free parameter into a choice with consequences on both sides.
 *
 * ## The model
 *
 * A Stribeck friction curve, which is the standard form in the robot
 * identification literature:
 *
 *     tau_f(qd) = -[ ( Fc + (Fs - Fc) exp(-|qd/vs|) ) sign(qd)  +  Fv qd ]
 *
 * with Fc the Coulomb level, Fs the higher static level at breakaway, vs the
 * Stribeck velocity over which the transition happens, and Fv the viscous
 * coefficient. The exponential term is what produces the characteristic dip
 * as the joint starts moving.
 *
 * The leading minus sign is not decoration. The bracket is the magnitude
 * opposing motion, which is the form usually written on the left of the
 * equation of motion and subtracted; qcollide_friction_torque() returns the
 * negated version, so a caller adds it to the joint torque. The two
 * conventions differ by exactly this sign and getting it backwards produces a
 * plant that gains energy, so the formula is written the way the code returns
 * it. The Stribeck exponent is 1 here, not a free parameter.
 *
 * The parameter values used in the experiments are representative of
 * collaborative-robot joints rather than identified from a specific arm, and
 * the README says so. What matters for the measurement is that the friction
 * is a realistic *shape* the observer does not know about, not that it
 * reproduces a particular robot.
 *
 * ## Why sign(qd) is regularised
 *
 * The true Coulomb term is discontinuous at zero velocity, which is both
 * physically questionable and numerically hostile: an integrator stepping
 * across zero velocity sees an instantaneous torque reversal and the step
 * size collapses. A hyperbolic tangent over a small velocity band gives the
 * same behaviour away from zero without the discontinuity. The band is a
 * parameter so its effect can be measured rather than assumed negligible.
 */
#ifndef QCOLLIDE_FRICTION_H
#define QCOLLIDE_FRICTION_H

#include <stddef.h>

#include "qcollide/dynamics.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double coulomb;    /**< Fc, Nm */
    double stiction;   /**< Fs, Nm; at least Fc */
    double viscous;    /**< Fv, Nm s/rad */
    double stribeck;   /**< vs, rad/s; velocity scale of the transition */
    double smoothing;  /**< rad/s; width of the tanh band around zero */
} qcollide_friction;

/**
 * Uniform friction on every joint, with values representative of a small
 * collaborative arm.
 *
 * @param scale multiplies all torque terms, so a single number sweeps the
 *              severity of the model error in the experiments.
 */
void qcollide_friction_typical(qcollide_friction *f, double scale);

/**
 * Friction torque opposing motion at one joint.
 *
 * Returns the torque the friction applies, so the caller adds it to the
 * joint torque rather than subtracting it. The sign convention is checked in
 * the tests, since getting it backwards produces a plant that gains energy
 * instead of losing it, which looks like an unstable controller rather than
 * a sign error.
 */
double qcollide_friction_torque(const qcollide_friction *f, double qd);

/** Friction torque for every joint of an n-joint arm. */
void qcollide_friction_vector(const qcollide_friction *f, size_t n,
                              const double *qd, double *tau_f);

/* ------------------------------------------------------------------ */
/* Inertial parameter error                                            */
/* ------------------------------------------------------------------ */

/**
 * Copy a model with every mass and inertia scaled by (1 + error).
 *
 * The second source of model error. Inertial parameters come from CAD or
 * from identification and are wrong by a few percent at best; the residual
 * carries that error whenever the arm accelerates. Applying a uniform
 * relative error is the simplest form that still produces the effect, and
 * makes the severity a single number.
 *
 * @param error relative error, e.g. 0.05 for five percent heavy.
 */
qcollide_status qcollide_model_perturb(qcollide_model *dst,
                                       const qcollide_model *src,
                                       double error);

#ifdef __cplusplus
}
#endif

#endif /* QCOLLIDE_FRICTION_H */
