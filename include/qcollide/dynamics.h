/**
 * @file dynamics.h
 * @brief Rigid-body dynamics of a planar serial arm.
 *
 * This is the plant the collision observer is tested against. It is written
 * here rather than pulled from a dynamics library for one reason: the
 * observer's correctness is judged against this model, so if the model were
 * a black box a disagreement would be untraceable. Everything here is
 * checked against Pinocchio on identical inputs (tools/check_dynamics.py),
 * which is the same arrangement qdsp has with CMSIS-DSP.
 *
 * ## Model
 *
 * An n-link planar arm under gravity, revolute joints, each link a uniform
 * rod with mass m_i, length l_i and centre of mass at distance c_i from its
 * proximal joint. The equations of motion are
 *
 *     M(q) qdd + C(q, qd) qd + g(q) = tau + tau_ext
 *
 * with M symmetric positive definite. Planar keeps the algebra small enough
 * to write out and check by hand, which matters because this file is the
 * reference everything else is measured against.
 *
 * ## Why C(q, qd) is returned as a matrix, not a vector
 *
 * Most implementations return the Coriolis *vector* C(q,qd)*qd directly,
 * which is cheaper. The momentum observer in observer.h needs the matrix
 * itself, because it uses the transpose: the skew-symmetry of Mdot - 2C is
 * what makes the observer work without joint acceleration. Returning the
 * vector would have forced the observer to reconstruct the matrix by finite
 * differences, defeating the point.
 */
#ifndef QCOLLIDE_DYNAMICS_H
#define QCOLLIDE_DYNAMICS_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Maximum number of links. Raising it costs static memory and nothing else. */
#ifndef QCOLLIDE_MAX_DOF
#define QCOLLIDE_MAX_DOF 6
#endif

typedef enum {
    QCOLLIDE_OK        =  0,
    QCOLLIDE_ERR_ARG   = -1,   /**< null pointer or out-of-range dimension */
    QCOLLIDE_ERR_RANGE = -2    /**< valid arguments, numerically unsafe */
} qcollide_status;

/** One link: uniform rod, revolute joint at its proximal end. */
typedef struct {
    double mass;       /**< kg */
    double length;     /**< m, joint to joint */
    double com;        /**< m, joint to centre of mass along the link */
    double inertia;    /**< kg m^2 about the centre of mass, out of plane */
} qcollide_link;

typedef struct {
    qcollide_link link[QCOLLIDE_MAX_DOF];
    size_t        n;
    double        gravity;   /**< m/s^2, positive downwards */
} qcollide_model;

/**
 * Initialise a model from link parameters.
 *
 * @return QCOLLIDE_OK, or QCOLLIDE_ERR_ARG on a null pointer, n outside
 *         1..QCOLLIDE_MAX_DOF, a non-positive mass or length, a negative
 *         inertia, or a centre of mass outside [0, length].
 */
qcollide_status qcollide_model_init(qcollide_model *m,
                                    const qcollide_link *links,
                                    size_t n, double gravity);

/**
 * Build a uniform-rod arm: every link identical, centre of mass at the
 * midpoint, inertia m*l^2/12. Convenient for tests, where what matters is
 * that the model is exactly known rather than realistic.
 */
qcollide_status qcollide_model_uniform(qcollide_model *m, size_t n,
                                       double mass, double length,
                                       double gravity);

/**
 * Mass matrix M(q), row-major, n*n.
 *
 * Symmetric and positive definite for any configuration of a physical arm;
 * both properties are asserted in the tests, because a violation means the
 * model is wrong in a way that would otherwise show up much later as an
 * unstable simulation.
 */
qcollide_status qcollide_mass_matrix(const qcollide_model *m,
                                     const double *q, double *M);

/**
 * Coriolis matrix C(q, qd), row-major, n*n, in the Christoffel form.
 *
 * That specific choice matters: it is the one for which Mdot(q,qd) - 2C is
 * skew-symmetric. The observer relies on that identity, and the tests check
 * it numerically rather than taking it on trust.
 */
qcollide_status qcollide_coriolis_matrix(const qcollide_model *m,
                                         const double *q, const double *qd,
                                         double *C);

/** Gravity torque vector g(q), n elements. */
qcollide_status qcollide_gravity(const qcollide_model *m,
                                 const double *q, double *g);

/** Generalized momentum p = M(q) qd, n elements. */
qcollide_status qcollide_momentum(const qcollide_model *m,
                                  const double *q, const double *qd,
                                  double *p);

/** Total mechanical energy, kinetic plus potential. */
double qcollide_energy(const qcollide_model *m,
                       const double *q, const double *qd);

/**
 * Forward dynamics: joint accelerations from applied torques.
 *
 *     qdd = M(q)^-1 ( tau + tau_ext - C(q,qd) qd - g(q) )
 *
 * @param tau_ext optional external joint torque, may be NULL for none.
 * @return QCOLLIDE_ERR_RANGE if M is numerically singular, which for a
 *         physical arm means the parameters are wrong rather than the
 *         configuration being unusual.
 */
qcollide_status qcollide_forward_dynamics(const qcollide_model *m,
                                          const double *q, const double *qd,
                                          const double *tau,
                                          const double *tau_ext,
                                          double *qdd);

/**
 * Inverse dynamics: torque required for a commanded acceleration.
 *
 *     tau = M(q) qdd + C(q,qd) qd + g(q)
 */
qcollide_status qcollide_inverse_dynamics(const qcollide_model *m,
                                          const double *q, const double *qd,
                                          const double *qdd, double *tau);

/**
 * One fixed-step integration of the plant, fourth-order Runge-Kutta.
 *
 * RK4 rather than Euler because the tests check energy conservation, and
 * Euler injects energy at a rate that would swamp the effect being measured.
 * The integration error is quantified in the tests instead of assumed small.
 */
qcollide_status qcollide_step(const qcollide_model *m, double *q, double *qd,
                              const double *tau, const double *tau_ext,
                              double dt);

/**
 * External joint torque as a function of state, for qcollide_step_ext().
 *
 * @param q, qd     state at the stage being evaluated.
 * @param tau_ext   output, n elements.
 * @param n         number of joints.
 * @param user_data whatever was handed to qcollide_step_ext().
 */
typedef void (*qcollide_ext_torque_fn)(const double *q, const double *qd,
                                       double *tau_ext, size_t n,
                                       void *user_data);

/**
 * One fixed-step integration with a state-dependent external torque.
 *
 * qcollide_step() holds tau_ext constant across all four Runge-Kutta stages,
 * which is right for a commanded torque under zero-order hold but wrong for
 * anything that depends on the state. Friction is the case that matters here:
 * evaluating it once at the start of the step and freezing it reduces that
 * term to first order, so the integrator stops being fourth order exactly
 * when friction is switched on, which is the configuration the sweep spends
 * most of its time in.
 *
 * This variant evaluates the callback at each stage instead.
 * qcollide_step() is implemented in terms of it, with a callback that returns
 * the same vector every time, so there is one integrator rather than two.
 *
 * @param ext_fn may be NULL, meaning no external torque.
 */
qcollide_status qcollide_step_ext(const qcollide_model *m, double *q, double *qd,
                                  const double *tau,
                                  qcollide_ext_torque_fn ext_fn, void *user_data,
                                  double dt);

#ifdef __cplusplus
}
#endif

#endif /* QCOLLIDE_DYNAMICS_H */
