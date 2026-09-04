#include "qcollide/dynamics.h"

#include <math.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* Model setup                                                         */
/* ------------------------------------------------------------------ */

qcollide_status qcollide_model_init(qcollide_model *m,
                                    const qcollide_link *links,
                                    size_t n, double gravity)
{
    size_t i;

    if (m == NULL || links == NULL || n == 0u || n > QCOLLIDE_MAX_DOF) {
        return QCOLLIDE_ERR_ARG;
    }
    for (i = 0u; i < n; ++i) {
        /* A non-positive mass or length is not an unusual robot, it is a
           typo. Catching it here avoids a mass matrix that is silently
           indefinite much later. */
        if (links[i].mass <= 0.0 || links[i].length <= 0.0 ||
            links[i].inertia < 0.0) {
            return QCOLLIDE_ERR_ARG;
        }
        /* A centre of mass outside the link is the same class of mistake as a
           negative mass: not an unusual robot, a typo. Caught here rather than
           as a mass matrix that is quietly wrong. */
        if (links[i].com < 0.0 || links[i].com > links[i].length) {
            return QCOLLIDE_ERR_ARG;
        }
        m->link[i] = links[i];
    }
    m->n = n;
    m->gravity = gravity;
    return QCOLLIDE_OK;
}

qcollide_status qcollide_model_uniform(qcollide_model *m, size_t n,
                                       double mass, double length,
                                       double gravity)
{
    qcollide_link links[QCOLLIDE_MAX_DOF];
    size_t i;

    if (n == 0u || n > QCOLLIDE_MAX_DOF) {
        return QCOLLIDE_ERR_ARG;
    }
    for (i = 0u; i < n; ++i) {
        links[i].mass = mass;
        links[i].length = length;
        links[i].com = 0.5 * length;
        links[i].inertia = mass * length * length / 12.0;
    }
    return qcollide_model_init(m, links, n, gravity);
}

/* ------------------------------------------------------------------ */
/* Kinematic helpers                                                   */
/* ------------------------------------------------------------------ */

/**
 * Absolute link angles: for a planar serial chain the orientation of link i
 * is the sum of joints 0..i. Computing these once and reusing them is what
 * keeps the mass matrix readable.
 */
static void absolute_angles(const qcollide_model *m, const double *q,
                            double *theta)
{
    size_t i;
    double acc = 0.0;
    for (i = 0u; i < m->n; ++i) {
        acc += q[i];
        theta[i] = acc;
    }
}

/**
 * Position of the centre of mass of link k, and of every joint before it.
 *
 * The mass matrix entry M(i,j) is a sum over links k >= max(i,j) of terms
 * built from the partial derivatives of the centre-of-mass position with
 * respect to q_i and q_j. Those partials are computed here.
 *
 * dp_k/dq_i for a planar chain is the vector from joint i to the centre of
 * mass of link k, rotated by 90 degrees, and it is zero for i > k. That
 * geometric fact is the whole derivation.
 */
static void com_jacobian(const qcollide_model *m, const double *theta,
                         size_t k, double *jx, double *jy)
{
    size_t i, j;

    for (i = 0u; i < m->n; ++i) {
        jx[i] = 0.0;
        jy[i] = 0.0;
    }

    for (i = 0u; i <= k; ++i) {
        /* Vector from joint i to the centre of mass of link k. */
        double dx = 0.0, dy = 0.0;
        for (j = i; j < k; ++j) {
            dx += m->link[j].length * cos(theta[j]);
            dy += m->link[j].length * sin(theta[j]);
        }
        dx += m->link[k].com * cos(theta[k]);
        dy += m->link[k].com * sin(theta[k]);

        /* A rotation of the joint turns that vector: the derivative is the
           vector rotated by +90 degrees. */
        jx[i] = -dy;
        jy[i] =  dx;
    }
}

/* ------------------------------------------------------------------ */
/* Mass matrix                                                         */
/* ------------------------------------------------------------------ */

qcollide_status qcollide_mass_matrix(const qcollide_model *m,
                                     const double *q, double *M)
{
    double theta[QCOLLIDE_MAX_DOF];
    double jx[QCOLLIDE_MAX_DOF], jy[QCOLLIDE_MAX_DOF];
    size_t i, j, k, n;

    if (m == NULL || q == NULL || M == NULL) {
        return QCOLLIDE_ERR_ARG;
    }
    n = m->n;
    absolute_angles(m, q, theta);

    for (i = 0u; i < n * n; ++i) {
        M[i] = 0.0;
    }

    for (k = 0u; k < n; ++k) {
        com_jacobian(m, theta, k, jx, jy);
        for (i = 0u; i < n; ++i) {
            for (j = 0u; j < n; ++j) {
                /* Translational part: m_k * Jv_i . Jv_j
                   Rotational part: I_k, present whenever both joints are at
                   or before link k, since every one of them rotates it. */
                double t = m->link[k].mass * (jx[i] * jx[j] + jy[i] * jy[j]);
                double r = (i <= k && j <= k) ? m->link[k].inertia : 0.0;
                M[i * n + j] += t + r;
            }
        }
    }
    return QCOLLIDE_OK;
}

/* ------------------------------------------------------------------ */
/* Coriolis matrix                                                     */
/* ------------------------------------------------------------------ */

/**
 * Partial derivatives of the mass matrix, dM/dq_k for every k, by central
 * differences.
 *
 * Two design points worth recording.
 *
 * First, why finite differences at all: the analytical derivative of M is
 * available but long, and an error in it would break the skew-symmetry of
 * Mdot - 2C without breaking anything visible. The observer built on top
 * would still run and still produce plausible residuals. Central differences
 * are exact to second order and the skew-symmetry is verified numerically in
 * the tests, so the cost buys a derivation that cannot be quietly wrong.
 *
 * Second, why all n derivatives are computed at once. The Christoffel sum
 * needs dM_ij/dq_k, dM_ik/dq_j and dM_jk/dq_i for every triple, and a naive
 * implementation calls a derivative helper inside the triple loop. That
 * repeats the same handful of distinct derivatives O(n^3) times: for a
 * 3-link arm it costs 162 mass matrix evaluations per Coriolis call, and
 * RK4 needs four of those per integration step. There are only n distinct
 * derivatives, so computing them once costs 2n evaluations instead --
 * twenty-seven times less work at n = 3, and the gap grows as n^2. The
 * first version of this file did it the naive way and the sanitizer build
 * of the test suite took long enough to time out, which is how it was
 * noticed.
 */
static void mass_derivatives(const qcollide_model *m, const double *q,
                             double *dMdq)
{
    double qp[QCOLLIDE_MAX_DOF], qm[QCOLLIDE_MAX_DOF];
    double Mp[QCOLLIDE_MAX_DOF * QCOLLIDE_MAX_DOF];
    double Mm[QCOLLIDE_MAX_DOF * QCOLLIDE_MAX_DOF];
    const double h = 1e-6;
    size_t a, k, idx, n = m->n;

    for (k = 0u; k < n; ++k) {
        for (a = 0u; a < n; ++a) {
            qp[a] = q[a];
            qm[a] = q[a];
        }
        qp[k] += h;
        qm[k] -= h;

        (void)qcollide_mass_matrix(m, qp, Mp);
        (void)qcollide_mass_matrix(m, qm, Mm);

        for (idx = 0u; idx < n * n; ++idx) {
            dMdq[k * n * n + idx] = (Mp[idx] - Mm[idx]) / (2.0 * h);
        }
    }
}

qcollide_status qcollide_coriolis_matrix(const qcollide_model *m,
                                         const double *q, const double *qd,
                                         double *C)
{
    double dMdq[QCOLLIDE_MAX_DOF * QCOLLIDE_MAX_DOF * QCOLLIDE_MAX_DOF];
    size_t i, j, k, n;

    if (m == NULL || q == NULL || qd == NULL || C == NULL) {
        return QCOLLIDE_ERR_ARG;
    }
    n = m->n;

    mass_derivatives(m, q, dMdq);

    for (i = 0u; i < n; ++i) {
        for (j = 0u; j < n; ++j) {
            double c = 0.0;
            for (k = 0u; k < n; ++k) {
                /* Christoffel form:
                   C(i,j) = 1/2 sum_k ( dM_ij/dq_k + dM_ik/dq_j - dM_jk/dq_i ) qd_k
                   This is the form for which Mdot - 2C is skew-symmetric. */
                const double t = dMdq[k * n * n + i * n + j]
                               + dMdq[j * n * n + i * n + k]
                               - dMdq[i * n * n + j * n + k];
                c += 0.5 * t * qd[k];
            }
            C[i * n + j] = c;
        }
    }
    return QCOLLIDE_OK;
}

/* ------------------------------------------------------------------ */
/* Gravity                                                             */
/* ------------------------------------------------------------------ */

qcollide_status qcollide_gravity(const qcollide_model *m,
                                 const double *q, double *g)
{
    double theta[QCOLLIDE_MAX_DOF];
    double jx[QCOLLIDE_MAX_DOF], jy[QCOLLIDE_MAX_DOF];
    size_t i, k, n;

    if (m == NULL || q == NULL || g == NULL) {
        return QCOLLIDE_ERR_ARG;
    }
    n = m->n;
    absolute_angles(m, q, theta);

    for (i = 0u; i < n; ++i) {
        g[i] = 0.0;
    }

    /* Gravity acts along -y, so the potential is sum_k m_k * gravity * y_k
       and the joint torque is its gradient. */
    for (k = 0u; k < n; ++k) {
        com_jacobian(m, theta, k, jx, jy);
        for (i = 0u; i < n; ++i) {
            g[i] += m->link[k].mass * m->gravity * jy[i];
        }
    }
    return QCOLLIDE_OK;
}

/* ------------------------------------------------------------------ */
/* Momentum and energy                                                 */
/* ------------------------------------------------------------------ */

qcollide_status qcollide_momentum(const qcollide_model *m,
                                  const double *q, const double *qd,
                                  double *p)
{
    double M[QCOLLIDE_MAX_DOF * QCOLLIDE_MAX_DOF];
    size_t i, j, n;
    qcollide_status st;

    if (m == NULL || q == NULL || qd == NULL || p == NULL) {
        return QCOLLIDE_ERR_ARG;
    }
    n = m->n;
    st = qcollide_mass_matrix(m, q, M);
    if (st != QCOLLIDE_OK) return st;

    for (i = 0u; i < n; ++i) {
        double s = 0.0;
        for (j = 0u; j < n; ++j) {
            s += M[i * n + j] * qd[j];
        }
        p[i] = s;
    }
    return QCOLLIDE_OK;
}

double qcollide_energy(const qcollide_model *m,
                       const double *q, const double *qd)
{
    double M[QCOLLIDE_MAX_DOF * QCOLLIDE_MAX_DOF];
    double theta[QCOLLIDE_MAX_DOF];
    double kinetic = 0.0, potential = 0.0;
    size_t i, j, k, n;

    if (m == NULL || q == NULL || qd == NULL) {
        return 0.0;
    }
    n = m->n;
    (void)qcollide_mass_matrix(m, q, M);

    for (i = 0u; i < n; ++i) {
        for (j = 0u; j < n; ++j) {
            kinetic += 0.5 * qd[i] * M[i * n + j] * qd[j];
        }
    }

    absolute_angles(m, q, theta);
    for (k = 0u; k < n; ++k) {
        double y = 0.0;
        for (j = 0u; j < k; ++j) {
            y += m->link[j].length * sin(theta[j]);
        }
        y += m->link[k].com * sin(theta[k]);
        potential += m->link[k].mass * m->gravity * y;
    }
    return kinetic + potential;
}

/* ------------------------------------------------------------------ */
/* Forward and inverse dynamics                                        */
/* ------------------------------------------------------------------ */

/**
 * Solve A x = b in place by Gauss-Jordan with partial pivoting.
 *
 * Named for what it does. It was called solve_spd, which promised a Cholesky
 * factorisation exploiting the symmetry and positive definiteness of the mass
 * matrix; this is the general method, and the name should not claim otherwise.
 */
static int solve_linear_system(double *A, double *b, size_t n)
{
    size_t i, j, k;
    double scale = 0.0;
    double pivot_floor;

    /* The singularity test has to be relative to the size of the entries. An
       absolute floor treats a model expressed in grams and millimetres, whose
       mass matrix is legitimately of order 1e-9, as singular. */
    for (i = 0u; i < n * n; ++i) {
        const double v = fabs(A[i]);
        if (v > scale) { scale = v; }
    }
    if (!(scale > 0.0)) {
        return 0;
    }
    pivot_floor = 1e-12 * scale;

    for (k = 0u; k < n; ++k) {
        size_t piv = k;
        double best = fabs(A[k * n + k]);
        for (i = k + 1u; i < n; ++i) {
            const double v = fabs(A[i * n + k]);
            if (v > best) { best = v; piv = i; }
        }
        /* The mass matrix of a physical arm is positive definite, so a tiny
           pivot means the parameters are wrong, not that the configuration
           is awkward. Reporting rather than continuing keeps that
           distinction visible. */
        if (best < pivot_floor) {
            return 0;
        }
        if (piv != k) {
            for (j = 0u; j < n; ++j) {
                const double t = A[k * n + j];
                A[k * n + j] = A[piv * n + j];
                A[piv * n + j] = t;
            }
            { const double t = b[k]; b[k] = b[piv]; b[piv] = t; }
        }
        {
            const double d = A[k * n + k];
            for (j = k; j < n; ++j) A[k * n + j] /= d;
            b[k] /= d;
        }
        for (i = 0u; i < n; ++i) {
            if (i == k) continue;
            {
                const double f = A[i * n + k];
                if (f == 0.0) continue;
                for (j = k; j < n; ++j) A[i * n + j] -= f * A[k * n + j];
                b[i] -= f * b[k];
            }
        }
    }
    return 1;
}

qcollide_status qcollide_forward_dynamics(const qcollide_model *m,
                                          const double *q, const double *qd,
                                          const double *tau,
                                          const double *tau_ext,
                                          double *qdd)
{
    double M[QCOLLIDE_MAX_DOF * QCOLLIDE_MAX_DOF];
    double C[QCOLLIDE_MAX_DOF * QCOLLIDE_MAX_DOF];
    double g[QCOLLIDE_MAX_DOF];
    double rhs[QCOLLIDE_MAX_DOF];
    size_t i, j, n;

    if (m == NULL || q == NULL || qd == NULL || tau == NULL || qdd == NULL) {
        return QCOLLIDE_ERR_ARG;
    }
    n = m->n;

    (void)qcollide_mass_matrix(m, q, M);
    (void)qcollide_coriolis_matrix(m, q, qd, C);
    (void)qcollide_gravity(m, q, g);

    for (i = 0u; i < n; ++i) {
        double cq = 0.0;
        for (j = 0u; j < n; ++j) {
            cq += C[i * n + j] * qd[j];
        }
        rhs[i] = tau[i] - cq - g[i] + ((tau_ext != NULL) ? tau_ext[i] : 0.0);
    }

    if (!solve_linear_system(M, rhs, n)) {
        return QCOLLIDE_ERR_RANGE;
    }
    for (i = 0u; i < n; ++i) {
        qdd[i] = rhs[i];
    }
    return QCOLLIDE_OK;
}

qcollide_status qcollide_inverse_dynamics(const qcollide_model *m,
                                          const double *q, const double *qd,
                                          const double *qdd, double *tau)
{
    double M[QCOLLIDE_MAX_DOF * QCOLLIDE_MAX_DOF];
    double C[QCOLLIDE_MAX_DOF * QCOLLIDE_MAX_DOF];
    double g[QCOLLIDE_MAX_DOF];
    size_t i, j, n;

    if (m == NULL || q == NULL || qd == NULL || qdd == NULL || tau == NULL) {
        return QCOLLIDE_ERR_ARG;
    }
    n = m->n;

    (void)qcollide_mass_matrix(m, q, M);
    (void)qcollide_coriolis_matrix(m, q, qd, C);
    (void)qcollide_gravity(m, q, g);

    for (i = 0u; i < n; ++i) {
        double s = g[i];
        for (j = 0u; j < n; ++j) {
            s += M[i * n + j] * qdd[j] + C[i * n + j] * qd[j];
        }
        tau[i] = s;
    }
    return QCOLLIDE_OK;
}

/* ------------------------------------------------------------------ */
/* Integration                                                         */
/* ------------------------------------------------------------------ */

/**
 * External torque held constant across the step, which is what the plain
 * qcollide_step() promises. Used to express it in terms of the callback form
 * so there is one integrator rather than two.
 */
typedef struct {
    const double *tau_ext;
} qcollide_const_ext;

static void qcollide_const_ext_fn(const double *q, const double *qd,
                                  double *tau_ext, size_t n, void *user_data)
{
    const qcollide_const_ext *c = (const qcollide_const_ext *)user_data;
    size_t i;
    (void)q;
    (void)qd;
    for (i = 0u; i < n; ++i) {
        tau_ext[i] = (c->tau_ext != NULL) ? c->tau_ext[i] : 0.0;
    }
}

qcollide_status qcollide_step_ext(const qcollide_model *m, double *q, double *qd,
                                  const double *tau,
                                  qcollide_ext_torque_fn ext_fn, void *user_data,
                                  double dt)
{
    /* Zero-initialised because the compiler cannot prove n > 0 from the model
       invariant, and warns at -O2. The invariant does hold -- model init
       rejects n == 0 -- so this costs a few stores and removes a warning that
       would otherwise have to be suppressed. */
    double k1q[QCOLLIDE_MAX_DOF] = { 0 }, k1v[QCOLLIDE_MAX_DOF] = { 0 };
    double k2q[QCOLLIDE_MAX_DOF] = { 0 }, k2v[QCOLLIDE_MAX_DOF] = { 0 };
    double k3q[QCOLLIDE_MAX_DOF] = { 0 }, k3v[QCOLLIDE_MAX_DOF] = { 0 };
    double k4q[QCOLLIDE_MAX_DOF] = { 0 }, k4v[QCOLLIDE_MAX_DOF] = { 0 };
    double tq[QCOLLIDE_MAX_DOF] = { 0 }, tv[QCOLLIDE_MAX_DOF] = { 0 };
    double ext[QCOLLIDE_MAX_DOF] = { 0 };
    size_t i, n;
    qcollide_status st;

    if (m == NULL || q == NULL || qd == NULL || tau == NULL) {
        return QCOLLIDE_ERR_ARG;
    }
    n = m->n;

    /* Stage 1 */
    if (ext_fn != NULL) { ext_fn(q, qd, ext, n, user_data); }
    st = qcollide_forward_dynamics(m, q, qd, tau, ext, k1v);
    if (st != QCOLLIDE_OK) return st;
    for (i = 0u; i < n; ++i) k1q[i] = qd[i];

    /* Stage 2 */
    for (i = 0u; i < n; ++i) {
        tq[i] = q[i] + 0.5 * dt * k1q[i];
        tv[i] = qd[i] + 0.5 * dt * k1v[i];
    }
    if (ext_fn != NULL) { ext_fn(tq, tv, ext, n, user_data); }
    st = qcollide_forward_dynamics(m, tq, tv, tau, ext, k2v);
    if (st != QCOLLIDE_OK) return st;
    for (i = 0u; i < n; ++i) k2q[i] = tv[i];

    /* Stage 3 */
    for (i = 0u; i < n; ++i) {
        tq[i] = q[i] + 0.5 * dt * k2q[i];
        tv[i] = qd[i] + 0.5 * dt * k2v[i];
    }
    if (ext_fn != NULL) { ext_fn(tq, tv, ext, n, user_data); }
    st = qcollide_forward_dynamics(m, tq, tv, tau, ext, k3v);
    if (st != QCOLLIDE_OK) return st;
    for (i = 0u; i < n; ++i) k3q[i] = tv[i];

    /* Stage 4 */
    for (i = 0u; i < n; ++i) {
        tq[i] = q[i] + dt * k3q[i];
        tv[i] = qd[i] + dt * k3v[i];
    }
    if (ext_fn != NULL) { ext_fn(tq, tv, ext, n, user_data); }
    st = qcollide_forward_dynamics(m, tq, tv, tau, ext, k4v);
    if (st != QCOLLIDE_OK) return st;
    for (i = 0u; i < n; ++i) k4q[i] = tv[i];

    for (i = 0u; i < n; ++i) {
        q[i]  += dt / 6.0 * (k1q[i] + 2.0 * k2q[i] + 2.0 * k3q[i] + k4q[i]);
        qd[i] += dt / 6.0 * (k1v[i] + 2.0 * k2v[i] + 2.0 * k3v[i] + k4v[i]);
    }
    return QCOLLIDE_OK;
}

qcollide_status qcollide_step(const qcollide_model *m, double *q, double *qd,
                              const double *tau, const double *tau_ext,
                              double dt)
{
    qcollide_const_ext holder;
    holder.tau_ext = tau_ext;
    return qcollide_step_ext(m, q, qd, tau, qcollide_const_ext_fn, &holder, dt);
}
