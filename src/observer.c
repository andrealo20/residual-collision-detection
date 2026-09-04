#include "qcollide/observer.h"

#include <math.h>

qcollide_status qcollide_observer_init(qcollide_observer *o,
                                       const qcollide_model *model,
                                       const double *gain)
{
    size_t i;

    if (o == NULL || model == NULL || gain == NULL) {
        return QCOLLIDE_ERR_ARG;
    }
    for (i = 0u; i < model->n; ++i) {
        /* A non-positive gain does not give a slow observer, it gives one
           that never converges or diverges. Rejecting here rather than
           producing a residual that drifts is the difference between an
           error and a wrong answer. */
        if (!(gain[i] > 0.0)) {
            return QCOLLIDE_ERR_ARG;
        }
    }

    o->model = model;
    /* Zero the whole arrays, not just the first n. qcollide_observer_residual
       returns the entire buffer, so a caller reading past n would otherwise
       see whatever was on the stack. */
    for (i = 0u; i < QCOLLIDE_MAX_DOF; ++i) {
        o->gain[i] = 0.0;
        o->integral[i] = 0.0;
        o->residual[i] = 0.0;
        o->p0[i] = 0.0;
        o->prev[i] = 0.0;
    }
    for (i = 0u; i < model->n; ++i) {
        o->gain[i] = gain[i];
        o->integral[i] = 0.0;
        o->residual[i] = 0.0;
        o->p0[i] = 0.0;
        o->prev[i] = 0.0;
    }
    o->started = 0;
    o->have_prev = 0;
    return QCOLLIDE_OK;
}

qcollide_status qcollide_observer_init_uniform(qcollide_observer *o,
                                               const qcollide_model *model,
                                               double gain)
{
    double g[QCOLLIDE_MAX_DOF];
    size_t i;

    if (model == NULL) {
        return QCOLLIDE_ERR_ARG;
    }
    for (i = 0u; i < model->n; ++i) {
        g[i] = gain;
    }
    return qcollide_observer_init(o, model, g);
}

qcollide_status qcollide_observer_reset(qcollide_observer *o,
                                        const double *q, const double *qd)
{
    size_t i;
    qcollide_status st;

    if (o == NULL || o->model == NULL || q == NULL || qd == NULL) {
        return QCOLLIDE_ERR_ARG;
    }

    st = qcollide_momentum(o->model, q, qd, o->p0);
    if (st != QCOLLIDE_OK) {
        return st;
    }
    for (i = 0u; i < o->model->n; ++i) {
        o->integral[i] = 0.0;
        o->residual[i] = 0.0;
        o->prev[i] = 0.0;
    }
    o->started = 1;
    o->have_prev = 0;
    return QCOLLIDE_OK;
}

qcollide_status qcollide_observer_update(qcollide_observer *o,
                                         const double *q, const double *qd,
                                         const double *tau, double dt)
{
    double p[QCOLLIDE_MAX_DOF];
    double C[QCOLLIDE_MAX_DOF * QCOLLIDE_MAX_DOF];
    double g[QCOLLIDE_MAX_DOF];
    size_t i, j, n;
    qcollide_status st;

    if (o == NULL || o->model == NULL || q == NULL || qd == NULL ||
        tau == NULL || !(dt > 0.0)) {
        return QCOLLIDE_ERR_ARG;
    }
    n = o->model->n;

    if (!o->started) {
        /* First sample: capture the momentum reference and report nothing.
           There is no history to integrate, and emitting a residual here
           would look exactly like a contact at start-up. */
        return qcollide_observer_reset(o, q, qd);
    }

    st = qcollide_momentum(o->model, q, qd, p);
    if (st != QCOLLIDE_OK) return st;
    st = qcollide_coriolis_matrix(o->model, q, qd, C);
    if (st != QCOLLIDE_OK) return st;
    st = qcollide_gravity(o->model, q, g);
    if (st != QCOLLIDE_OK) return st;

    for (i = 0u; i < n; ++i) {
        /* C' qd, the transpose being what the momentum identity calls for.
           Using C qd instead would give a residual that looks plausible and
           is wrong whenever the arm is moving. */
        double ctqd = 0.0;
        for (j = 0u; j < n; ++j) {
            ctqd += C[j * n + i] * qd[j];
        }

        /* The residual from the previous sample appears inside its own
           integral. That is what closes the loop and makes rdot = K(tau_ext
           - r) rather than an open integration that would drift on any model
           error at all. */
        const double integrand = tau[i] + ctqd - g[i] + o->residual[i];

        /* Trapezoidal rather than forward Euler. The integrand tracks pdot,
           which grows with joint velocity, so a first-order rule leaves an
           error proportional to dt times the speed of the arm -- measured at
           0.33 Nm of spurious residual on a joint moving at 19 rad/s, which
           is within the range of a real light contact. */
        if (o->have_prev) {
            o->integral[i] += 0.5 * (integrand + o->prev[i]) * dt;
        } else {
            o->integral[i] += integrand * dt;
        }
        o->prev[i] = integrand;

        o->residual[i] = o->gain[i] * (p[i] - o->integral[i] - o->p0[i]);
    }
    o->have_prev = 1;
    return QCOLLIDE_OK;
}

const double *qcollide_observer_residual(const qcollide_observer *o)
{
    return (o != NULL) ? o->residual : NULL;
}

double qcollide_observer_magnitude(const qcollide_observer *o)
{
    double worst = 0.0;
    size_t i;

    if (o == NULL || o->model == NULL) {
        return 0.0;
    }
    for (i = 0u; i < o->model->n; ++i) {
        const double a = fabs(o->residual[i]);
        if (a > worst) worst = a;
    }
    return worst;
}
