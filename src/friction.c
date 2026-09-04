#include "qcollide/friction.h"

#include <math.h>

void qcollide_friction_typical(qcollide_friction *f, double scale)
{
    if (f == NULL) return;
    /* Representative of a small collaborative-arm joint rather than
       identified from a particular robot: Coulomb of the order of a newton
       metre, stiction roughly half again as large, viscous of the order of a
       newton metre per rad/s. The README states this rather than implying
       the numbers came from hardware. */
    f->coulomb   = 0.8 * scale;
    f->stiction  = 1.2 * scale;
    f->viscous   = 0.5 * scale;
    f->stribeck  = 0.05;
    f->smoothing = 0.01;
}

double qcollide_friction_torque(const qcollide_friction *f, double qd)
{
    double s, dip, magnitude;

    if (f == NULL) return 0.0;

    /* tanh over a narrow band instead of a hard sign. Away from zero the two
       agree to machine precision; near zero the discontinuous version makes
       the integrator step across an instantaneous torque reversal. */
    s = tanh(qd / ((f->smoothing > 0.0) ? f->smoothing : 1e-6));

    /* Stribeck dip: friction falls from the static level towards Coulomb as
       the joint picks up speed. */
    dip = (f->stiction - f->coulomb)
        * exp(-fabs(qd / ((f->stribeck > 0.0) ? f->stribeck : 1e-6)));

    magnitude = f->coulomb + dip;

    /* Opposes motion, hence the negative sign; the viscous term likewise. */
    return -(magnitude * s + f->viscous * qd);
}

void qcollide_friction_vector(const qcollide_friction *f, size_t n,
                              const double *qd, double *tau_f)
{
    size_t i;
    if (qd == NULL || tau_f == NULL) return;
    if (n > QCOLLIDE_MAX_DOF) return;   /* the caller's n, as in sensors.c */
    for (i = 0u; i < n; ++i) {
        tau_f[i] = qcollide_friction_torque(f, qd[i]);
    }
}

qcollide_status qcollide_model_perturb(qcollide_model *dst,
                                       const qcollide_model *src,
                                       double error)
{
    size_t i;

    if (dst == NULL || src == NULL) {
        return QCOLLIDE_ERR_ARG;
    }
    /* Below -100 % the masses turn negative, which is not a wrong model but
       a nonsensical one. */
    if (error <= -1.0) {
        return QCOLLIDE_ERR_ARG;
    }

    *dst = *src;
    for (i = 0u; i < src->n; ++i) {
        dst->link[i].mass    = src->link[i].mass * (1.0 + error);
        dst->link[i].inertia = src->link[i].inertia * (1.0 + error);
    }
    return QCOLLIDE_OK;
}
