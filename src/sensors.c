#include "qcollide/sensors.h"

#include <math.h>

/* ------------------------------------------------------------------ */
/* Deterministic noise                                                 */
/* ------------------------------------------------------------------ */

void qcollide_noise_seed(qcollide_noise *n, uint64_t seed)
{
    if (n == NULL) return;
    /* Any non-zero state works; xorshift is degenerate at zero. */
    n->state = (seed != 0u) ? seed : 0x9E3779B97F4A7C15ull;
    n->spare = 0.0;
    n->have_spare = 0;
}

static double uniform01(qcollide_noise *n)
{
    /* xorshift64*, adequate for injecting measurement noise and small enough
       to be obviously deterministic across platforms, which rand() is not. */
    uint64_t x = n->state;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    n->state = x;
    x *= 0x2545F4914F6CDD1Dull;
    /* Top 53 bits into [0,1), avoiding exactly zero for the logarithm. */
    return ((double)(x >> 11) + 0.5) * (1.0 / 9007199254740992.0);
}

double qcollide_noise_gaussian(qcollide_noise *n)
{
    double u1, u2, r, theta;

    if (n == NULL) return 0.0;
    if (n->have_spare) {
        n->have_spare = 0;
        return n->spare;
    }

    u1 = uniform01(n);
    u2 = uniform01(n);
    r = sqrt(-2.0 * log(u1));
    theta = 6.283185307179586 * u2;

    n->spare = r * sin(theta);
    n->have_spare = 1;
    return r * cos(theta);
}

/* ------------------------------------------------------------------ */
/* Sensors                                                             */
/* ------------------------------------------------------------------ */

void qcollide_sensors_typical(qcollide_sensor_config *cfg, double torque_scale)
{
    if (cfg == NULL) return;
    /* 17-bit encoder and torque noise of the order of one percent of a small
       joint's rated torque. Representative of an industrial joint rather
       than measured from one, and the README says so. */
    cfg->encoder_counts = 131072u;
    cfg->torque_noise = 0.05 * torque_scale;
    cfg->velocity_mode = QCOLLIDE_VEL_FILTERED;
    cfg->velocity_cutoff = 50.0;
}

qcollide_status qcollide_sensors_init(qcollide_sensors *s,
                                      const qcollide_sensor_config *cfg,
                                      uint64_t seed)
{
    if (s == NULL || cfg == NULL) {
        return QCOLLIDE_ERR_ARG;
    }
    if (cfg->torque_noise < 0.0) {
        return QCOLLIDE_ERR_ARG;
    }
    if (cfg->velocity_mode == QCOLLIDE_VEL_FILTERED &&
        !(cfg->velocity_cutoff > 0.0)) {
        return QCOLLIDE_ERR_ARG;
    }

    s->cfg = *cfg;
    qcollide_noise_seed(&s->noise, seed);
    qcollide_sensors_reset(s);
    return QCOLLIDE_OK;
}

void qcollide_sensors_reset(qcollide_sensors *s)
{
    size_t i;
    if (s == NULL) return;
    for (i = 0u; i < QCOLLIDE_MAX_DOF; ++i) {
        s->last_q[i] = 0.0;
        s->vel[i] = 0.0;
    }
    s->have_last = 0;
}

/** Round to the nearest encoder count. */
static double quantise(double angle, uint32_t counts)
{
    double step;
    if (counts == 0u) return angle;
    step = 6.283185307179586 / (double)counts;
    return floor(angle / step + 0.5) * step;
}

qcollide_status qcollide_sensors_measure(qcollide_sensors *s, size_t n,
                                         const double *q_true,
                                         const double *qd_true,
                                         const double *tau_true,
                                         double *q_meas, double *qd_meas,
                                         double *tau_meas, double dt)
{
    size_t i;
    double alpha = 0.0;

    if (s == NULL || q_true == NULL || qd_true == NULL || tau_true == NULL ||
        q_meas == NULL || qd_meas == NULL || tau_meas == NULL || !(dt > 0.0)) {
        return QCOLLIDE_ERR_ARG;
    }
    /* n is the caller's, unlike everywhere else in the library where it comes
       from a validated model, and it indexes s->last_q and s->vel. */
    if (n > QCOLLIDE_MAX_DOF) {
        return QCOLLIDE_ERR_ARG;
    }

    if (s->cfg.velocity_mode == QCOLLIDE_VEL_FILTERED) {
        /* First-order low pass, discretised as a one-pole recursion. The
           filter is what trades velocity noise for lag, and the lag appears
           in the residual as a torque that is not there -- which is exactly
           what M4 measures. */
        const double tau_f = 1.0 / (6.283185307179586 * s->cfg.velocity_cutoff);
        alpha = dt / (tau_f + dt);
    }

    for (i = 0u; i < n; ++i) {
        const double q_q = quantise(q_true[i], s->cfg.encoder_counts);
        q_meas[i] = q_q;

        switch (s->cfg.velocity_mode) {
        case QCOLLIDE_VEL_EXACT:
            s->vel[i] = qd_true[i];
            break;

        case QCOLLIDE_VEL_DIFF:
        case QCOLLIDE_VEL_FILTERED:
        default:
            if (!s->have_last) {
                /* No history: report the true velocity for one sample
                   rather than a difference against an undefined previous
                   position, which would be a large spurious transient at
                   exactly the moment the observer initialises. */
                s->vel[i] = qd_true[i];
            } else {
                const double raw = (q_q - s->last_q[i]) / dt;
                if (s->cfg.velocity_mode == QCOLLIDE_VEL_DIFF) {
                    s->vel[i] = raw;
                } else {
                    s->vel[i] += alpha * (raw - s->vel[i]);
                }
            }
            break;
        }

        s->last_q[i] = q_q;
        qd_meas[i] = s->vel[i];

        tau_meas[i] = tau_true[i] +
                      s->cfg.torque_noise * qcollide_noise_gaussian(&s->noise);
    }

    s->have_last = 1;
    return QCOLLIDE_OK;
}
