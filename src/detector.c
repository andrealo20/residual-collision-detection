#include "qcollide/detector.h"

qcollide_status qcollide_detector_init(qcollide_detector *d, double threshold,
                                       double release, unsigned confirm)
{
    if (d == NULL || !(threshold > 0.0) || confirm == 0u) {
        return QCOLLIDE_ERR_ARG;
    }
    /* release == threshold is rejected rather than allowed as "no
       hysteresis": a residual hovering at the threshold would then produce a
       stream of trigger and release events, which downstream is worse than
       either decision held consistently. */
    if (release < 0.0 || release >= threshold) {
        return QCOLLIDE_ERR_ARG;
    }

    d->threshold = threshold;
    d->release = release;
    d->confirm = confirm;
    qcollide_detector_reset(d);
    return QCOLLIDE_OK;
}

void qcollide_detector_reset(qcollide_detector *d)
{
    if (d == NULL) return;
    d->above = 0u;
    d->triggered = 0;
    d->trigger_count = 0u;
}

int qcollide_detector_update(qcollide_detector *d, double magnitude)
{
    if (d == NULL) return 0;

    if (d->triggered) {
        /* Release on the lower threshold, immediately: holding a stale
           trigger after the contact has clearly gone adds latency to the
           recovery without adding robustness. */
        if (magnitude < d->release) {
            d->triggered = 0;
            d->above = 0u;
        }
        return d->triggered;
    }

    if (magnitude >= d->threshold) {
        d->above++;
        if (d->above >= d->confirm) {
            d->triggered = 1;
            d->trigger_count++;
        }
    } else {
        /* The counter resets on a single sample below. Requiring N
           consecutive samples is the point; allowing gaps would let noise
           accumulate towards a trigger over an arbitrarily long window. */
        d->above = 0u;
    }
    return d->triggered;
}

int qcollide_detector_triggered(const qcollide_detector *d)
{
    return (d != NULL) ? d->triggered : 0;
}

unsigned qcollide_detector_events(const qcollide_detector *d)
{
    return (d != NULL) ? d->trigger_count : 0u;
}

void qcollide_contact_torque(const qcollide_contact *c, size_t n,
                             double t, double *tau_ext)
{
    size_t i;
    int active;

    if (tau_ext == NULL) return;
    for (i = 0u; i < n; ++i) {
        tau_ext[i] = 0.0;
    }
    if (c == NULL) return;

    active = (t >= c->start) &&
             ((c->duration <= 0.0) || (t < c->start + c->duration));
    if (!active) return;

    for (i = 0u; i < n; ++i) {
        tau_ext[i] = c->torque[i];
    }
}
