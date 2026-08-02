/**
 * @file detector.h
 * @brief Turning a residual into a collision decision, and measuring how well.
 *
 * The observer produces a continuous estimate of external joint torque.
 * Deciding whether a contact has occurred means comparing it against a
 * threshold, and the useful part of the problem is everything around that
 * comparison.
 *
 * ## Why a bare threshold is not enough
 *
 * A single sample above the threshold is not evidence of contact. Torque
 * measurement noise, quantisation, and any transient model error all produce
 * brief excursions, and a detector that fires on the first one produces
 * nuisance stops. Requiring the residual to stay above the threshold for a
 * number of consecutive samples trades detection latency for false alarms,
 * and the exchange rate is measurable rather than a matter of taste.
 *
 * Once triggered, the detector releases only when the residual falls below a
 * lower threshold. Without that hysteresis a residual sitting near the
 * threshold produces a stream of trigger-release events, which downstream is
 * worse than either decision taken consistently.
 *
 * ## What is measured
 *
 * Three quantities, all of them properties of a detector rather than of the
 * residual:
 *
 *   - **detection latency**: time from the contact starting to the detector
 *     firing. Bounded below by the observer time constant 1/K, since the
 *     residual cannot reach the threshold faster than it converges.
 *   - **false alarm rate**: firings per unit time with no contact present.
 *   - **minimum detectable torque**: the smallest external torque that fires
 *     the detector at all, found by bisection.
 *
 * These are the numbers a threshold has to be chosen against, and M3 measures
 * how unmodelled friction moves all three.
 */
#ifndef QCOLLIDE_DETECTOR_H
#define QCOLLIDE_DETECTOR_H

#include <stddef.h>

#include "qcollide/observer.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double threshold;        /**< Nm, residual magnitude to trigger */
    double release;          /**< Nm, magnitude to clear; must be < threshold */
    unsigned confirm;        /**< consecutive samples required to trigger */

    unsigned above;          /**< consecutive samples currently above */
    int      triggered;      /**< current decision */
    unsigned trigger_count;  /**< rising edges since reset */
} qcollide_detector;

/**
 * Configure a detector.
 *
 * @param threshold  strictly positive.
 * @param release    in [0, threshold). Equal to threshold means no
 *                   hysteresis, which is rejected: a detector that chatters
 *                   is worse than one biased in either direction.
 * @param confirm    consecutive samples required; at least 1.
 */
qcollide_status qcollide_detector_init(qcollide_detector *d, double threshold,
                                       double release, unsigned confirm);

/** Clear the decision and the sample counter, keeping the configuration. */
void qcollide_detector_reset(qcollide_detector *d);

/**
 * Feed one residual magnitude.
 * @return non-zero while the detector is triggered.
 */
int qcollide_detector_update(qcollide_detector *d, double magnitude);

/** Current decision. */
int qcollide_detector_triggered(const qcollide_detector *d);

/** Rising edges since the last reset. The false-alarm count when idle. */
unsigned qcollide_detector_events(const qcollide_detector *d);

/* ------------------------------------------------------------------ */
/* Contact model                                                       */
/* ------------------------------------------------------------------ */

/**
 * A contact, expressed as an external joint torque over a time window.
 *
 * Real contacts are forces at a point on a link, which map to joint torques
 * through the transpose of a Jacobian. Specifying the joint torque directly
 * skips that step deliberately: it makes the ground truth exact, so a
 * detection result measures the detector rather than the contact model. A
 * force-and-location interface would add a second thing that could be wrong.
 */
typedef struct {
    double start;                     /**< s */
    double duration;                  /**< s; 0 or less means until the end */
    double torque[QCOLLIDE_MAX_DOF];  /**< Nm applied at the joints */
} qcollide_contact;

/** External torque at time t, zero outside the contact window. */
void qcollide_contact_torque(const qcollide_contact *c, size_t n,
                             double t, double *tau_ext);

#ifdef __cplusplus
}
#endif

#endif /* QCOLLIDE_DETECTOR_H */
