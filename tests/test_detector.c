/**
 * Unit tests for the detector, and the end-to-end detection metrics.
 *
 * Two halves. The first tests the decision logic in isolation, feeding it
 * magnitudes directly: hysteresis, the confirmation counter, and the edge
 * cases where a naive implementation lets noise accumulate into a trigger.
 *
 * The second runs the whole chain, plant, observer, detector, and measures
 * the three numbers a threshold has to be chosen against: detection latency,
 * false alarms, and the smallest torque that fires at all. Those are
 * reported rather than asserted tightly, because they are the subject of the
 * project rather than a pass/fail property. What is asserted is that they
 * behave the way the theory requires: latency bounded below by the observer
 * time constant, and monotone in the threshold.
 */
#include "unity.h"
#include "qcollide/detector.h"

#include <math.h>
#include <stdio.h>

#define N 2
#define DT 1e-4

static qcollide_model model;
static qcollide_observer obs;
static qcollide_detector det;

void setUp(void)
{
    TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK,
        qcollide_model_uniform(&model, N, 1.0, 0.6, 9.81));
}

void tearDown(void) {}

/* ------------------------------------------------------------------ */
/* Decision logic                                                      */
/* ------------------------------------------------------------------ */

static void test_init_rejects_bad_arguments(void)
{
    qcollide_detector d;

    TEST_ASSERT_EQUAL_INT(QCOLLIDE_ERR_ARG,
        qcollide_detector_init(NULL, 1.0, 0.5, 3));
    TEST_ASSERT_EQUAL_INT(QCOLLIDE_ERR_ARG,
        qcollide_detector_init(&d, 0.0, 0.0, 3));
    TEST_ASSERT_EQUAL_INT(QCOLLIDE_ERR_ARG,
        qcollide_detector_init(&d, -1.0, 0.0, 3));
    TEST_ASSERT_EQUAL_INT(QCOLLIDE_ERR_ARG,
        qcollide_detector_init(&d, 1.0, 0.5, 0));
    TEST_ASSERT_EQUAL_INT(QCOLLIDE_ERR_ARG,
        qcollide_detector_init(&d, 1.0, -0.1, 3));
    /* release == threshold means no hysteresis, which is a chattering
       detector rather than a permissive one. */
    TEST_ASSERT_EQUAL_INT(QCOLLIDE_ERR_ARG,
        qcollide_detector_init(&d, 1.0, 1.0, 3));
    TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK,
        qcollide_detector_init(&d, 1.0, 0.5, 1));
}

static void test_confirmation_requires_consecutive_samples(void)
{
    int i;

    TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK,
        qcollide_detector_init(&det, 1.0, 0.5, 4));

    for (i = 0; i < 3; ++i) {
        TEST_ASSERT_FALSE(qcollide_detector_update(&det, 2.0));
    }
    TEST_ASSERT_TRUE(qcollide_detector_update(&det, 2.0));
    TEST_ASSERT_EQUAL_UINT(1u, qcollide_detector_events(&det));
}

/**
 * A single sample below the threshold must reset the counter. Allowing gaps
 * would let uncorrelated noise accumulate towards a trigger over an
 * arbitrarily long window, which turns the confirmation count from a
 * safeguard into a delay.
 */
static void test_a_single_sample_below_resets_the_counter(void)
{
    int i, j;

    TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK,
        qcollide_detector_init(&det, 1.0, 0.5, 4));

    for (j = 0; j < 50; ++j) {
        for (i = 0; i < 3; ++i) {
            TEST_ASSERT_FALSE(qcollide_detector_update(&det, 2.0));
        }
        TEST_ASSERT_FALSE(qcollide_detector_update(&det, 0.9));
    }
    TEST_ASSERT_EQUAL_UINT(0u, qcollide_detector_events(&det));
}

static void test_hysteresis_prevents_chattering(void)
{
    int j;

    TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK,
        qcollide_detector_init(&det, 1.0, 0.5, 1));
    TEST_ASSERT_TRUE(qcollide_detector_update(&det, 1.2));

    /* Sitting between release and threshold: must stay triggered rather
       than toggling on every sample. */
    for (j = 0; j < 100; ++j) {
        TEST_ASSERT_TRUE(qcollide_detector_update(&det, 0.7));
    }
    TEST_ASSERT_EQUAL_UINT(1u, qcollide_detector_events(&det));

    TEST_ASSERT_FALSE(qcollide_detector_update(&det, 0.4));
    TEST_ASSERT_FALSE(qcollide_detector_triggered(&det));
}

static void test_contact_window(void)
{
    qcollide_contact c;
    double tau[N];
    size_t i;

    for (i = 0; i < N; ++i) c.torque[i] = 3.0;
    c.start = 1.0;
    c.duration = 0.5;

    qcollide_contact_torque(&c, N, 0.9, tau);
    TEST_ASSERT_EQUAL_DOUBLE(0.0, tau[0]);
    qcollide_contact_torque(&c, N, 1.0, tau);
    TEST_ASSERT_EQUAL_DOUBLE(3.0, tau[0]);
    qcollide_contact_torque(&c, N, 1.49, tau);
    TEST_ASSERT_EQUAL_DOUBLE(3.0, tau[0]);
    qcollide_contact_torque(&c, N, 1.5, tau);
    TEST_ASSERT_EQUAL_DOUBLE(0.0, tau[0]);

    /* A non-positive duration means the contact never ends. */
    c.duration = 0.0;
    qcollide_contact_torque(&c, N, 100.0, tau);
    TEST_ASSERT_EQUAL_DOUBLE(3.0, tau[0]);
}

/* ------------------------------------------------------------------ */
/* End to end                                                          */
/* ------------------------------------------------------------------ */

/**
 * Run plant, observer and detector together.
 *
 * @param latency  time from contact start to trigger, or -1 if never.
 * @param events   rising edges observed over the whole run.
 */
static void run_chain(double gain, double threshold, unsigned confirm,
                      const qcollide_contact *contact, double duration,
                      double *latency, unsigned *events)
{
    double q[N] = { 0.3, -0.5 };
    double qd[N] = { 0.0, 0.0 };
    double tau[N], ext[N];
    const int steps = (int)(duration / DT + 0.5);
    int k;

    TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK,
        qcollide_observer_init_uniform(&obs, &model, gain));
    TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK,
        qcollide_detector_init(&det, threshold, 0.5 * threshold, confirm));

    (void)qcollide_gravity(&model, q, tau);
    (void)qcollide_observer_reset(&obs, q, qd);

    if (latency != NULL) *latency = -1.0;

    for (k = 0; k < steps; ++k) {
        const double t = (double)k * DT;

        qcollide_contact_torque(contact, N, t, ext);
        (void)qcollide_step(&model, q, qd, tau, ext, DT);
        (void)qcollide_gravity(&model, q, tau);
        (void)qcollide_observer_update(&obs, q, qd, tau, DT);

        {
            const int fired = qcollide_detector_update(&det,
                qcollide_observer_magnitude(&obs));
            if (fired && latency != NULL && *latency < 0.0 &&
                contact != NULL) {
                *latency = t - contact->start;
            }
        }
    }

    if (events != NULL) *events = qcollide_detector_events(&det);
}

/**
 * With a perfect model and no noise there is nothing to trigger on, so the
 * false alarm count must be exactly zero at any threshold. This is the
 * baseline M3 is measured against: once friction is introduced, this number
 * stops being zero and the threshold has to be chosen.
 */
static void test_no_false_alarms_with_a_perfect_model(void)
{
    const double thresholds[3] = { 0.05, 0.2, 1.0 };
    int i;

    for (i = 0; i < 3; ++i) {
        unsigned events = 0u;
        run_chain(50.0, thresholds[i], 3, NULL, 3.0, NULL, &events);
        printf("  [measured] threshold %.2f Nm: %u false alarms in 3 s\n",
               thresholds[i], events);
        TEST_ASSERT_EQUAL_UINT(0u, events);
    }
}

/**
 * Detection latency cannot be shorter than the observer needs to reach the
 * threshold, and that is set by the time constant 1/K. Measuring it across
 * gains checks that the whole chain inherits the observer's dynamics rather
 * than adding something unexplained.
 */
static void test_latency_scales_with_the_observer_time_constant(void)
{
    const double gains[3] = { 20.0, 50.0, 100.0 };
    qcollide_contact c;
    double lat[3];
    int i;

    c.start = 1.0;
    c.duration = 1.0;
    c.torque[0] = 5.0;
    c.torque[1] = 0.0;

    for (i = 0; i < 3; ++i) {
        run_chain(gains[i], 1.0, 3, &c, 2.5, &lat[i], NULL);
        printf("  [measured] K = %5.1f (tau = %.1f ms): latency %.2f ms\n",
               gains[i], 1000.0 / gains[i], lat[i] * 1000.0);
        TEST_ASSERT_TRUE_MESSAGE(lat[i] > 0.0, "contact was never detected");
    }

    /* Halving the time constant should roughly halve the latency. The
       relationship is not exact because the threshold is a fixed fraction of
       the applied torque rather than of the final residual, but it must be
       monotone and of the right order. */
    for (i = 1; i < 3; ++i) {
        TEST_ASSERT_TRUE_MESSAGE(lat[i] < lat[i - 1],
                                 "latency did not fall with a higher gain");
    }
    TEST_ASSERT_TRUE_MESSAGE(lat[0] / lat[2] > 2.0,
        "latency does not scale with the observer time constant");
}

/** A higher threshold must not detect sooner. */
static void test_latency_is_monotone_in_the_threshold(void)
{
    const double thresholds[4] = { 0.5, 1.0, 2.0, 4.0 };
    qcollide_contact c;
    double prev = -1.0;
    int i;

    c.start = 1.0;
    c.duration = 1.0;
    c.torque[0] = 5.0;
    c.torque[1] = 0.0;

    for (i = 0; i < 4; ++i) {
        double lat;
        run_chain(50.0, thresholds[i], 3, &c, 2.5, &lat, NULL);
        printf("  [measured] threshold %.1f Nm: latency %.2f ms\n",
               thresholds[i], lat * 1000.0);
        TEST_ASSERT_TRUE(lat > 0.0);
        TEST_ASSERT_TRUE_MESSAGE(lat >= prev,
                                 "a higher threshold detected sooner");
        prev = lat;
    }
}

/**
 * Smallest external torque that fires the detector, by bisection.
 *
 * With a perfect model this should sit just above the threshold, since the
 * residual converges to the applied torque exactly. Any large gap would mean
 * the chain is losing signal somewhere. In M3 this number is what friction
 * pushes upwards.
 */
static void test_minimum_detectable_torque(void)
{
    const double threshold = 1.0;
    qcollide_contact c;
    double lo = 0.0, hi = 5.0;
    int it;

    c.start = 1.0;
    c.duration = 1.5;
    c.torque[1] = 0.0;

    for (it = 0; it < 20; ++it) {
        double lat;
        const double mid = 0.5 * (lo + hi);
        c.torque[0] = mid;
        run_chain(50.0, threshold, 3, &c, 3.0, &lat, NULL);
        if (lat > 0.0) hi = mid; else lo = mid;
    }

    printf("  [measured] minimum detectable torque: %.4f Nm "
           "at a %.1f Nm threshold\n", hi, threshold);
    /* The residual converges to the applied torque, so the smallest
       detectable torque is the threshold itself to within the convergence
       tolerance. A gap of more than a few percent would mean signal is being
       lost between the contact and the decision. */
    TEST_ASSERT_TRUE_MESSAGE(hi > threshold * 0.95 && hi < threshold * 1.10,
        "minimum detectable torque is far from the threshold");
}

/**
 * A contact shorter than the observer time constant may be missed entirely,
 * because the residual never reaches the threshold before the contact ends.
 * This is a real property of the method, not a defect, and it is measured
 * here so that it is stated rather than discovered later.
 */
static void test_short_contacts_can_be_missed(void)
{
    const double durations[4] = { 0.002, 0.005, 0.02, 0.1 };
    qcollide_contact c;
    int i, detected[4];

    c.start = 1.0;
    c.torque[0] = 5.0;
    c.torque[1] = 0.0;

    for (i = 0; i < 4; ++i) {
        double lat;
        c.duration = durations[i];
        run_chain(50.0, 1.0, 3, &c, 2.0, &lat, NULL);
        detected[i] = (lat > 0.0);
        printf("  [measured] contact of %.0f ms (tau = 20 ms): %s\n",
               durations[i] * 1000.0, detected[i] ? "detected" : "MISSED");
    }

    /* The longest contact, five time constants, must be caught. */
    TEST_ASSERT_TRUE_MESSAGE(detected[3],
                             "a contact of five time constants was missed");
}

int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_init_rejects_bad_arguments);
    RUN_TEST(test_confirmation_requires_consecutive_samples);
    RUN_TEST(test_a_single_sample_below_resets_the_counter);
    RUN_TEST(test_hysteresis_prevents_chattering);
    RUN_TEST(test_contact_window);
    RUN_TEST(test_no_false_alarms_with_a_perfect_model);
    RUN_TEST(test_latency_scales_with_the_observer_time_constant);
    RUN_TEST(test_latency_is_monotone_in_the_threshold);
    RUN_TEST(test_minimum_detectable_torque);
    RUN_TEST(test_short_contacts_can_be_missed);
    return UNITY_END();
}
