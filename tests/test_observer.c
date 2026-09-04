/**
 * Unit tests for the momentum-based residual observer.
 *
 * The observer claims one thing: the residual converges to the external
 * joint torque, with a first-order response whose time constant is 1/K. Both
 * halves of that claim are tested directly, because both can fail silently.
 * A residual that converges to the wrong value still looks like a working
 * detector, and a time constant that does not match the gain means the
 * detection latency cannot be predicted from the tuning.
 *
 * The model given to the observer is the same one the plant uses, so these
 * tests measure the observer in isolation. Model error is introduced
 * deliberately in M3, where its effect can be attributed.
 */
#include "unity.h"
#include "qcollide/observer.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#define N 2
#define DT 1e-4

static qcollide_model model;
static qcollide_observer obs;

void setUp(void)
{
    TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK,
        qcollide_model_uniform(&model, N, 1.0, 0.6, 9.81));
}

void tearDown(void) {}

/**
 * Run the plant and observer together, holding the arm against gravity so it
 * stays near its starting pose. Returns the residual at the end.
 *
 * Gravity compensation is recomputed each step from the current pose, which
 * is what a real controller does; it is not exact station-keeping, and the
 * arm drifts slightly. That is deliberate: an observer tested only at rest
 * would not exercise the Coriolis term at all.
 */
static void run(double duration, const double *tau_ext,
                double contact_time, double *residual_out,
                double *first_cross_time, double cross_level)
{
    double q[N] = { 0.3, -0.5 };
    double qd[N] = { 0.0, 0.0 };
    double tau[N], ext[N];
    const int steps = (int)(duration / DT + 0.5);
    int k, i;

    for (i = 0; i < N; ++i) ext[i] = 0.0;

    (void)qcollide_gravity(&model, q, tau);
    TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK, qcollide_observer_reset(&obs, q, qd));

    if (first_cross_time != NULL) *first_cross_time = -1.0;

    for (k = 0; k < steps; ++k) {
        const double t = (double)k * DT;

        if (tau_ext != NULL && t >= contact_time) {
            for (i = 0; i < N; ++i) ext[i] = tau_ext[i];
        }

        TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK,
            qcollide_step(&model, q, qd, tau, ext, DT));
        (void)qcollide_gravity(&model, q, tau);
        TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK,
            qcollide_observer_update(&obs, q, qd, tau, DT));

        if (first_cross_time != NULL && *first_cross_time < 0.0 &&
            qcollide_observer_magnitude(&obs) >= cross_level) {
            *first_cross_time = t - contact_time;
        }
    }

    if (residual_out != NULL) {
        for (i = 0; i < N; ++i) {
            residual_out[i] = qcollide_observer_residual(&obs)[i];
        }
    }
}

/* ------------------------------------------------------------------ */
/* Contract                                                            */
/* ------------------------------------------------------------------ */

static void test_init_rejects_bad_arguments(void)
{
    qcollide_observer o;
    double g[N] = { 10.0, 10.0 };

    TEST_ASSERT_EQUAL_INT(QCOLLIDE_ERR_ARG,
        qcollide_observer_init(NULL, &model, g));
    TEST_ASSERT_EQUAL_INT(QCOLLIDE_ERR_ARG,
        qcollide_observer_init(&o, NULL, g));
    TEST_ASSERT_EQUAL_INT(QCOLLIDE_ERR_ARG,
        qcollide_observer_init(&o, &model, NULL));

    /* A gain of zero does not give a slow observer, it gives one whose
       residual is identically zero and which therefore never detects
       anything. A negative gain diverges. Both are rejected. */
    g[0] = 0.0;
    TEST_ASSERT_EQUAL_INT(QCOLLIDE_ERR_ARG,
        qcollide_observer_init(&o, &model, g));
    g[0] = -5.0;
    TEST_ASSERT_EQUAL_INT(QCOLLIDE_ERR_ARG,
        qcollide_observer_init(&o, &model, g));
}

static void test_update_rejects_bad_arguments(void)
{
    double q[N] = { 0.0, 0.0 }, qd[N] = { 0.0, 0.0 }, tau[N] = { 0.0, 0.0 };

    TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK,
        qcollide_observer_init_uniform(&obs, &model, 20.0));
    TEST_ASSERT_EQUAL_INT(QCOLLIDE_ERR_ARG,
        qcollide_observer_update(NULL, q, qd, tau, DT));
    TEST_ASSERT_EQUAL_INT(QCOLLIDE_ERR_ARG,
        qcollide_observer_update(&obs, NULL, qd, tau, DT));
    TEST_ASSERT_EQUAL_INT(QCOLLIDE_ERR_ARG,
        qcollide_observer_update(&obs, q, qd, tau, 0.0));
    TEST_ASSERT_EQUAL_INT(QCOLLIDE_ERR_ARG,
        qcollide_observer_update(&obs, q, qd, tau, -DT));
}

/**
 * The first update must produce no residual. There is no history to
 * integrate yet, and a spike at start-up is indistinguishable from a contact
 * at start-up, which on a real robot means an unexplained protective stop
 * every time it powers on.
 */
static void test_first_update_produces_no_residual(void)
{
    double q[N] = { 0.3, -0.5 }, qd[N] = { 0.4, 0.2 }, tau[N] = { 5.0, 1.0 };

    TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK,
        qcollide_observer_init_uniform(&obs, &model, 50.0));
    TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK,
        qcollide_observer_update(&obs, q, qd, tau, DT));
    TEST_ASSERT_EQUAL_DOUBLE(0.0, qcollide_observer_magnitude(&obs));
}

/* ------------------------------------------------------------------ */
/* The two claims                                                      */
/* ------------------------------------------------------------------ */

static void test_residual_stays_zero_without_contact(void)
{
    double r[N];

    TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK,
        qcollide_observer_init_uniform(&obs, &model, 50.0));
    run(2.0, NULL, 0.0, r, NULL, 0.0);

    printf("  [measured] residual after 2 s with no contact: %.3e, %.3e\n",
           r[0], r[1]);
    TEST_ASSERT_TRUE_MESSAGE(fabs(r[0]) < 1e-3 && fabs(r[1]) < 1e-3,
                             "residual drifts without contact");
}

static void test_residual_converges_to_external_torque(void)
{
    const double ext[N] = { 4.0, -2.0 };
    double r[N];

    TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK,
        qcollide_observer_init_uniform(&obs, &model, 50.0));
    /* One second of contact is twenty time constants at K = 50. */
    run(2.0, ext, 1.0, r, NULL, 0.0);

    printf("  [measured] residual %.4f, %.4f against applied %.1f, %.1f\n",
           r[0], r[1], ext[0], ext[1]);
    TEST_ASSERT_TRUE_MESSAGE(fabs(r[0] - ext[0]) < 0.05,
                             "residual did not converge on joint 0");
    TEST_ASSERT_TRUE_MESSAGE(fabs(r[1] - ext[1]) < 0.05,
                             "residual did not converge on joint 1");
}

/**
 * The response should be first order with time constant 1/K, which is what
 * makes detection latency predictable from the gain alone. Checked by the
 * time to reach 63.2 % of the final value, one time constant by definition.
 */
static void test_time_constant_matches_the_gain(void)
{
    const double gains[3] = { 20.0, 50.0, 100.0 };
    const double ext[N] = { 5.0, 0.0 };
    int g;

    for (g = 0; g < 3; ++g) {
        double cross, r[N];
        const double expected = 1.0 / gains[g];

        TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK,
            qcollide_observer_init_uniform(&obs, &model, gains[g]));
        /* 63.2 % of the applied torque is one time constant. */
        run(1.0 + 10.0 * expected, ext, 0.5, r, &cross, 0.632 * ext[0]);

        printf("  [measured] K = %5.1f: tau = %.4f s, expected %.4f s "
               "(ratio %.2f)\n", gains[g], cross, expected, cross / expected);
        TEST_ASSERT_TRUE_MESSAGE(cross > 0.0, "residual never crossed 63 %");
        TEST_ASSERT_TRUE_MESSAGE(cross / expected > 0.7 &&
                                 cross / expected < 1.4,
                                 "time constant does not match 1/K");
    }
}

/**
 * A contact on one joint must not produce a residual on the other. The two
 * joints are dynamically coupled, so this is not automatic: it holds only
 * because the residual equation decouples, and if the Coriolis transpose
 * were wrong it would not.
 */
static void test_residual_isolates_the_loaded_joint(void)
{
    const double ext[N] = { 6.0, 0.0 };
    double r[N];

    TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK,
        qcollide_observer_init_uniform(&obs, &model, 50.0));
    run(2.0, ext, 1.0, r, NULL, 0.0);

    printf("  [measured] loaded joint %.4f, unloaded joint %.4f\n", r[0], r[1]);
    TEST_ASSERT_TRUE_MESSAGE(fabs(r[0] - ext[0]) < 0.05,
                             "loaded joint did not converge");
    TEST_ASSERT_TRUE_MESSAGE(fabs(r[1]) < 0.1,
                             "unloaded joint shows a spurious residual");
}

/**
 * Using C instead of C transpose in the residual is the single most likely
 * implementation error, and it produces a detector that works at rest and
 * fails whenever the arm moves. This test drives the arm hard enough that
 * the two differ, and would fail if the transpose were dropped.
 */
static void test_residual_is_correct_while_moving(void)
{
    double q[N] = { 0.3, -0.5 };
    double qd[N] = { 2.0, -1.5 };   /* moving from the start */
    double tau[N] = { 0.0, 0.0 };
    double ext[N] = { 0.0, 0.0 };
    const int steps = (int)(1.0 / DT + 0.5);
    double worst = 0.0;
    int k;

    TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK,
        qcollide_observer_init_uniform(&obs, &model, 50.0));
    TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK, qcollide_observer_reset(&obs, q, qd));

    /* Free motion under gravity: no applied torque, no contact, so the
       residual must stay at zero throughout despite large Coriolis terms. */
    for (k = 0; k < steps; ++k) {
        TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK,
            qcollide_step(&model, q, qd, tau, ext, DT));
        TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK,
            qcollide_observer_update(&obs, q, qd, tau, DT));
        {
            const double a = qcollide_observer_magnitude(&obs);
            if (a > worst) worst = a;
        }
    }

    printf("  [measured] worst residual during free swing: %.4f Nm\n", worst);
    TEST_ASSERT_TRUE_MESSAGE(worst < 0.05,
        "residual is non-zero during free motion: check the Coriolis transpose");
}

/* ------------------------------------------------------------------ */
/* Reset                                                               */
/* ------------------------------------------------------------------ */

static void test_reset_clears_an_established_residual(void)
{
    const double ext[N] = { 5.0, 0.0 };
    double q[N] = { 0.3, -0.5 }, qd[N] = { 0.0, 0.0 };
    double r[N];

    TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK,
        qcollide_observer_init_uniform(&obs, &model, 50.0));
    run(2.0, ext, 1.0, r, NULL, 0.0);
    TEST_ASSERT_TRUE(fabs(r[0]) > 1.0);

    TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK, qcollide_observer_reset(&obs, q, qd));
    TEST_ASSERT_EQUAL_DOUBLE(0.0, qcollide_observer_magnitude(&obs));
}

int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_init_rejects_bad_arguments);
    RUN_TEST(test_update_rejects_bad_arguments);
    RUN_TEST(test_first_update_produces_no_residual);
    RUN_TEST(test_residual_stays_zero_without_contact);
    RUN_TEST(test_residual_converges_to_external_torque);
    RUN_TEST(test_time_constant_matches_the_gain);
    RUN_TEST(test_residual_isolates_the_loaded_joint);
    RUN_TEST(test_residual_is_correct_while_moving);
    RUN_TEST(test_reset_clears_an_established_residual);
    return UNITY_END();
}
