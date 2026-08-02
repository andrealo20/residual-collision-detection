/**
 * Tests for the sensor model.
 *
 * The interesting property here is not that quantisation quantises, but the
 * chain of consequences: encoder resolution sets velocity noise, velocity
 * noise reaches the residual, and filtering trades that noise for lag which
 * also reaches the residual. Each link is measured.
 *
 * The noise generator is tested too, because a Box-Muller implementation
 * with the wrong constant produces samples that look Gaussian and have the
 * wrong variance, and every noise figure downstream would then be wrong by a
 * factor nobody would notice.
 */
#include "unity.h"
#include "qcollide/observer.h"
#include "qcollide/sensors.h"

#include <math.h>
#include <stdio.h>

#define N 2
#define DT 1e-4

static qcollide_model plant;
static const double PI = 3.14159265358979323846;

void setUp(void)
{
    TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK,
        qcollide_model_uniform(&plant, N, 1.0, 0.6, 9.81));
}

void tearDown(void) {}

/* ------------------------------------------------------------------ */
/* Noise generator                                                     */
/* ------------------------------------------------------------------ */

static void test_gaussian_has_the_right_moments(void)
{
    qcollide_noise n;
    double sum = 0.0, sum_sq = 0.0, mean, var;
    const int count = 200000;
    int i;

    qcollide_noise_seed(&n, 42u);
    for (i = 0; i < count; ++i) {
        const double x = qcollide_noise_gaussian(&n);
        sum += x;
        sum_sq += x * x;
    }
    mean = sum / (double)count;
    var = sum_sq / (double)count - mean * mean;

    printf("  [measured] mean %.4f, variance %.4f over %d samples\n",
           mean, var, count);
    /* A Box-Muller with a wrong constant still looks Gaussian and has the
       wrong variance, which would silently rescale every noise figure in
       the project. */
    TEST_ASSERT_TRUE_MESSAGE(fabs(mean) < 0.01, "noise mean is not zero");
    TEST_ASSERT_TRUE_MESSAGE(fabs(var - 1.0) < 0.02,
                             "noise variance is not one");
}

static void test_noise_is_deterministic(void)
{
    qcollide_noise a, b;
    int i;

    qcollide_noise_seed(&a, 7u);
    qcollide_noise_seed(&b, 7u);
    for (i = 0; i < 1000; ++i) {
        /* The values are pulled into locals first. Calling a stateful
           generator inside an assertion macro is a mistake: the macro may
           reference its arguments more than once, and each reference
           advances the stream. The first version of this test did exactly
           that and failed while the generator was correct. */
        const double va = qcollide_noise_gaussian(&a);
        const double vb = qcollide_noise_gaussian(&b);
        TEST_ASSERT_EQUAL_DOUBLE(va, vb);
    }

    /* Different seeds must actually differ, or a "deterministic" generator
       would be a constant one. */
    qcollide_noise_seed(&b, 8u);
    {
        double diff = 0.0;
        qcollide_noise_seed(&a, 7u);
        for (i = 0; i < 100; ++i) {
            diff += fabs(qcollide_noise_gaussian(&a) -
                         qcollide_noise_gaussian(&b));
        }
        TEST_ASSERT_TRUE_MESSAGE(diff > 1.0, "different seeds give the same stream");
    }
}

/* ------------------------------------------------------------------ */
/* Quantisation                                                        */
/* ------------------------------------------------------------------ */

static void test_position_is_quantised_to_encoder_counts(void)
{
    qcollide_sensors s;
    qcollide_sensor_config cfg;
    double q[N] = { 0.123456789, -0.987654321 };
    double qd[N] = { 0.0, 0.0 }, tau[N] = { 0.0, 0.0 };
    double qm[N], qdm[N], taum[N];
    const double step = 2.0 * PI / 131072.0;
    int i;

    qcollide_sensors_typical(&cfg, 1.0);
    cfg.torque_noise = 0.0;
    TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK, qcollide_sensors_init(&s, &cfg, 1u));
    TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK,
        qcollide_sensors_measure(&s, N, q, qd, tau, qm, qdm, taum, DT));

    for (i = 0; i < N; ++i) {
        const double counts = qm[i] / step;
        TEST_ASSERT_DOUBLE_WITHIN(1e-9, floor(counts + 0.5), counts);
        /* And within half a count of the truth. */
        TEST_ASSERT_TRUE(fabs(qm[i] - q[i]) <= 0.5 * step + 1e-12);
    }
}

static void test_zero_encoder_counts_disables_quantisation(void)
{
    qcollide_sensors s;
    qcollide_sensor_config cfg;
    double q[N] = { 0.123456789, -0.987654321 };
    double qd[N] = { 0.0, 0.0 }, tau[N] = { 0.0, 0.0 };
    double qm[N], qdm[N], taum[N];

    qcollide_sensors_typical(&cfg, 1.0);
    cfg.encoder_counts = 0u;
    cfg.torque_noise = 0.0;
    TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK, qcollide_sensors_init(&s, &cfg, 1u));
    TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK,
        qcollide_sensors_measure(&s, N, q, qd, tau, qm, qdm, taum, DT));

    TEST_ASSERT_EQUAL_DOUBLE(q[0], qm[0]);
    TEST_ASSERT_EQUAL_DOUBLE(q[1], qm[1]);
}

/* ------------------------------------------------------------------ */
/* Velocity estimation                                                 */
/* ------------------------------------------------------------------ */

/**
 * Velocity from differencing quantised position, on a known constant-velocity
 * signal. The noise should match the theoretical bound: one encoder count
 * over one sample period.
 */
static void test_differenced_velocity_noise_matches_theory(void)
{
    qcollide_sensors s;
    qcollide_sensor_config cfg;
    const double true_vel = 0.5;
    const double step = 2.0 * PI / 131072.0;
    const double bound = step / DT;
    double q[N], qd[N] = { 0.5, 0.5 }, tau[N] = { 0.0, 0.0 };
    double qm[N], qdm[N], taum[N];
    double worst = 0.0;
    int k;

    qcollide_sensors_typical(&cfg, 1.0);
    cfg.torque_noise = 0.0;
    cfg.velocity_mode = QCOLLIDE_VEL_DIFF;
    TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK, qcollide_sensors_init(&s, &cfg, 1u));

    for (k = 0; k < 2000; ++k) {
        q[0] = true_vel * (double)k * DT;
        q[1] = q[0];
        (void)qcollide_sensors_measure(&s, N, q, qd, tau, qm, qdm, taum, DT);
        if (k > 2) {
            const double e = fabs(qdm[0] - true_vel);
            if (e > worst) worst = e;
        }
    }

    printf("  [measured] velocity error from differencing: %.3f rad/s "
           "(one count per sample is %.3f)\n", worst, bound);
    /* The difference is either zero or one count per sample, so the error
       cannot exceed one count per sample plus the true velocity's own
       fractional part. */
    TEST_ASSERT_TRUE_MESSAGE(worst <= bound + 1e-9,
                             "velocity error exceeds the quantisation bound");
    TEST_ASSERT_TRUE_MESSAGE(worst > 0.1 * bound,
        "no quantisation noise appeared: is quantisation actually enabled?");
}

static void test_filter_reduces_velocity_noise(void)
{
    const qcollide_velocity_mode modes[2] =
        { QCOLLIDE_VEL_DIFF, QCOLLIDE_VEL_FILTERED };
    double worst[2];
    int m, k;

    for (m = 0; m < 2; ++m) {
        qcollide_sensors s;
        qcollide_sensor_config cfg;
        double q[N], qd[N] = { 0.5, 0.5 }, tau[N] = { 0.0, 0.0 };
        double qm[N], qdm[N], taum[N];

        qcollide_sensors_typical(&cfg, 1.0);
        cfg.torque_noise = 0.0;
        cfg.velocity_mode = modes[m];
        (void)qcollide_sensors_init(&s, &cfg, 1u);

        worst[m] = 0.0;
        for (k = 0; k < 4000; ++k) {
            q[0] = 0.5 * (double)k * DT;
            q[1] = q[0];
            (void)qcollide_sensors_measure(&s, N, q, qd, tau, qm, qdm, taum, DT);
            if (k > 1000) {   /* let the filter settle */
                const double e = fabs(qdm[0] - 0.5);
                if (e > worst[m]) worst[m] = e;
            }
        }
    }

    printf("  [measured] velocity error: %.4f unfiltered, %.4f filtered\n",
           worst[0], worst[1]);
    TEST_ASSERT_TRUE_MESSAGE(worst[1] < 0.2 * worst[0],
                             "the filter does not reduce velocity noise");
}

/**
 * The filter buys noise reduction with lag, and lag is not free: on an
 * accelerating joint the estimate trails the truth. This measures the cost
 * so it is a stated trade rather than an unstated one.
 */
static void test_filter_introduces_lag(void)
{
    qcollide_sensors s;
    qcollide_sensor_config cfg;
    double q[N], qd[N] = { 0.0, 0.0 }, tau[N] = { 0.0, 0.0 };
    double qm[N], qdm[N], taum[N];
    double worst = 0.0;
    const double a = 4.0;      /* rad/s^2 */
    int k;

    qcollide_sensors_typical(&cfg, 1.0);
    cfg.torque_noise = 0.0;
    cfg.velocity_mode = QCOLLIDE_VEL_FILTERED;
    (void)qcollide_sensors_init(&s, &cfg, 1u);

    for (k = 0; k < 5000; ++k) {
        const double t = (double)k * DT;
        q[0] = 0.5 * a * t * t;
        q[1] = q[0];
        qd[0] = a * t;
        qd[1] = qd[0];
        (void)qcollide_sensors_measure(&s, N, q, qd, tau, qm, qdm, taum, DT);
        if (k > 2000) {
            const double e = fabs(qdm[0] - qd[0]);
            if (e > worst) worst = e;
        }
    }

    /* A first-order filter of time constant T trailing a ramp of slope a
       settles at a lag of a*T. At 50 Hz, T = 3.2 ms, so about 0.013 rad/s. */
    printf("  [measured] velocity lag under %.1f rad/s^2: %.4f rad/s "
           "(a*T predicts %.4f)\n", a, worst, a / (2.0 * PI * 50.0));
    TEST_ASSERT_TRUE_MESSAGE(worst > 0.5 * a / (2.0 * PI * 50.0) &&
                             worst < 3.0 * a / (2.0 * PI * 50.0),
        "velocity lag does not match the filter time constant");
}

/* ------------------------------------------------------------------ */
/* The M4 effect                                                       */
/* ------------------------------------------------------------------ */

static double residual_floor(qcollide_velocity_mode mode, double torque_noise)
{
    qcollide_observer obs;
    qcollide_sensors sens;
    qcollide_sensor_config cfg;
    double q[N] = { 0.4, -0.3 }, qd[N] = { 0.0, 0.0 };
    double tau[N], qdd[N], qm[N], qdm[N], taum[N];
    double worst = 0.0;
    const double w = 2.0 * PI * 0.4;
    int k, i;

    qcollide_sensors_typical(&cfg, 1.0);
    cfg.torque_noise = torque_noise;
    cfg.velocity_mode = mode;
    (void)qcollide_sensors_init(&sens, &cfg, 20260801u);
    (void)qcollide_observer_init_uniform(&obs, &plant, 50.0);

    for (k = 0; k < 20000; ++k) {
        const double t = (double)k * DT;
        for (i = 0; i < N; ++i) {
            qdd[i] = -0.6 * w * w * sin(w * t + (double)i * 0.7);
        }
        (void)qcollide_inverse_dynamics(&plant, q, qd, qdd, tau);
        (void)qcollide_step(&plant, q, qd, tau, NULL, DT);
        (void)qcollide_sensors_measure(&sens, N, q, qd, tau, qm, qdm, taum, DT);

        if (k == 0) {
            (void)qcollide_observer_reset(&obs, qm, qdm);
            continue;
        }
        (void)qcollide_observer_update(&obs, qm, qdm, taum, DT);
        if (t > 0.3) {
            const double m = qcollide_observer_magnitude(&obs);
            if (m > worst) worst = m;
        }
    }
    return worst;
}

/**
 * The M4 result: velocity estimation dominates the sensor contribution to
 * the residual, and torque noise is almost irrelevant beside it. That is not
 * the intuitive ordering -- the torque measurement is the one that appears
 * explicitly in the residual equation -- so it is worth measuring rather
 * than assuming.
 */
static void test_velocity_estimation_dominates_sensor_error(void)
{
    const double raw_diff   = residual_floor(QCOLLIDE_VEL_DIFF, 0.0);
    const double filtered   = residual_floor(QCOLLIDE_VEL_FILTERED, 0.0);
    const double exact_only = residual_floor(QCOLLIDE_VEL_EXACT, 0.05);
    const double both       = residual_floor(QCOLLIDE_VEL_FILTERED, 0.05);

    printf("  [measured] residual floor:\n");
    printf("    raw encoder difference        %8.3f Nm\n", raw_diff);
    printf("    filtered difference           %8.3f Nm\n", filtered);
    printf("    exact velocity + torque noise %8.3f Nm\n", exact_only);
    printf("    filtered + torque noise       %8.3f Nm\n", both);

    TEST_ASSERT_TRUE_MESSAGE(raw_diff > 10.0 * filtered,
        "filtering does not substantially reduce the residual floor");
    TEST_ASSERT_TRUE_MESSAGE(exact_only < 0.2 * filtered,
        "torque noise contributes more than velocity error, unexpectedly");
    /* Adding torque noise on top of filtered velocity changes almost
       nothing, which is the point. */
    TEST_ASSERT_TRUE_MESSAGE(both < 1.3 * filtered,
        "torque noise adds more than expected on top of velocity error");
}

int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_gaussian_has_the_right_moments);
    RUN_TEST(test_noise_is_deterministic);
    RUN_TEST(test_position_is_quantised_to_encoder_counts);
    RUN_TEST(test_zero_encoder_counts_disables_quantisation);
    RUN_TEST(test_differenced_velocity_noise_matches_theory);
    RUN_TEST(test_filter_reduces_velocity_noise);
    RUN_TEST(test_filter_introduces_lag);
    RUN_TEST(test_velocity_estimation_dominates_sensor_error);
    return UNITY_END();
}
