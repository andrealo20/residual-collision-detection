/**
 * Tests for the friction model, parameter perturbation, and the effect both
 * have on the residual.
 *
 * The friction tests check physical properties rather than particular
 * numbers: friction must dissipate energy, oppose motion, and show the
 * Stribeck dip. Getting the sign backwards would produce a plant that gains
 * energy, which looks like an unstable controller rather than a sign error,
 * so it is worth an explicit test.
 *
 * The last group is the M3 result in miniature: model error raises the
 * residual floor, and it does so monotonically in the severity. That is the
 * property the whole threshold trade-off rests on, and if it did not hold
 * the sweep in bench/sweep.c would be measuring noise.
 */
#include "unity.h"
#include "qcollide/friction.h"
#include "qcollide/observer.h"

#include <math.h>
#include <stdio.h>

#define N 2
#define DT 1e-4

static qcollide_model plant;

void setUp(void)
{
    TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK,
        qcollide_model_uniform(&plant, N, 1.0, 0.6, 9.81));
}

void tearDown(void) {}

/* ------------------------------------------------------------------ */
/* Friction shape                                                      */
/* ------------------------------------------------------------------ */

/**
 * Friction must oppose motion at every velocity. A sign error here gives a
 * plant that gains energy, which presents as an unstable simulation and
 * sends the search in entirely the wrong direction.
 */
static void test_friction_opposes_motion(void)
{
    qcollide_friction f;
    int i;

    qcollide_friction_typical(&f, 1.0);

    for (i = -40; i <= 40; ++i) {
        const double qd = 0.1 * (double)i;
        const double tau = qcollide_friction_torque(&f, qd);
        if (fabs(qd) < 1e-9) continue;
        /* Opposing motion means the product with velocity is negative,
           i.e. friction removes power from the joint. */
        TEST_ASSERT_TRUE_MESSAGE(tau * qd < 0.0,
                                 "friction does not oppose motion");
    }
}

static void test_friction_is_zero_at_zero_velocity(void)
{
    qcollide_friction f;
    qcollide_friction_typical(&f, 1.0);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.0, qcollide_friction_torque(&f, 0.0));
}

static void test_friction_is_odd_in_velocity(void)
{
    qcollide_friction f;
    int i;

    qcollide_friction_typical(&f, 1.0);
    for (i = 1; i <= 40; ++i) {
        const double qd = 0.1 * (double)i;
        TEST_ASSERT_DOUBLE_WITHIN(1e-12,
            -qcollide_friction_torque(&f, qd),
             qcollide_friction_torque(&f, -qd));
    }
}

/**
 * The Stribeck dip: friction falls from the static level towards Coulomb as
 * the joint picks up speed, before viscous drag takes over and it rises
 * again. Without the dip the model is just Coulomb plus viscous, which is a
 * different and easier problem.
 */
static void test_stribeck_dip_is_present(void)
{
    qcollide_friction f;
    double near_zero, dip, fast;

    qcollide_friction_typical(&f, 1.0);
    near_zero = fabs(qcollide_friction_torque(&f, 0.02));
    dip       = fabs(qcollide_friction_torque(&f, 0.30));
    fast      = fabs(qcollide_friction_torque(&f, 3.00));

    printf("  [measured] friction magnitude: %.3f at 0.02, %.3f at 0.30, "
           "%.3f at 3.0 rad/s\n", near_zero, dip, fast);
    TEST_ASSERT_TRUE_MESSAGE(dip < near_zero,
                             "no Stribeck dip: friction should fall first");
    TEST_ASSERT_TRUE_MESSAGE(fast > dip,
                             "viscous term does not dominate at speed");
}

static void test_friction_scales_linearly(void)
{
    qcollide_friction a, b;
    int i;

    qcollide_friction_typical(&a, 1.0);
    qcollide_friction_typical(&b, 2.0);
    for (i = 1; i <= 20; ++i) {
        const double qd = 0.15 * (double)i;
        TEST_ASSERT_DOUBLE_WITHIN(1e-9,
            2.0 * qcollide_friction_torque(&a, qd),
                  qcollide_friction_torque(&b, qd));
    }
}

/** Friction must remove mechanical energy from a freely swinging arm. */
static void test_friction_dissipates_energy(void)
{
    qcollide_friction f;
    double q[N] = { 0.9, -0.4 }, qd[N] = { 0.0, 0.0 };
    double tau[N] = { 0.0, 0.0 }, fr[N];
    double e0, e1;
    int k;

    qcollide_friction_typical(&f, 1.0);
    e0 = qcollide_energy(&plant, q, qd);

    for (k = 0; k < 30000; ++k) {
        qcollide_friction_vector(&f, N, qd, fr);
        TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK,
            qcollide_step(&plant, q, qd, tau, fr, DT));
    }
    e1 = qcollide_energy(&plant, q, qd);

    printf("  [measured] energy %.3f -> %.3f J over 3 s of free swing\n", e0, e1);
    TEST_ASSERT_TRUE_MESSAGE(e1 < e0, "friction did not dissipate energy");
}

/* ------------------------------------------------------------------ */
/* Parameter perturbation                                              */
/* ------------------------------------------------------------------ */

static void test_perturb_scales_mass_and_inertia(void)
{
    qcollide_model p;
    size_t i;

    TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK,
        qcollide_model_perturb(&p, &plant, 0.05));
    for (i = 0; i < plant.n; ++i) {
        TEST_ASSERT_DOUBLE_WITHIN(1e-12,
            plant.link[i].mass * 1.05, p.link[i].mass);
        TEST_ASSERT_DOUBLE_WITHIN(1e-12,
            plant.link[i].inertia * 1.05, p.link[i].inertia);
        /* Geometry is unchanged: this models a mass estimate that is wrong,
           not a robot that is a different shape. */
        TEST_ASSERT_DOUBLE_WITHIN(1e-12,
            plant.link[i].length, p.link[i].length);
    }
}

static void test_perturb_rejects_negative_mass(void)
{
    qcollide_model p;
    TEST_ASSERT_EQUAL_INT(QCOLLIDE_ERR_ARG,
        qcollide_model_perturb(&p, &plant, -1.0));
    TEST_ASSERT_EQUAL_INT(QCOLLIDE_ERR_ARG,
        qcollide_model_perturb(&p, &plant, -1.5));
    TEST_ASSERT_EQUAL_INT(QCOLLIDE_ERR_ARG,
        qcollide_model_perturb(NULL, &plant, 0.05));
}

/* ------------------------------------------------------------------ */
/* The M3 effect                                                       */
/* ------------------------------------------------------------------ */

/**
 * Run the plant with a given model error and return the peak residual with
 * no contact present. This is the floor any threshold has to clear.
 */
static double residual_floor(double friction_scale, double param_error)
{
    qcollide_model belief;
    qcollide_friction f;
    qcollide_observer obs;
    double q[N] = { 0.4, -0.3 }, qd[N] = { 0.0, 0.0 };
    double tau[N], fr[N], qdd[N];
    double worst = 0.0;
    const double w = 2.0 * 3.14159265358979323846 * 0.4;
    int k, i;

    (void)qcollide_model_perturb(&belief, &plant, param_error);
    qcollide_friction_typical(&f, friction_scale);
    (void)qcollide_observer_init_uniform(&obs, &belief, 50.0);

    for (k = 0; k < 20000; ++k) {
        const double t = (double)k * DT;
        for (i = 0; i < N; ++i) {
            qdd[i] = -0.6 * w * w * sin(w * t + (double)i * 0.7);
        }
        (void)qcollide_inverse_dynamics(&plant, q, qd, qdd, tau);
        qcollide_friction_vector(&f, N, qd, fr);
        (void)qcollide_step(&plant, q, qd, tau, fr, DT);

        if (k == 0) {
            (void)qcollide_observer_reset(&obs, q, qd);
            continue;
        }
        (void)qcollide_observer_update(&obs, q, qd, tau, DT);
        if (t > 0.3) {
            const double m = qcollide_observer_magnitude(&obs);
            if (m > worst) worst = m;
        }
    }
    return worst;
}

/**
 * The central claim of M3: model error raises the residual floor, and it
 * rises monotonically with the severity. If it did not, the threshold
 * trade-off measured in bench/sweep.c would be measuring noise.
 */
static void test_friction_raises_the_residual_floor(void)
{
    const double scales[4] = { 0.0, 0.25, 0.5, 1.0 };
    double floors[4];
    int i;

    for (i = 0; i < 4; ++i) {
        floors[i] = residual_floor(scales[i], 0.0);
        printf("  [measured] friction x%.2f: residual floor %.3f Nm\n",
               scales[i], floors[i]);
    }

    TEST_ASSERT_TRUE_MESSAGE(floors[0] < 0.01,
        "residual is not near zero with a perfect model");
    for (i = 1; i < 4; ++i) {
        TEST_ASSERT_TRUE_MESSAGE(floors[i] > floors[i - 1],
            "residual floor did not rise with friction");
    }
}

static void test_parameter_error_raises_the_residual_floor(void)
{
    const double errors[4] = { 0.0, 0.02, 0.05, 0.10 };
    double floors[4];
    int i;

    for (i = 0; i < 4; ++i) {
        floors[i] = residual_floor(0.0, errors[i]);
        printf("  [measured] inertia +%.0f%%: residual floor %.3f Nm\n",
               errors[i] * 100.0, floors[i]);
    }

    for (i = 1; i < 4; ++i) {
        TEST_ASSERT_TRUE_MESSAGE(floors[i] > floors[i - 1],
            "residual floor did not rise with parameter error");
    }
    /* Roughly proportional: doubling the error should roughly double the
       floor, since the residual carries the error linearly. */
    printf("  [measured] floor ratio 10%%/5%%: %.2f (linear would be 2.0)\n",
           floors[3] / floors[2]);
    TEST_ASSERT_TRUE_MESSAGE(floors[3] / floors[2] > 1.5 &&
                             floors[3] / floors[2] < 2.5,
        "parameter error does not enter the residual roughly linearly");
}

int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_friction_opposes_motion);
    RUN_TEST(test_friction_is_zero_at_zero_velocity);
    RUN_TEST(test_friction_is_odd_in_velocity);
    RUN_TEST(test_stribeck_dip_is_present);
    RUN_TEST(test_friction_scales_linearly);
    RUN_TEST(test_friction_dissipates_energy);
    RUN_TEST(test_perturb_scales_mass_and_inertia);
    RUN_TEST(test_perturb_rejects_negative_mass);
    RUN_TEST(test_friction_raises_the_residual_floor);
    RUN_TEST(test_parameter_error_raises_the_residual_floor);
    return UNITY_END();
}
