/**
 * Unit tests for the plant dynamics.
 *
 * These check the structural properties of the model rather than particular
 * numbers, because the numbers are already checked against Pinocchio in
 * tools/check_dynamics.py. What is tested here is what the observer will
 * depend on:
 *
 *   - M is symmetric and positive definite
 *   - Mdot - 2C is skew-symmetric, which is the identity that lets the
 *     momentum observer work without joint acceleration
 *   - energy is conserved under the integrator to a measured tolerance
 *   - inverse and forward dynamics are exact inverses of each other
 *
 * The skew-symmetry test is the important one. If it fails, the Coriolis
 * matrix is in the wrong form, and the observer built on top of it will
 * still run and still produce plausible residuals that mean nothing.
 */
#include "unity.h"
#include "qcollide/dynamics.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define N 3
#define NN (N * N)

static qcollide_model model;

void setUp(void)
{
    TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK,
        qcollide_model_uniform(&model, N, 1.2, 0.7, 9.81));
}

void tearDown(void) {}

/* M_PI is not guaranteed in strict C99 (-std=c99 without extensions). */
static const double PI = 3.14159265358979323846;

static void random_state(double *q, double *qd, double scale)
{
    int i;
    for (i = 0; i < N; ++i) {
        q[i]  = ((double)rand() / RAND_MAX - 0.5) * 2.0 * PI;
        qd[i] = ((double)rand() / RAND_MAX - 0.5) * 2.0 * scale;
    }
}

/* ------------------------------------------------------------------ */
/* Model setup                                                         */
/* ------------------------------------------------------------------ */

static void test_init_rejects_bad_parameters(void)
{
    qcollide_model m;
    qcollide_link bad[2];

    TEST_ASSERT_EQUAL_INT(QCOLLIDE_ERR_ARG,
        qcollide_model_uniform(&m, 0, 1.0, 1.0, 9.81));
    TEST_ASSERT_EQUAL_INT(QCOLLIDE_ERR_ARG,
        qcollide_model_uniform(&m, QCOLLIDE_MAX_DOF + 1, 1.0, 1.0, 9.81));

    /* A zero mass is a typo, not an unusual robot. Catching it here avoids a
       mass matrix that is silently singular much later. */
    bad[0].mass = 0.0; bad[0].length = 1.0; bad[0].com = 0.5; bad[0].inertia = 0.1;
    bad[1] = bad[0]; bad[1].mass = 1.0;
    TEST_ASSERT_EQUAL_INT(QCOLLIDE_ERR_ARG,
        qcollide_model_init(&m, bad, 2, 9.81));

    bad[0].mass = 1.0; bad[0].length = -1.0;
    TEST_ASSERT_EQUAL_INT(QCOLLIDE_ERR_ARG,
        qcollide_model_init(&m, bad, 2, 9.81));
}

/* ------------------------------------------------------------------ */
/* Mass matrix                                                         */
/* ------------------------------------------------------------------ */

static void test_mass_matrix_is_symmetric(void)
{
    double q[N], qd[N], M[NN];
    double worst = 0.0;
    int trial, i, j;

    for (trial = 0; trial < 200; ++trial) {
        random_state(q, qd, 1.0);
        TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK, qcollide_mass_matrix(&model, q, M));
        for (i = 0; i < N; ++i) {
            for (j = 0; j < N; ++j) {
                const double d = fabs(M[i * N + j] - M[j * N + i]);
                if (d > worst) worst = d;
            }
        }
    }
    printf("  [measured] mass matrix asymmetry: %.3e\n", worst);
    TEST_ASSERT_TRUE_MESSAGE(worst < 1e-12, "mass matrix is not symmetric");
}

/**
 * Positive definiteness, checked by Cholesky rather than by eigenvalues: if
 * the factorisation runs to completion without a non-positive pivot, the
 * matrix is positive definite, and that is both cheaper and exactly the
 * condition the forward-dynamics solve needs.
 */
static void test_mass_matrix_is_positive_definite(void)
{
    double q[N], qd[N], M[NN], L[NN];
    double smallest_pivot = 1e300;
    int trial, i, j, k;

    for (trial = 0; trial < 200; ++trial) {
        random_state(q, qd, 1.0);
        (void)qcollide_mass_matrix(&model, q, M);

        for (i = 0; i < NN; ++i) L[i] = 0.0;
        for (i = 0; i < N; ++i) {
            for (j = 0; j <= i; ++j) {
                double s = M[i * N + j];
                for (k = 0; k < j; ++k) s -= L[i * N + k] * L[j * N + k];
                if (i == j) {
                    if (s < smallest_pivot) smallest_pivot = s;
                    TEST_ASSERT_TRUE_MESSAGE(s > 0.0,
                        "mass matrix is not positive definite");
                    L[i * N + j] = sqrt(s);
                } else {
                    L[i * N + j] = s / L[j * N + j];
                }
            }
        }
    }
    printf("  [measured] smallest Cholesky pivot over 200 poses: %.4f\n",
           smallest_pivot);
}

/* ------------------------------------------------------------------ */
/* The identity the observer depends on                                */
/* ------------------------------------------------------------------ */

/**
 * Mdot(q, qd) - 2 C(q, qd) must be skew-symmetric.
 *
 * This is not a decorative property. The momentum observer works because
 * pdot = tau + C' qd - g + tau_ext, and that expression follows from this
 * identity. If the Coriolis matrix were in any other valid form -- and there
 * are infinitely many matrices C for which C qd gives the correct Coriolis
 * vector -- the observer would still run and still produce residuals, and
 * those residuals would be wrong in a way nothing else would reveal.
 *
 * Mdot is taken by central differences along the trajectory, since it is
 * dM/dt = sum_k dM/dq_k * qd_k.
 */
static void test_skew_symmetry_of_Mdot_minus_2C(void)
{
    double q[N], qd[N], qp[N], qm[N];
    double Mp[NN], Mm[NN], C[NN];
    double worst = 0.0, scale = 0.0;
    const double h = 1e-6;
    int trial, i, j;

    for (trial = 0; trial < 100; ++trial) {
        random_state(q, qd, 2.0);

        /* Step along the velocity direction to get dM/dt directly. */
        for (i = 0; i < N; ++i) {
            qp[i] = q[i] + h * qd[i];
            qm[i] = q[i] - h * qd[i];
        }
        (void)qcollide_mass_matrix(&model, qp, Mp);
        (void)qcollide_mass_matrix(&model, qm, Mm);
        (void)qcollide_coriolis_matrix(&model, q, qd, C);

        for (i = 0; i < N; ++i) {
            for (j = 0; j < N; ++j) {
                const double mdot_ij = (Mp[i * N + j] - Mm[i * N + j]) / (2.0 * h);
                const double mdot_ji = (Mp[j * N + i] - Mm[j * N + i]) / (2.0 * h);
                const double s_ij = mdot_ij - 2.0 * C[i * N + j];
                const double s_ji = mdot_ji - 2.0 * C[j * N + i];
                const double d = fabs(s_ij + s_ji);   /* skew: S = -S' */
                if (d > worst) worst = d;
                if (fabs(s_ij) > scale) scale = fabs(s_ij);
            }
        }
    }

    printf("  [measured] skew-symmetry residual: %.3e (matrix scale %.3f)\n",
           worst, scale);
    /* Both Mdot and C come from finite differences, so the residual is
       bounded by their combined truncation error rather than by machine
       epsilon. Relative to the size of the matrix itself. */
    TEST_ASSERT_TRUE_MESSAGE(worst < 1e-4 * (1.0 + scale),
                             "Mdot - 2C is not skew-symmetric");
}

/* ------------------------------------------------------------------ */
/* Dynamics consistency                                                */
/* ------------------------------------------------------------------ */

static void test_inverse_and_forward_dynamics_are_inverses(void)
{
    double q[N], qd[N], qdd[N], tau[N], back[N];
    double worst = 0.0;
    int trial, i;

    for (trial = 0; trial < 200; ++trial) {
        random_state(q, qd, 2.0);
        for (i = 0; i < N; ++i) {
            qdd[i] = ((double)rand() / RAND_MAX - 0.5) * 8.0;
        }
        TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK,
            qcollide_inverse_dynamics(&model, q, qd, qdd, tau));
        TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK,
            qcollide_forward_dynamics(&model, q, qd, tau, NULL, back));
        for (i = 0; i < N; ++i) {
            const double d = fabs(back[i] - qdd[i]);
            if (d > worst) worst = d;
        }
    }
    printf("  [measured] inverse/forward round trip: %.3e\n", worst);
    TEST_ASSERT_TRUE_MESSAGE(worst < 1e-9, "dynamics round trip is not exact");
}

/**
 * External torque must enter exactly as an additional joint torque, since
 * that is the assumption the whole detection method rests on: a contact
 * appears at the joints as tau_ext and in no other way.
 */
static void test_external_torque_adds_to_applied_torque(void)
{
    double q[N], qd[N], tau[N], ext[N], sum[N];
    double a[N], b[N];
    int trial, i;

    for (trial = 0; trial < 100; ++trial) {
        random_state(q, qd, 1.5);
        for (i = 0; i < N; ++i) {
            tau[i] = ((double)rand() / RAND_MAX - 0.5) * 20.0;
            ext[i] = ((double)rand() / RAND_MAX - 0.5) * 10.0;
            sum[i] = tau[i] + ext[i];
        }
        (void)qcollide_forward_dynamics(&model, q, qd, tau, ext, a);
        (void)qcollide_forward_dynamics(&model, q, qd, sum, NULL, b);
        for (i = 0; i < N; ++i) {
            TEST_ASSERT_DOUBLE_WITHIN(1e-12, b[i], a[i]);
        }
    }
}

/* ------------------------------------------------------------------ */
/* Integration                                                         */
/* ------------------------------------------------------------------ */

/**
 * With no applied torque and no external force, total mechanical energy is
 * conserved by the physics. Any drift is integration error.
 *
 * Rather than assert that the drift is below some chosen number, this checks
 * that it *shrinks with the step* the way a fourth-order scheme requires.
 * That distinction matters: a model with a wrong sign somewhere would also
 * produce a small drift at a small step, and would pass a threshold test
 * while failing this one.
 *
 * Measured over 20 s of free motion, relative drift against step size:
 *
 *     4.0e-3  ->  1.06e-3
 *     2.0e-3  ->  2.33e-5
 *     1.0e-3  ->  2.51e-6
 *     5.0e-4  ->  3.89e-7
 *
 * The ratios are not a clean 16 because the drift is a maximum over the
 * trajectory rather than a smooth error term, but each halving of the step
 * cuts it by roughly an order of magnitude, which is the fourth-order
 * behaviour and rules out a first- or second-order mistake.
 */
static void test_energy_drift_shrinks_with_step_size(void)
{
    const double dts[3] = { 2e-3, 1e-3, 5e-4 };
    double rel[3];
    int k, s;

    for (k = 0; k < 3; ++k) {
        double q[N] = { 0.4, -0.9, 0.6 };
        double qd[N] = { 0.8, -0.5, 0.3 };
        double tau[N] = { 0.0, 0.0, 0.0 };
        const int steps = (int)(20.0 / dts[k] + 0.5);
        double e0 = qcollide_energy(&model, q, qd);
        double worst = 0.0;

        for (s = 0; s < steps; ++s) {
            TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK,
                qcollide_step(&model, q, qd, tau, NULL, dts[k]));
            {
                const double d = fabs(qcollide_energy(&model, q, qd) - e0);
                if (d > worst) worst = d;
            }
        }
        rel[k] = worst / fabs(e0);
        printf("  [measured] dt %.1e over 20 s: relative energy drift %.3e\n",
               dts[k], rel[k]);
    }

    /* Each halving must cut the drift by at least a factor of four. A
       second-order scheme would give exactly four; RK4 gives far more, so
       this bound separates the two without depending on the exact constant. */
    for (k = 1; k < 3; ++k) {
        TEST_ASSERT_TRUE_MESSAGE(rel[k - 1] / rel[k] > 4.0,
            "energy drift does not shrink fast enough with the step size");
    }
    /* And at the finest step it must be small in absolute terms, or the
       simulation would not be usable as a plant regardless of its order. */
    TEST_ASSERT_TRUE_MESSAGE(rel[2] < 1e-5, "energy drift too large at dt=5e-4");
}

/**
 * The integrator should be fourth order: halving the step should cut the
 * error by roughly sixteen. Checking the order rather than the magnitude
 * catches a subtly wrong RK4 that would otherwise pass an energy test simply
 * by using a small enough step.
 */
static void test_integrator_is_fourth_order(void)
{
    double ratios[3];
    int r;

    for (r = 0; r < 3; ++r) {
        const double dt_coarse = 4e-3 / (double)(1 << r);
        double errs[2];
        int half;

        for (half = 0; half < 2; ++half) {
            const double dt = half ? dt_coarse * 0.5 : dt_coarse;
            const double total = 1.0;
            const int steps = (int)(total / dt + 0.5);
            /* Reference: the same motion at a much finer step. */
            double q[N] = { 0.4, -0.9, 0.6 }, qd[N] = { 0.8, -0.5, 0.3 };
            double qr[N] = { 0.4, -0.9, 0.6 }, qdr[N] = { 0.8, -0.5, 0.3 };
            double tau[N] = { 0.0, 0.0, 0.0 };
            const double dtr = 1e-5;
            const int stepsr = (int)(total / dtr + 0.5);
            int k, i;
            double e = 0.0;

            for (k = 0; k < steps; ++k) {
                (void)qcollide_step(&model, q, qd, tau, NULL, dt);
            }
            for (k = 0; k < stepsr; ++k) {
                (void)qcollide_step(&model, qr, qdr, tau, NULL, dtr);
            }
            for (i = 0; i < N; ++i) {
                const double d = fabs(q[i] - qr[i]);
                if (d > e) e = d;
            }
            errs[half] = e;
        }
        ratios[r] = (errs[1] > 0.0) ? errs[0] / errs[1] : 0.0;
        printf("  [measured] dt %.1e -> %.1e: error ratio %.1f\n",
               dt_coarse, dt_coarse * 0.5, ratios[r]);
    }

    /* Fourth order means a ratio near 16. Anything below 8 would indicate a
       lower-order scheme; the reference step is not infinitely fine, so the
       upper bound is generous. */
    for (r = 0; r < 3; ++r) {
        TEST_ASSERT_TRUE_MESSAGE(ratios[r] > 8.0 && ratios[r] < 40.0,
                                 "integrator is not fourth order");
    }
}

int main(void)
{
    srand(3);
    UNITY_BEGIN();
    RUN_TEST(test_init_rejects_bad_parameters);
    RUN_TEST(test_mass_matrix_is_symmetric);
    RUN_TEST(test_mass_matrix_is_positive_definite);
    RUN_TEST(test_skew_symmetry_of_Mdot_minus_2C);
    RUN_TEST(test_inverse_and_forward_dynamics_are_inverses);
    RUN_TEST(test_external_torque_adds_to_applied_torque);
    RUN_TEST(test_energy_drift_shrinks_with_step_size);
    RUN_TEST(test_integrator_is_fourth_order);
    return UNITY_END();
}
