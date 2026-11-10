/**
 * Tests for the benchmark arm in bench/arm.h: does it follow the trajectory
 * the README describes?
 *
 * The sweep and the trace used to drive the arm open loop from rest. Nothing
 * checked the resulting motion, and it was not the described one: without
 * friction the joints wound up through several radians, with friction the
 * second joint never moved. Every measured number sat on that motion. These
 * tests pin the motion itself, so a change to the controller or the starting
 * state cannot quietly move the benchmark onto a different experiment again.
 */
#include "unity.h"
#include "qcollide/friction.h"

#include "../bench/arm.h"

#include <math.h>

#define DT 1e-4

static qcollide_model plant;

void setUp(void)
{
    TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK,
        qcollide_model_uniform(&plant, ARM_N, 1.0, 0.6, 9.81));
}

void tearDown(void) {}

static void friction_torque(const double *q, const double *qd, double *tau_ext,
                            size_t n, void *user_data)
{
    (void)q;
    qcollide_friction_vector((const qcollide_friction *)user_data, n, qd,
                             tau_ext);
}

/**
 * Track for 5 s, two periods, and report the worst position error and the
 * peak speed of each joint.
 */
static void track(double friction_scale, double *err, double *speed)
{
    qcollide_friction fric;
    double q[ARM_N], qd[ARM_N], qdd[ARM_N], tau[ARM_N];
    double qref[ARM_N], qdref[ARM_N], qddref[ARM_N];
    int k, i;

    qcollide_friction_typical(&fric, friction_scale);
    arm_desired(0.0, q, qd, qdd);
    for (i = 0; i < ARM_N; ++i) {
        err[i] = 0.0;
        speed[i] = 0.0;
    }

    for (k = 0; k < 50000; ++k) {
        const double t = (double)k * DT;
        arm_command(&plant, t, q, qd, tau);
        TEST_ASSERT_EQUAL_INT(QCOLLIDE_OK,
            qcollide_step_ext(&plant, q, qd, tau, friction_torque, &fric, DT));
        arm_desired(t + DT, qref, qdref, qddref);
        for (i = 0; i < ARM_N; ++i) {
            if (fabs(q[i] - qref[i]) > err[i]) err[i] = fabs(q[i] - qref[i]);
            if (fabs(qd[i]) > speed[i]) speed[i] = fabs(qd[i]);
        }
    }
}

/** With no disturbance the arm sits on the trajectory, up to the hold of
 *  the torque across each step. */
static void test_tracks_exactly_without_friction(void)
{
    double err[ARM_N], speed[ARM_N];
    track(0.0, err, speed);
    TEST_ASSERT_DOUBLE_WITHIN(1e-4, 0.0, err[0]);
    TEST_ASSERT_DOUBLE_WITHIN(1e-4, 0.0, err[1]);
}

/** Under full friction both joints still swing through the trajectory.
 *  Stiction of 1.2 Nm against the PD stiffness leaves a few hundredths of a
 *  radian of tracking error, which does not change the experiment. */
static void test_both_joints_move_under_friction(void)
{
    const double peak = ARM_AMP * 2.0 * 3.14159265358979323846 * ARM_FREQ_HZ;
    double err[ARM_N], speed[ARM_N];
    track(1.0, err, speed);
    TEST_ASSERT_DOUBLE_WITHIN(0.05, 0.0, err[0]);
    TEST_ASSERT_DOUBLE_WITHIN(0.05, 0.0, err[1]);
    TEST_ASSERT_DOUBLE_WITHIN(0.1 * peak, peak, speed[0]);
    TEST_ASSERT_DOUBLE_WITHIN(0.1 * peak, peak, speed[1]);
}

int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_tracks_exactly_without_friction);
    RUN_TEST(test_both_joints_move_under_friction);
    return UNITY_END();
}
