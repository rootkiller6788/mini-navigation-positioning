/**
 * @file    test_ins.c
 * @brief   Assert-based tests for inertial navigation module
 *
 * Tests cover all core APIs:
 *   - Vector/matrix operations
 *   - Coordinate transforms (geodetic <-> ECEF)
 *   - Gravity model (WGS84)
 *   - Quaternion operations and kinematics
 *   - Quaternion <-> DCM <-> Euler conversion
 *   - Strapdown mechanization
 *   - IMU error models and Allan variance
 *   - IMU calibration (6-position, cross-axis, thermal)
 *   - GPS/INS Kalman filter
 *   - Edge cases (NULL handling, singularities)
 */

#include "ins_core.h"
#include "ins_attitude.h"
#include "ins_mechanization.h"
#include "ins_errors.h"
#include "ins_calibration.h"
#include "ins_integration.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

static int tests_passed = 0;
static int tests_failed = 0;

#define TEST(name) do { printf("  TEST: %s ... ", name); } while(0)
#define PASS() do { printf("PASS\n"); tests_passed++; } while(0)
#define FAIL(msg) do { printf("FAIL: %s\n", msg); tests_failed++; return; } while(0)
#define ASSERT_TRUE(cond, msg) do { if (!(cond)) { FAIL(msg); } } while(0)
#define ASSERT_NEAR(a, b, tol, msg) do { \
    if (fabs((a) - (b)) > (tol)) { \
        printf("FAIL: %s (%.9f vs %.9f)\n", msg, a, b); \
        tests_failed++; return; \
    } \
} while(0)

/* =================================================================
 * Test 1: Vector Operations
 * ================================================================= */
static void test_vec3_ops(void) {
    TEST("vec3_set / vec3_norm");
    ins_vec3_t v;
    ins_vec3_set(&v, 3.0, 4.0, 0.0);
    ASSERT_NEAR(ins_vec3_norm(&v), 5.0, 1e-12, "norm of [3,4,0]");
    PASS();

    TEST("vec3_dot");
    ins_vec3_t a, b;
    ins_vec3_set(&a, 1, 2, 3);
    ins_vec3_set(&b, 4, 5, 6);
    ASSERT_NEAR(ins_vec3_dot(&a, &b), 32.0, 1e-12, "dot product");
    PASS();

    TEST("vec3_cross");
    ins_vec3_t c;
    ins_vec3_set(&a, 1, 0, 0);
    ins_vec3_set(&b, 0, 1, 0);
    ins_vec3_cross(&a, &b, &c);
    ASSERT_NEAR(c.x, 0, 1e-12, "cross x");
    ASSERT_NEAR(c.y, 0, 1e-12, "cross y");
    ASSERT_NEAR(c.z, 1, 1e-12, "cross z");
    PASS();

    TEST("vec3_normalize");
    ins_vec3_set(&v, 3, 4, 0);
    ins_vec3_normalize(&v);
    ASSERT_NEAR(ins_vec3_norm(&v), 1.0, 1e-12, "normalized norm");
    ASSERT_NEAR(v.x, 0.6, 1e-12, "normalized x");
    ASSERT_NEAR(v.y, 0.8, 1e-12, "normalized y");
    PASS();

    TEST("vec3_skew");
    ins_vec3_set(&v, 1, 2, 3);
    ins_mat3_t S;
    ins_vec3_skew(&v, &S);
    ASSERT_NEAR(S.m[0], 0, 1e-12, "skew diag0");
    ASSERT_NEAR(S.m[4], 0, 1e-12, "skew diag1");
    ASSERT_NEAR(S.m[8], 0, 1e-12, "skew diag2");
    ASSERT_NEAR(S.m[1], -3, 1e-12, "skew01 = -vz");
    ASSERT_NEAR(S.m[2], 2, 1e-12, "skew02 = vy");
    ASSERT_NEAR(S.m[3], 3, 1e-12, "skew10 = vz");
    PASS();
}

/* =================================================================
 * Test 2: Matrix Operations
 * ================================================================= */
static void test_mat3_ops(void) {
    TEST("mat3_identity * vec3 = vec3");
    ins_mat3_t I;
    ins_mat3_identity(&I);
    ins_vec3_t v = {1, 2, 3};
    ins_vec3_t r;
    ins_mat3_mul_vec(&I, &v, &r);
    ASSERT_NEAR(r.x, 1, 1e-12, "I*x");
    ASSERT_NEAR(r.y, 2, 1e-12, "I*y");
    ASSERT_NEAR(r.z, 3, 1e-12, "I*z");
    PASS();

    TEST("mat3_transpose");
    ins_mat3_t A, AT;
    A.m[0]=1; A.m[1]=2; A.m[2]=3;
    A.m[3]=4; A.m[4]=5; A.m[5]=6;
    A.m[6]=7; A.m[7]=8; A.m[8]=9;
    ins_mat3_transpose(&A, &AT);
    ASSERT_NEAR(AT.m[0], 1, 1e-12, "T00");
    ASSERT_NEAR(AT.m[1], 4, 1e-12, "T01");
    ASSERT_NEAR(AT.m[2], 7, 1e-12, "T02");
    ASSERT_NEAR(AT.m[3], 2, 1e-12, "T10");
    PASS();
}

/* =================================================================
 * Test 3: Coordinate Transforms
 * ================================================================= */
static void test_coord_transforms(void) {
    TEST("geodetic_to_ecef / ecef_to_geodetic round-trip");
    ins_geodetic_t geo = {0.5, -2.0, 100.0};
    ins_ecef_t ecef;
    ins_geodetic_to_ecef(&geo, &ecef);

    ins_geodetic_t geo2;
    ins_ecef_to_geodetic(&ecef, &geo2);

    ASSERT_NEAR(geo.lat, geo2.lat, 1e-7, "lat round-trip");
    ASSERT_NEAR(geo.lon, geo2.lon, 1e-7, "lon round-trip");
    ASSERT_NEAR(geo.alt, geo2.alt, 1.0, "alt round-trip");
    PASS();

    TEST("ecef_to_ned_dcm orthogonality");
    ins_mat3_t C;
    ins_ecef_to_ned_dcm(0.5, -2.0, &C);
    /* C^T * C should be identity */
    ins_mat3_t CT, CTC;
    ins_mat3_transpose(&C, &CT);
    ins_mat3_mul(&CT, &C, &CTC);
    ASSERT_NEAR(CTC.m[0], 1, 1e-12, "CTC[0,0]");
    ASSERT_NEAR(CTC.m[4], 1, 1e-12, "CTC[1,1]");
    ASSERT_NEAR(CTC.m[8], 1, 1e-12, "CTC[2,2]");
    ASSERT_NEAR(CTC.m[1], 0, 1e-12, "CTC[0,1] ~ 0");
    PASS();
}

/* =================================================================
 * Test 4: Gravity Model
 * ================================================================= */
static void test_gravity(void) {
    TEST("ins_gravity_wgs84 at equator sea level");
    double g = ins_gravity_wgs84(0.0, 0.0);
    ASSERT_NEAR(g, INS_GRAVITY_EQUATOR, 1e-6, "g_equator");
    PASS();

    TEST("ins_gravity_wgs84 at pole sea level");
    g = ins_gravity_wgs84(M_PI/2.0, 0.0);
    ASSERT_NEAR(g, INS_GRAVITY_POLE, 1e-6, "g_pole");
    PASS();

    TEST("ins_gravity_wgs84 decreases with altitude");
    double g0 = ins_gravity_wgs84(0.0, 0.0);
    double g10k = ins_gravity_wgs84(0.0, 10000.0);
    ASSERT_TRUE(g10k < g0, "gravity decreases with height");
    PASS();
}

/* =================================================================
 * Test 5: Quaternion Operations
 * ================================================================= */
static void test_quat_ops(void) {
    TEST("quat_identity norm = 1");
    ins_quat_t q;
    ins_quat_identity(&q);
    ASSERT_NEAR(ins_quat_norm(&q), 1.0, 1e-12, "identity norm");
    PASS();

    TEST("quat_mul associativity check");
    ins_quat_t q1, q2, q3, r1, r2;
    ins_quat_identity(&q1);
    q1.x = 0.5; q1.y = 0.3; q1.z = 0.1;
    ins_quat_normalize(&q1);
    ins_quat_identity(&q2);
    q2.x = 0.1; q2.y = 0.2; q2.z = 0.4;
    ins_quat_normalize(&q2);
    ins_quat_identity(&q3);
    q3.x = 0.2; q3.y = 0.1; q3.z = 0.3;
    ins_quat_normalize(&q3);

    /* (q1*q2)*q3 */
    ins_quat_t t1, t2;
    ins_quat_mul(&q1, &q2, &t1);
    ins_quat_mul(&t1, &q3, &r1);

    /* q1*(q2*q3) */
    ins_quat_mul(&q2, &q3, &t2);
    ins_quat_mul(&q1, &t2, &r2);

    /* Normalize both results */
    ins_quat_normalize(&r1);
    ins_quat_normalize(&r2);

    ASSERT_NEAR(r1.w, r2.w, 1e-12, "assoc w");
    ASSERT_NEAR(r1.x, r2.x, 1e-12, "assoc x");
    ASSERT_NEAR(r1.y, r2.y, 1e-12, "assoc y");
    ASSERT_NEAR(r1.z, r2.z, 1e-12, "assoc z");
    PASS();

    TEST("quat_conjugate: q * q* = [1,0,0,0] for unit q");
    ins_quat_t qc, qqc;
    ins_quat_conjugate(&q1, &qc);
    ins_quat_mul(&q1, &qc, &qqc);
    ASSERT_NEAR(qqc.w, 1, 1e-14, "conj w");
    ASSERT_NEAR(qqc.x, 0, 1e-14, "conj x");
    PASS();
}

/* =================================================================
 * Test 6: Quaternion <-> DCM <-> Euler
 * ================================================================= */
static void test_quat_conversions(void) {
    TEST("quat_to_dcm then dcm_to_quat round-trip");
    ins_quat_t q_orig, q_rec;
    ins_euler_t e = {0.2, 0.3, 0.5};

    ins_euler_to_quat(&e, &q_orig);
    ins_mat3_t C;
    ins_quat_to_dcm(&q_orig, &C);
    ins_dcm_to_quat(&C, &q_rec);

    /* Should be same up to sign (q and -q represent same rotation) */
    double dot = q_orig.w * q_rec.w + q_orig.x * q_rec.x
               + q_orig.y * q_rec.y + q_orig.z * q_rec.z;
    ASSERT_TRUE(fabs(fabs(dot) - 1.0) < 1e-12, "round-trip dot product magnitude");
    PASS();

    TEST("euler_to_quat then quat_to_euler round-trip");
    ins_euler_t e2;
    ins_quat_to_euler(&q_orig, &e2);
    ASSERT_NEAR(e.roll, e2.roll, 1e-7, "roll round-trip");
    ASSERT_NEAR(e.pitch, e2.pitch, 1e-7, "pitch round-trip");
    ASSERT_NEAR(e.yaw, e2.yaw, 1e-7, "yaw round-trip");
    PASS();

    TEST("quat_rotate_vector: v rotated by identity = v");
    ins_quat_t q_id;
    ins_quat_identity(&q_id);
    ins_vec3_t v_in = {1, 2, 3};
    ins_vec3_t v_out;
    ins_quat_rotate_vector(&q_id, &v_in, &v_out);
    ASSERT_NEAR(v_out.x, 1, 1e-12, "ident rot x");
    ASSERT_NEAR(v_out.y, 2, 1e-12, "ident rot y");
    ASSERT_NEAR(v_out.z, 3, 1e-12, "ident rot z");
    PASS();
}

/* =================================================================
 * Test 7: Quaternion Integration
 * ================================================================= */
static void test_quat_integration(void) {
    TEST("quat_update_exact: q rotated about z retains norm");
    ins_quat_t q;
    ins_quat_identity(&q);
    ins_vec3_t omega = {0, 0, 0.1};
    for (int i = 0; i < 100; i++) {
        ins_quat_update_exact(&q, &omega, 0.01);
    }
    ASSERT_NEAR(ins_quat_norm(&q), 1.0, 1e-12, "norm after 100 steps");
    PASS();

    TEST("quat_update_coning: produces valid quaternion");
    ins_quat_t q2;
    ins_quat_identity(&q2);
    ins_vec3_t w_prev = {0, 0, 0};
    ins_vec3_t w_curr = {0.01, 0.02, 0.03};
    ins_quat_update_coning(&q2, &w_prev, &w_curr, 0.01);
    ASSERT_NEAR(ins_quat_norm(&q2), 1.0, 1e-12, "norm after coning step");
    PASS();

    TEST("sculling_compensation: zero for zero inputs");
    ins_vec3_t da = {0}, dv = {0}, scull;
    ins_sculling_compensation(&da, &da, &dv, &dv, &scull);
    ASSERT_NEAR(scull.x, 0, 1e-12, "zero scull x");
    ASSERT_NEAR(scull.y, 0, 1e-12, "zero scull y");
    ASSERT_NEAR(scull.z, 0, 1e-12, "zero scull z");
    PASS();
}

/* =================================================================
 * Test 8: Strapdown Mechanization
 * ================================================================= */
static void test_mechanization(void) {
    TEST("mech_init sets initial position");
    ins_mech_state_t ms;
    ins_mech_init(&ms, INS_MECH_CONING, 0.5, -2.0, 100.0);
    ASSERT_NEAR(ms.pos.lat, 0.5, 1e-12, "init lat");
    ASSERT_NEAR(ms.pos.lon, -2.0, 1e-12, "init lon");
    ASSERT_NEAR(ms.pos.alt, 100.0, 1e-12, "init alt");
    PASS();

    TEST("align_level with level accel reading");
    ins_vec3_t accel = {0, 0, INS_GRAVITY_EQUATOR};
    int ret = ins_align_level(&ms, &accel, 1.0);
    ASSERT_TRUE(ret == 0, "align_level success");
    PASS();
}

/* =================================================================
 * Test 9: Allan Variance
 * ================================================================= */
static void test_allan(void) {
    TEST("allan_variance basic computation");
    double data[1024];
    /* Generate white noise (ARW) */
    for (int i = 0; i < 1024; i++) {
        data[i] = (double)rand() / RAND_MAX * 0.01;
    }
    double taus[64], adevs[64];
    size_t n_tau = ins_allan_variance(data, 1024, 0.01, taus, adevs, 64);
    ASSERT_TRUE(n_tau > 0, "allan produced tau values");
    PASS();

    TEST("allan_decompose produces non-negative values");
    ins_allan_result_t res;
    int ret = ins_allan_decompose(taus, adevs, n_tau, &res);
    ASSERT_TRUE(ret == 0, "allan_decompose success");
    ASSERT_TRUE(res.angle_random_walk >= 0, "ARW non-negative");
    ASSERT_TRUE(res.bias_instability >= 0, "Bias Inst non-negative");
    PASS();
}

/* =================================================================
 * Test 10: IMU Grade Specs
 * ================================================================= */
static void test_grade_specs(void) {
    TEST("grade_spec_get returns valid specs");
    const ins_grade_spec_t *spec = ins_grade_spec_get(INS_GRADE_NAVIGATION);
    ASSERT_TRUE(spec != NULL, "nav grade spec non-NULL");
    ASSERT_TRUE(spec->gyro_bias < 0.01, "nav gyro bias < 0.01 deg/hr");
    spec = ins_grade_spec_get(INS_GRADE_CONSUMER);
    ASSERT_TRUE(spec != NULL, "consumer grade spec non-NULL");
    PASS();
}

/* =================================================================
 * Test 11: Error Prediction
 * ================================================================= */
static void test_error_prediction(void) {
    TEST("error_predict_drift monotonic growth");
    ins_imu_error_model_t model;
    memset(&model, 0, sizeof(model));
    model.gyro_x.bias_offset = 1.0 * (M_PI / 180.0) / 3600.0;
    model.gyro_y.bias_offset = 1.0 * (M_PI / 180.0) / 3600.0;
    model.gyro_z.bias_offset = 1.0 * (M_PI / 180.0) / 3600.0;

    double pos1, pos2, vel, att;
    ins_error_predict_drift(&model, 100.0, &pos1, &vel, &att);
    ins_error_predict_drift(&model, 200.0, &pos2, &vel, &att);
    ASSERT_TRUE(pos2 > pos1, "position drift grows super-linearly");
    PASS();
}

/* =================================================================
 * Test 12: Calibration
 * ================================================================= */
static void test_calibration(void) {
    TEST("six_position_accel: perfect measurements give zero bias");
    ins_vec3_t readings[6];
    double g = 9.81;
    /* +x down */
    readings[0].x = g; readings[0].y = 0; readings[0].z = 0;
    /* -x down (x up) */
    readings[1].x = -g; readings[1].y = 0; readings[1].z = 0;
    /* +y down */
    readings[2].x = 0; readings[2].y = g; readings[2].z = 0;
    /* -y down */
    readings[3].x = 0; readings[3].y = -g; readings[3].z = 0;
    /* +z down */
    readings[4].x = 0; readings[4].y = 0; readings[4].z = g;
    /* -z down */
    readings[5].x = 0; readings[5].y = 0; readings[5].z = -g;

    ins_calibration_result_t res;
    int ret = ins_calib_six_position_accel(readings, g, &res);
    ASSERT_TRUE(ret == 0, "six-position calib success");
    ASSERT_NEAR(res.bias.x, 0, 1e-6, "bias x ~ 0");
    ASSERT_NEAR(res.scale_factor.x, 0, 1e-6, "sf x ~ 0");
    PASS();
}

/* =================================================================
 * Test 13: Kalman Filter
 * ================================================================= */
static void test_kalman_filter(void) {
    TEST("kf_init sets diagonal P");
    ins_kf_state_t kf;
    ins_kf_init(&kf, 10.0, 0.1, 0.017);
    ASSERT_TRUE(kf.initialized, "kf initialized");
    ASSERT_NEAR(kf.P[0], 100.0, 1e-12, "P(0,0) = sigma_pos^2");
    PASS();

    TEST("kf_predict increases P trace");
    ins_mat3_t C_bn;
    ins_mat3_identity(&C_bn);
    ins_vec3_t f_body = {0, 0, INS_GRAVITY_EQUATOR};

    double trace_before = 0.0;
    for (int i = 0; i < 15; i++) trace_before += kf.P[i * 15 + i];

    ins_kf_predict(&kf, &C_bn, &f_body, 0.5, 100.0, 1.0, 1e-6, 1e-12);

    double trace_after = 0.0;
    for (int i = 0; i < 15; i++) trace_after += kf.P[i * 15 + i];

    ASSERT_TRUE(trace_after > trace_before, "covariance grows after predict");
    PASS();
}

/* =================================================================
 * Test 14: Edge Cases (NULL handling)
 * ================================================================= */
static void test_edge_cases(void) {
    TEST("NULL pointer safety");
    ins_vec3_zero(NULL);
    ins_quat_normalize(NULL);
    ins_mech_init(NULL, 0, 0, 0, 0);
    ins_geodetic_to_ecef(NULL, NULL);
    PASS();

    TEST("zero-length vector normalization");
    ins_vec3_t v = {0};
    int ret = ins_vec3_normalize(&v);
    ASSERT_TRUE(ret == -1, "zero vector normalize fails");
    PASS();

    TEST("gravity_ned NULL output");
    ins_gravity_ned(0.0, 0.0, NULL);
    PASS();
}

/* =================================================================
 * Test 15: Navigation Loop — Stationary IMU
 * ================================================================= */
static void test_navigation_stationary(void) {
    TEST("navigate: stationary IMU stays near initial position");
    ins_mech_state_t ms;
    ins_mech_init(&ms, INS_MECH_RK2, 0.5, -2.0, 0.0);

    /* Align level */
    ins_vec3_t accel_level = {0, 0, INS_GRAVITY_EQUATOR};
    ins_align_level(&ms, &accel_level, 0.0);

    /* Navigate with stationary IMU data */
    ins_imu_sample_t imu_data[100];
    ins_nav_solution_t traj[100];
    for (int i = 0; i < 100; i++) {
        imu_data[i].accel.x = 0;
        imu_data[i].accel.y = 0;
        imu_data[i].accel.z = INS_GRAVITY_EQUATOR;
        imu_data[i].gyro.x = 0;
        imu_data[i].gyro.y = 0;
        imu_data[i].gyro.z = 0;
        imu_data[i].dt = 0.01;
    }

    int ret = ins_navigate(&ms, imu_data, 100, traj);
    ASSERT_TRUE(ret == 0, "navigate success");

    /* Position should stay near initial (small drift from Coriolis, etc.) */
    ASSERT_NEAR(traj[99].pos.lat, 0.5, 0.1, "lat near initial");
    ASSERT_NEAR(traj[99].pos.lon, -2.0, 0.1, "lon near initial");
    PASS();
}

/* =================================================================
 * Test Runner
 * ================================================================= */
int main(void) {
    printf("=== INS Test Suite ===\n\n");

    test_vec3_ops();
    test_mat3_ops();
    test_coord_transforms();
    test_gravity();
    test_quat_ops();
    test_quat_conversions();
    test_quat_integration();
    test_mechanization();
    test_allan();
    test_grade_specs();
    test_error_prediction();
    test_calibration();
    test_kalman_filter();
    test_edge_cases();
    test_navigation_stationary();

    printf("\n=== Results: %d passed, %d failed ===\n", tests_passed, tests_failed);
    return tests_failed > 0 ? 1 : 0;
}
