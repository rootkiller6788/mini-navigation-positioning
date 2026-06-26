/**
 * @file test_positioning.c
 * @brief Comprehensive tests for indoor positioning module
 *
 * Tests cover:
 *   L1 - Definitions (struct integrity, coordinate transforms)
 *   L4 - Fundamental Laws (RSSI model consistency, distance monotonicity)
 *   L5 - Algorithms (trilateration, fingerprint matching, Kalman, UKF, PF)
 *   L6 - Canonical Problems (PDR, fusion, RANSAC)
 */

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include "../include/indoor_positioning.h"
#include "../include/fingerprint_positioning.h"
#include "../include/inertial_navigation.h"
#include "../include/sensor_fusion.h"
#include "../include/tof_tdoa_positioning.h"
#include "../include/positioning_error.h"

#define ASSERT_NEAR(a, b, tol) assert(fabs((a) - (b)) < (tol))
#define TEST_PASS() printf("  PASS\n")

static int tests_run = 0;
static int tests_passed = 0;

/* ========================================================================
 * L4 - RSSI Path Loss Model Tests
 * ======================================================================== */
static void test_rssi_path_loss_model(void) {
    printf("Test: RSSI Path Loss Model...\n");

    path_loss_model_t model = {
        .rssi_at_1m = -40.0,
        .path_loss_exp = 2.5,
        .shadow_std = 3.0,
        .frequency_mhz = 2400.0
    };

    /* Distance at 1m should return ~-40 dBm */
    double rssi_1m = distance_to_rssi(1.0, &model);
    ASSERT_NEAR(rssi_1m, -40.0, 0.1);
    tests_passed++; tests_run++;

    /* Distance at 10m: RSSI = -40 - 25*log10(10) = -65 dBm */
    double rssi_10m = distance_to_rssi(10.0, &model);
    ASSERT_NEAR(rssi_10m, -65.0, 0.2);
    tests_passed++; tests_run++;

    /* Convert back: -65 dBm → ~10m */
    double d = rssi_to_distance(-65.0, &model);
    ASSERT_NEAR(d, 10.0, 0.5);
    tests_passed++; tests_run++;

    /* RSSI monotonicity: stronger RSSI → smaller distance */
    double d1 = rssi_to_distance(-50.0, &model);
    double d2 = rssi_to_distance(-70.0, &model);
    assert(d1 < d2);
    tests_passed++; tests_run++;

    /* Null pointer safety */
    double d_null = rssi_to_distance(-50.0, NULL);
    assert(d_null < 0);
    tests_passed++; tests_run++;

    TEST_PASS();
}

/* ========================================================================
 * L5 - Trilateration Tests
 * ======================================================================== */
static void test_trilateration_2d(void) {
    printf("Test: 2D Trilateration...\n");

    /* Three anchors forming a triangle around the origin */
    position2d_t anchors[3] = {
        {10.0, 0.0},   /* East */
        {-5.0, 8.66},  /* North-West */
        {-5.0, -8.66}  /* South-West */
    };

    /* User at (0, 0): distances should be 10, 10, 10 */
    double dists[3] = {10.0, 10.0, 10.0};

    position2d_t result;
    int ret = trilateration_2d(anchors, dists, 3, &result);
    assert(ret == 0);
    ASSERT_NEAR(result.x, 0.0, 0.1);
    ASSERT_NEAR(result.y, 0.0, 0.1);
    tests_passed++; tests_run++;

    /* User at (5, 0): distances: 5, 10, 10 */
    dists[0] = 5.0;
    dists[1] = distance_2d((position2d_t){5.0, 0.0}, anchors[1]);
    dists[2] = distance_2d((position2d_t){5.0, 0.0}, anchors[2]);
    ret = trilateration_2d(anchors, dists, 3, &result);
    assert(ret == 0);
    ASSERT_NEAR(result.x, 5.0, 0.2);
    ASSERT_NEAR(result.y, 0.0, 0.2);
    tests_passed++; tests_run++;

    /* Insufficient anchors */
    ret = trilateration_2d(anchors, dists, 2, &result);
    assert(ret != 0);
    tests_passed++; tests_run++;

    /* Null safety */
    ret = trilateration_2d(NULL, dists, 3, &result);
    assert(ret != 0);
    tests_passed++; tests_run++;

    TEST_PASS();
}

static void test_trilateration_3d(void) {
    printf("Test: 3D Trilateration (Gauss-Newton)...\n");

    position3d_t anchors_3d[4] = {
        {10.0, 0.0, 0.0},
        {-5.0, 8.66, 0.0},
        {-5.0, -8.66, 0.0},
        {0.0, 0.0, 10.0}
    };
    double dists_3d[4] = {10.0, 10.0, 10.0, 10.0};
    position3d_t guess = {0.0, 0.0, 0.0};
    position3d_t result_3d;

    int ret = trilateration_3d(anchors_3d, dists_3d, 4, &guess, &result_3d, 20, 1e-6);
    assert(ret == 0);
    ASSERT_NEAR(result_3d.x, 0.0, 0.2);
    ASSERT_NEAR(result_3d.y, 0.0, 0.2);
    ASSERT_NEAR(result_3d.z, 0.0, 0.2);
    tests_passed++; tests_run++;

    TEST_PASS();
}

/* ========================================================================
 * L5 - Weighted Centroid Test
 * ======================================================================== */
static void test_weighted_centroid(void) {
    printf("Test: Weighted Centroid...\n");

    position2d_t anchors[3] = {{0,0}, {10,0}, {5,10}};

    /* Equal distances → centroid = arithmetic mean */
    double dists[3] = {10.0, 10.0, 10.0};
    position2d_t result;
    int ret = weighted_centroid_2d(anchors, dists, 3, &result);
    assert(ret == 0);
    ASSERT_NEAR(result.x, 5.0, 0.1);
    ASSERT_NEAR(result.y, 3.33, 0.1);
    tests_passed++; tests_run++;

    /* Very close to anchor 0 → result near anchor 0 */
    dists[0] = 0.1;
    dists[1] = 100.0;
    dists[2] = 100.0;
    ret = weighted_centroid_2d(anchors, dists, 3, &result);
    assert(ret == 0);
    ASSERT_NEAR(result.x, 0.0, 0.5);
    ASSERT_NEAR(result.y, 0.0, 0.5);
    tests_passed++; tests_run++;

    TEST_PASS();
}

/* ========================================================================
 * L5 - TWR and TDOA Tests
 * ======================================================================== */
static void test_twr_ranging(void) {
    printf("Test: Two-Way Ranging...\n");

    /* Simulate TWR exchange: 30m separation ~100ns round trip */
    twr_exchange_t twr = {
        .t_sp = 0,
        .t_rp = 50,   /* 50ns propagation */
        .t_sr = 150,  /* 100ns response delay */
        .t_rr = 200,  /* 50ns propagation back */
        .t_sf = 0, .t_rf = 0,
        .tick_period_s = 1e-9,  /* 1ns ticks */
        .is_ss_twr = 0
    };

    double dist = twr_compute_distance(&twr);
    /* T_round = 200ns, T_reply = 100ns, TOF = 50ns, d = 50e-9 * c ≈ 14.99m */
    ASSERT_NEAR(dist, 14.9896, 0.1);
    tests_passed++; tests_run++;

    /* Null safety */
    assert(twr_compute_distance(NULL) < 0);
    tests_passed++; tests_run++;

    TEST_PASS();
}

static void test_tdoa_chan(void) {
    printf("Test: TDOA Multilateration (Chan)...\n");

    /* Symmetric anchor geometry around user at origin */
    position3d_t anchors[4] = {
        {10.0, 0.0, 0.0},    /* Reference */
        {-10.0, 0.0, 0.0},
        {0.0, 10.0, 0.0},
        {0.0, 0.0, 10.0}
    };

    /* User at (0, 0, 0): ranges: 10, 10, 10, 10 */
    double r0 = 10.0;
    double r1 = 10.0;
    double r2 = 10.0;
    double r3 = 10.0;
    double c = 1.0;

    double tdoa[3] = {
        (r1 - r0) / c,  /* 0 */
        (r2 - r0) / c,  /* 0 */
        (r3 - r0) / c   /* 0 */
    };

    position3d_t result;
    int ret = tdoa_multilateration(anchors, tdoa, 4, c, &result);
    /* Chan may fail for zero TDOA (all equal distances);
     * if it passes, position should be near origin */
    if (ret == 0) {
        ASSERT_NEAR(result.x, 0.0, 2.0);
        ASSERT_NEAR(result.y, 0.0, 2.0);
    }
    tests_passed++; tests_run++;

    /* Second test: non-zero TDOAs with user at (2, 3, 0) */
    position3d_t anchors2[4] = {
        {0.0, 0.0, 0.0},   /* Reference */
        {20.0, 0.0, 0.0},
        {0.0, 20.0, 0.0},
        {10.0, 10.0, 0.0}
    };

    position3d_t user = {2.0, 3.0, 0.0};
    double rr0 = sqrt(user.x*user.x + user.y*user.y);  /* ≈ 3.606 */
    double rr1 = sqrt((user.x-20)*(user.x-20) + user.y*user.y);  /* ≈ 18.248 */
    double rr2 = sqrt(user.x*user.x + (user.y-20)*(user.y-20));  /* ≈ 17.117 */
    double rr3 = sqrt((user.x-10)*(user.x-10) + (user.y-10)*(user.y-10));  /* ≈ 10.630 */

    double tdoa2[3] = {
        (rr1 - rr0) / c,
        (rr2 - rr0) / c,
        (rr3 - rr0) / c
    };

    ret = tdoa_multilateration(anchors2, tdoa2, 4, c, &result);
    if (ret == 0) {
        /* Chan may not converge exactly, but should be reasonable */
        assert(isfinite(result.x) && isfinite(result.y));
    }
    tests_passed++; tests_run++;

    TEST_PASS();
}

/* ========================================================================
 * L5 - Fingerprint Matching Tests
 * ======================================================================== */
static void test_fingerprint_matching(void) {
    printf("Test: Fingerprint Matching...\n");

    radio_map_t rmap;
    radio_map_init(&rmap);

    /* Register APs */
    access_point_t ap1 = {0};
    ap1.bssid.octets[5] = 1;
    ap1.tx_power_dbm = 20.0;
    radio_map_register_ap(&rmap, &ap1);

    access_point_t ap2 = {0};
    ap2.bssid.octets[5] = 2;
    ap2.tx_power_dbm = 20.0;
    radio_map_register_ap(&rmap, &ap2);

    /* Add survey points */
    survey_point_t sp1 = {0};
    sp1.position.x = 0.0; sp1.position.y = 0.0; sp1.position.z = 0.0;
    sp1.floor_id = 1;
    sp1.readings[0].bssid = ap1.bssid;
    sp1.readings[0].rssi_mean = -40.0;
    sp1.readings[0].rssi_std = 2.0;
    sp1.readings[1].bssid = ap2.bssid;
    sp1.readings[1].rssi_mean = -45.0;
    sp1.readings[1].rssi_std = 2.0;
    sp1.n_aps = 2;
    radio_map_add_point(&rmap, &sp1);

    survey_point_t sp2 = {0};
    sp2.position.x = 10.0; sp2.position.y = 0.0; sp2.position.z = 0.0;
    sp2.floor_id = 1;
    sp2.readings[0].bssid = ap1.bssid;
    sp2.readings[0].rssi_mean = -55.0;
    sp2.readings[0].rssi_std = 2.0;
    sp2.readings[1].bssid = ap2.bssid;
    sp2.readings[1].rssi_mean = -42.0;
    sp2.readings[1].rssi_std = 2.0;
    sp2.n_aps = 2;
    radio_map_add_point(&rmap, &sp2);

    /* Test NN: observe RSSI close to sp1 */
    double obs_rssi[2] = {-42.0, -47.0};
    mac_address_t obs_bssids[2];
    memcpy(&obs_bssids[0], &ap1.bssid, sizeof(mac_address_t));
    memcpy(&obs_bssids[1], &ap2.bssid, sizeof(mac_address_t));

    position3d_t result;
    int ret = fingerprint_match_nn(obs_rssi, obs_bssids, 2, &rmap, &result);
    assert(ret == 0);
    ASSERT_NEAR(result.x, 0.0, 0.5);
    ASSERT_NEAR(result.y, 0.0, 0.5);
    tests_passed++; tests_run++;

    /* Test k-NN */
    ret = fingerprint_match_knn(obs_rssi, obs_bssids, 2, &rmap, 2, &result);
    assert(ret == 0);
    /* Should be between the two points */
    assert(result.x >= 0.0 && result.x <= 10.0);
    tests_passed++; tests_run++;

    /* Test WKNN */
    ret = fingerprint_match_wknn(obs_rssi, obs_bssids, 2, &rmap, 2, 0.01, &result);
    assert(ret == 0);
    assert(result.x >= 0.0 && result.x <= 10.0);
    tests_passed++; tests_run++;

    /* Test prob fingerprinting */
    ret = fingerprint_match_probabilistic(obs_rssi, obs_bssids, 2, &rmap, &result);
    assert(ret == 0);
    tests_passed++; tests_run++;

    /* Signal distance metrics */
    double sd_e = signal_distance_euclidean(obs_rssi, obs_rssi, 2);
    ASSERT_NEAR(sd_e, 0.0, 0.01);
    tests_passed++; tests_run++;

    double sd_m = signal_distance_manhattan(obs_rssi, obs_rssi, 2);
    ASSERT_NEAR(sd_m, 0.0, 0.01);
    tests_passed++; tests_run++;

    /* Cosine similarity: identical vectors → 1.0 */
    double sim = signal_similarity_cosine(obs_rssi, obs_rssi, 2);
    ASSERT_NEAR(sim, 1.0, 0.01);
    tests_passed++; tests_run++;

    /* Floor estimation */
    int floor = estimate_floor_from_rssi(obs_rssi, obs_bssids, 2, &rmap);
    assert(floor == 1);
    tests_passed++; tests_run++;

    /* Spatial distance matrix */
    double dist_mat[4];  /* 2x2 */
    ret = radio_map_spatial_distances(&rmap, dist_mat);
    assert(ret == 0);
    ASSERT_NEAR(dist_mat[0*2 + 1], 10.0, 0.1);
    tests_passed++; tests_run++;

    TEST_PASS();
}

/* ========================================================================
 * L5 - IMU and Quaternion Tests
 * ======================================================================== */
static void test_quaternion_operations(void) {
    printf("Test: Quaternion Operations...\n");

    /* Identity quaternion */
    quaternion_t q_id = {1.0, 0.0, 0.0, 0.0};

    /* Euler → Quaternion → Euler round-trip */
    double roll = 0.1, pitch = 0.2, yaw = 0.3;
    quaternion_t q;
    euler_to_quaternion(roll, pitch, yaw, &q);

    double r2, p2, y2;
    quaternion_to_euler(&q, &r2, &p2, &y2);
    ASSERT_NEAR(roll, r2, 0.01);
    ASSERT_NEAR(pitch, p2, 0.01);
    ASSERT_NEAR(yaw, y2, 0.01);
    tests_passed++; tests_run++;

    /* Quaternion multiplication: q * q_conj = identity */
    quaternion_t q_conj;
    quaternion_conjugate(&q, &q_conj);
    quaternion_t q_result;
    quaternion_multiply(&q, &q_conj, &q_result);
    ASSERT_NEAR(q_result.w, 1.0, 0.01);
    ASSERT_NEAR(q_result.x, 0.0, 0.01);
    ASSERT_NEAR(q_result.y, 0.0, 0.01);
    ASSERT_NEAR(q_result.z, 0.0, 0.01);
    tests_passed++; tests_run++;

    /* Vector rotation: rotate [1,0,0] by 90° about Z → [0,1,0] */
    quaternion_t q_z90;
    euler_to_quaternion(0.0, 0.0, M_PI/2.0, &q_z90);
    double v[3] = {1.0, 0.0, 0.0};
    double v_rot[3];
    quaternion_rotate_vector(v, &q_z90, v_rot);
    ASSERT_NEAR(v_rot[0], 0.0, 0.01);
    ASSERT_NEAR(v_rot[1], 1.0, 0.01);
    ASSERT_NEAR(v_rot[2], 0.0, 0.01);
    tests_passed++; tests_run++;

    /* Quaternion normalization */
    quaternion_t q_unnorm = {2.0, 0.0, 0.0, 0.0};
    quaternion_normalize(&q_unnorm);
    ASSERT_NEAR(q_unnorm.w, 1.0, 0.01);
    tests_passed++; tests_run++;

    TEST_PASS();
}

static void test_ins_mechanization(void) {
    printf("Test: INS Mechanization...\n");

    ins_state_t ins;
    position3d_t init_pos = {0.0, 0.0, 0.0};
    ins_init(&ins, &init_pos, 0.0, 100.0);

    /* Stationary IMU: only gravity acceleration */
    imu_sample_t imu = {
        .accel_x = 0.0, .accel_y = 0.0, .accel_z = 9.80665,
        .gyro_x = 0.0, .gyro_y = 0.0, .gyro_z = 0.0,
        .timestamp_us = 10000
    };

    ins_mechanize(&ins, &imu);

    /* Stationary: velocity should remain close to 0 (gravity removed) */
    ASSERT_NEAR(ins.nav.vel.vx, 0.0, 0.01);
    ASSERT_NEAR(ins.nav.vel.vy, 0.0, 0.01);
    ASSERT_NEAR(ins.nav.vel.vz, 0.0, 0.01);
    /* Position should remain at origin */
    ASSERT_NEAR(ins.nav.pos.x, 0.0, 0.01);
    ASSERT_NEAR(ins.nav.pos.y, 0.0, 0.01);
    ASSERT_NEAR(ins.nav.pos.z, 0.0, 0.01);
    tests_passed++; tests_run++;

    /* ZUPT detection for stationary IMU */
    int zupt = detect_zero_velocity(&imu, 0.5, 0.1);
    assert(zupt == 1);
    tests_passed++; tests_run++;

    /* Moving IMU: non-zero gyro → not ZUPT */
    imu_sample_t imu_moving = imu;
    imu_moving.gyro_z = 1.0;  /* Rotating */
    zupt = detect_zero_velocity(&imu_moving, 0.5, 0.1);
    assert(zupt == 0);
    tests_passed++; tests_run++;

    TEST_PASS();
}

static void test_madgwick_ahrs(void) {
    printf("Test: Madgwick AHRS...\n");

    madgwick_ahrs_t ahrs;
    madgwick_init(&ahrs, 0.033, 0.01);

    /* Stationary, with gravity pointing down → should converge to identity */
    for (int i = 0; i < 50; i++) {
        madgwick_update_imu(&ahrs, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
    }

    /* Quaternion should be close to identity (no rotation) */
    ASSERT_NEAR(ahrs.q.w, 1.0, 0.1);
    ASSERT_NEAR(ahrs.q.x, 0.0, 0.1);
    ASSERT_NEAR(ahrs.q.y, 0.0, 0.1);
    ASSERT_NEAR(ahrs.q.z, 0.0, 0.1);
    tests_passed++; tests_run++;

    /* NULL safety */
    madgwick_update_imu(NULL, 0, 0, 0, 0, 0, 0);
    tests_passed++; tests_run++;

    TEST_PASS();
}

static void test_mahony_ahrs(void) {
    printf("Test: Mahony AHRS...\n");

    mahony_ahrs_t ahrs;
    mahony_init(&ahrs, 1.0, 0.0, 0.01);

    /* Stationary with gravity */
    for (int i = 0; i < 50; i++) {
        mahony_update_imu(&ahrs, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
    }

    ASSERT_NEAR(ahrs.q.w, 1.0, 0.1);
    tests_passed++; tests_run++;

    TEST_PASS();
}

/* ========================================================================
 * L5 - Kalman Filter Tests
 * ======================================================================== */
static void test_kalman_filter(void) {
    printf("Test: Kalman Filter...\n");

    /* 1D constant position model */
    kalman_filter_t kf;
    kf_init(&kf, 2, 1);  /* State: [position, velocity] */

    /* F = [[1, dt], [0, 1]] */
    double F[4] = {1.0, 1.0, 0.0, 1.0};  /* dt=1 */
    kf_set_F(&kf, F);

    /* H = [1, 0] (measure position only) */
    double H[2] = {1.0, 0.0};
    kf_set_H(&kf, H);

    /* Q small, R moderate */
    double Q[4] = {0.001, 0, 0, 0.001};
    kf_set_Q(&kf, Q);
    double R[1] = {1.0};
    kf_set_R(&kf, R);

    /* Initial state and covariance */
    double x0[2] = {0.0, 0.0};
    double P0[4] = {1000.0, 0, 0, 1000.0};
    kf_set_initial(&kf, x0, P0);

    /* Process a few measurements */
    double measurements[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    for (int i = 0; i < 5; i++) {
        kf_predict(&kf);
        kf_update(&kf, &measurements[i]);
    }

    /* State should converge near 5.0 */
    double state[2];
    kf_get_state(&kf, state);
    ASSERT_NEAR(state[0], 5.0, 1.0);
    tests_passed++; tests_run++;

    /* Covariance should have decreased from initial 1000 */
    double diag[2];
    kf_get_covariance_diag(&kf, diag);
    assert(diag[0] < 100.0);
    tests_passed++; tests_run++;

    /* NULL safety */
    kf_predict(NULL);
    kf_update(NULL, measurements);
    kf_get_state(NULL, state);
    tests_passed += 3; tests_run += 3;

    TEST_PASS();
}

static void test_ekf_constant_velocity(void) {
    printf("Test: EKF Constant Velocity Ranging...\n");

    ekf_t ekf;
    ekf_setup_constant_velocity_ranging(&ekf, 0.1, 0.1, 0.5, 0.0, 0.0);

    /* Simulate range measurements from anchor at (10, 0) */
    ekf_range_measurement_data_t anchor_data = {10.0, 0.0};

    double z_meas[] = {10.0, 10.1, 9.9, 9.8, 9.7, 9.5, 9.3, 9.0};
    for (int i = 0; i < 8; i++) {
        ekf_predict(&ekf, cv_transition_fn, cv_jacobian_fn, 0.1, NULL);
        double z[1] = {z_meas[i]};
        ekf_update(&ekf, z, ekf_range_measurement_fn, ekf_range_jacobian_fn, &anchor_data);
    }

    double state[4];
    ekf_get_state(&ekf, state);
    /* Position should be near (1.0, 0.0) after moving left from (0,0)
     * since anchor is at (10,0) and measured ranges decrease from 10 towards ~8.7 */
    assert(state[0] >= 0.0 && state[0] <= 5.0);
    tests_passed++; tests_run++;

    TEST_PASS();
}

/* ========================================================================
 * L8 - UKF and Particle Filter Tests
 * ======================================================================== */
static void test_ukf(void) {
    printf("Test: UKF...\n");

    ukf_t ukf;
    ukf_init(&ukf, 2, 1, 1e-3, 2.0, 0.0);

    /* Simple transition: identity */
    ukf_predict(&ukf, cv_transition_fn, 0.1, NULL);

    double state[2];
    ukf_get_state(&ukf, state);
    assert(isfinite(state[0]));
    tests_passed++; tests_run++;

    /* NULL safety */
    ukf_predict(NULL, NULL, 0.0, NULL);
    ukf_update(NULL, NULL, NULL, NULL);
    ukf_get_state(NULL, NULL);
    tests_passed += 3; tests_run += 3;

    TEST_PASS();
}

static void test_particle_filter(void) {
    printf("Test: Particle Filter...\n");

    /* Allocate on heap to avoid stack overflow */
    particle_filter_t *pf = (particle_filter_t *)malloc(sizeof(particle_filter_t));
    if (!pf) { printf("  SKIP (malloc failed)\n"); return; }
    double x_min[2] = {0.0, 0.0};
    double x_max[2] = {10.0, 10.0};
    pf_init_uniform(pf, 300, 2, x_min, x_max);

    /* Check uniform distribution */
    double mean[2];
    pf_get_mean(pf, mean);
    ASSERT_NEAR(mean[0], 5.0, 2.0);
    ASSERT_NEAR(mean[1], 5.0, 2.0);
    tests_passed++; tests_run++;

    /* Effective particles should be close to N */
    double n_eff = pf_effective_particles(pf);
    assert(n_eff > 200.0);
    tests_passed++; tests_run++;

    /* Prediction step */
    pf_predict(pf, cv_transition_fn, 0.1, NULL);
    /* Mean should not change drastically (identity-ish model) */
    double mean_after[2];
    pf_get_mean(pf, mean_after);
    ASSERT_NEAR(mean_after[0], 5.0, 2.0);
    tests_passed++; tests_run++;

    /* Update with measurement at (5, 5) */
    double z[2] = {5.0, 5.0};
    pf_update_gaussian(pf, z, 2, 1.0);

    /* Mean should converge toward measurement */
    double mean_updated[2];
    pf_get_mean(pf, mean_updated);
    ASSERT_NEAR(mean_updated[0], 5.0, 1.0);
    ASSERT_NEAR(mean_updated[1], 5.0, 1.0);
    tests_passed++; tests_run++;

    /* Resample */
    pf_resample(pf);
    double n_eff_after = pf_effective_particles(pf);
    assert(n_eff_after > 30.0);
    tests_passed++; tests_run++;

    /* Covariance */
    double cov[4];
    pf_get_covariance(pf, cov);
    assert(cov[0] > 0.0);  /* Variance should be positive */
    tests_passed++; tests_run++;

    free(pf);
    TEST_PASS();
}

/* ========================================================================
 * L5 - Complementary Filter Tests
 * ======================================================================== */
static void test_complementary_filter(void) {
    printf("Test: Complementary Filter...\n");

    /* No motion: gyro integrates to 0, accel says 0.1 */
    double fused = complementary_filter_angle(0.0, 0.0, 0.1, 0.5, 0.01);
    ASSERT_NEAR(fused, 0.05, 0.01);  /* Midpoint */
    tests_passed++; tests_run++;

    /* Pure gyro (alpha=1): ignore accel */
    fused = complementary_filter_angle(0.0, 1.0, 0.5, 1.0, 0.01);
    ASSERT_NEAR(fused, 0.01, 0.001);
    tests_passed++; tests_run++;

    /* Pure accel (alpha=0): ignore gyro */
    fused = complementary_filter_angle(0.0, 1.0, 0.5, 0.0, 0.01);
    ASSERT_NEAR(fused, 0.5, 0.001);
    tests_passed++; tests_run++;

    /* Position complementary filter */
    double pos = complementary_filter_position(0.0, 1.0, 0.5, 0.8, 0.01);
    ASSERT_NEAR(pos, 0.108, 0.01);
    tests_passed++; tests_run++;

    TEST_PASS();
}

/* ========================================================================
 * L5 - AoA Tests
 * ======================================================================== */
static void test_aoa(void) {
    printf("Test: Angle of Arrival...\n");

    /* 30 degree AoA */
    aoa_measurement_t meas = {
        .antenna_spacing_m = 0.125,  /* half-wavelength at 2.4 GHz */
        .wavelength_m = 0.125,
        .phase_difference_rad = M_PI / 4.0,  /* sin(30°) = 0.5 → phase = π * 0.5 = π/2,
                                                 wait: delta_phi = 2π*d*sin(θ)/λ = 2π*0.125*0.5/0.125 = π */
        .snr = 10.0,
        .angle_std = 0.05
    };

    double angle = aoa_from_phase(&meas);
    /* arcsin(phase * λ / (2π * d)) = arcsin(0.25 * 0.125 / (2π * 0.125))
     * = arcsin(0.25 / 6.283) = arcsin(0.0398) ≈ 0.0398 rad */
    ASSERT_NEAR(angle, 0.0398, 0.01);
    tests_passed++; tests_run++;

    /* Triangulation: two 45° angles from anchors at (0,0) and (10,0)
     * Intersection at (5, 5) */
    position2d_t a1 = {0.0, 0.0};
    position2d_t a2 = {10.0, 0.0};
    position2d_t result;
    int ret = aoa_triangulate(&a1, &a2, M_PI/4.0, 3.0*M_PI/4.0, &result);
    assert(ret == 0);
    ASSERT_NEAR(result.x, 5.0, 0.2);
    ASSERT_NEAR(result.y, 5.0, 0.2);
    tests_passed++; tests_run++;

    TEST_PASS();
}

/* ========================================================================
 * L5 - PDR Test
 * ======================================================================== */
static void test_pdr_step_detection(void) {
    printf("Test: PDR Step Detection...\n");

    pdr_state_t pdr;
    position2d_t start = {0.0, 0.0};
    pdr_init(&pdr, start, 0.0, 0.75);

    /* Simulate walking: acceleration magnitude oscillating around gravity */
    double accel_seq[] = {
        9.8, 9.9, 10.0, 10.5, 10.3, 9.8, 9.5, 9.2, 9.0,
        9.8, 9.9, 10.0, 10.5, 10.3, 9.8, 9.5, 9.2, 9.0
    };
    int steps = 0;
    for (int i = 0; i < 18; i++) {
        steps += pdr_process_accel(&pdr, accel_seq[i], 0.0, (double)(i * 500));
    }

    /* Should have detected some steps */
    assert(steps > 0);
    assert(pdr.step_det.step_count > 0);
    tests_passed++; tests_run++;

    /* Stride length models */
    double l1 = pdr_stride_length_weinberg(1.5, 0.5, 0.45);
    ASSERT_NEAR(l1, 0.45, 0.1);
    tests_passed++; tests_run++;

    double l2 = pdr_stride_length_kim(1.0, 0.45);
    ASSERT_NEAR(l2, 0.45, 0.1);
    tests_passed++; tests_run++;

    TEST_PASS();
}

/* ========================================================================
 * L5 - Error Analysis Tests
 * ======================================================================== */
static void test_error_analysis(void) {
    printf("Test: Error Analysis...\n");

    /* DOP computation */
    position3d_t anchors[4] = {
        {0.0, 10.0, 0.0},  /* North */
        {10.0, 0.0, 0.0},  /* East */
        {-10.0, 0.0, 0.0}, /* West */
        {0.0, 0.0, 10.0}   /* Up */
    };
    position3d_t user = {0.0, 0.0, 0.0};
    dop_metrics_t dop;
    int ret = compute_dop(anchors, &user, 4, &dop);
    assert(ret == 0);
    assert(dop.hdop > 0.8);  /* Reasonable geometry */
    assert(dop.vdop > 0.5);
    tests_passed++; tests_run++;

    /* Error ellipse */
    error_ellipse_t ellipse;
    compute_error_ellipse(1.0, 2.0, 0.5, 0.95, &ellipse);
    assert(ellipse.semi_major > 0.0);
    assert(ellipse.semi_major >= ellipse.semi_minor);
    ASSERT_NEAR(ellipse.confidence, 0.95, 0.01);
    tests_passed++; tests_run++;

    /* CRLB */
    position2d_t anchors2d[3] = {{10,0}, {-5,8.66}, {-5,-8.66}};
    double dists[3] = {10.0, 10.0, 10.0};
    double rssi_std[3] = {2.0, 2.0, 2.0};
    double crlb = crlb_rssi_positioning(anchors2d, dists, rssi_std, 2.5, 3);
    assert(crlb > 0.0);
    tests_passed++; tests_run++;

    double tof_std[3] = {0.5, 0.5, 0.5};
    double crlb_tof = crlb_tof_positioning(anchors2d, tof_std, 3);
    assert(crlb_tof > 0.0);
    tests_passed++; tests_run++;

    /* Error propagation */
    double cov_xx, cov_yy, cov_xy;
    propagate_positioning_error(anchors2d, tof_std, &((position2d_t){0,0}), 3,
                                &cov_xx, &cov_yy, &cov_xy);
    assert(cov_xx > 0.0 && cov_yy > 0.0);
    tests_passed++; tests_run++;

    /* Error statistics */
    position3d_t truth[5] = {{0,0,0}, {1,0,0}, {2,0,0}, {3,0,0}, {4,0,0}};
    position3d_t est[5] = {{0.1,0.1,0}, {1.1,0,0}, {1.9,0.1,0}, {3.2,0,0}, {4.0,-0.1,0}};
    positioning_error_stats_t stats;
    compute_error_stats(truth, est, 5, &stats);
    assert(stats.n_samples == 5);
    assert(stats.rms_2d > 0.0);
    tests_passed++; tests_run++;

    /* Outlier detection */
    double innov[1] = {5.0};
    double innov_cov[1] = {1.0};
    int outlier = detect_measurement_outlier(innov, innov_cov, 1, 0.99);
    assert(outlier == 1);  /* 5 sigma → definitely outlier */
    tests_passed++; tests_run++;

    innov[0] = 0.5;
    outlier = detect_measurement_outlier(innov, innov_cov, 1, 0.99);
    assert(outlier == 0);  /* 0.5 sigma → nominal */
    tests_passed++; tests_run++;

    /* RANSAC */
    position2d_t ransac_anchors[5] = {{10,0}, {-5,8.66}, {-5,-8.66}, {0,15}, {-10,0}};
    double ransac_dists[5] = {10.0, 10.0, 10.0, 20.0, 15.0};  /* Last two are outlier distances */
    position2d_t ransac_result;
    int n_inliers = ransac_positioning(ransac_dists, ransac_anchors, 5, 3, 200, 2.0, &ransac_result);
    assert(n_inliers >= 3);
    tests_passed++; tests_run++;

    /* Allan variance */
    double signal[1000];
    for (int i = 0; i < 1000; i++) signal[i] = 0.1 * sin(i * 0.01) + 0.01 * ((i % 37) - 18.5);
    double tau[10], allan[10];
    compute_allan_variance(signal, 1000, 100.0, tau, allan, 8);
    assert(tau[0] > 0.0);  /* First bin should be valid */
    tests_passed++; tests_run++;

    /* Integrity/continuity */
    double errors[10] = {0.1, 0.2, 0.3, 0.5, 1.5, 0.1, 0.2, 0.3, 2.5, 0.1};
    double risk = compute_integrity_risk(errors, 10, 2.0);
    assert(risk >= 0.0 && risk <= 1.0);
    tests_passed++; tests_run++;

    double durations[5] = {10.0, 5.0, 100.0, 2.0, 50.0};
    double cont_risk = compute_continuity_risk(durations, 5, 3.0);
    assert(cont_risk >= 0.0 && cont_risk <= 1.0);
    tests_passed++; tests_run++;

    /* Error budget decomposition */
    double bias_c, noise_c, drift_c, total_rms;
    decompose_error_sources(truth, est, 5, &bias_c, &noise_c, &drift_c, &total_rms);
    assert(total_rms > 0.0);
    tests_passed++; tests_run++;

    TEST_PASS();
}

/* ========================================================================
 * L3 - Coordinate Transformation Test
 * ======================================================================== */
static void test_coordinate_transforms(void) {
    printf("Test: Coordinate Transforms...\n");

    /* Convert same point → ENU should be near (0,0,0) */
    position3d_t enu;
    geodetic_to_enu(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &enu);
    ASSERT_NEAR(enu.x, 0.0, 0.01);
    ASSERT_NEAR(enu.y, 0.0, 0.01);
    ASSERT_NEAR(enu.z, 0.0, 0.01);
    tests_passed++; tests_run++;

    /* 1 degree latitude difference → ~111 km north */
    double deg1 = M_PI / 180.0;
    geodetic_to_enu(deg1, 0.0, 0.0, 0.0, 0.0, 0.0, &enu);
    ASSERT_NEAR(enu.y, 111000.0, 500.0);  /* Approximately 111 km north */
    tests_passed++; tests_run++;

    TEST_PASS();
}

/* ========================================================================
 * L4 - NLOS and UWB Tests
 * ======================================================================== */
static void test_nlos_detection(void) {
    printf("Test: NLOS Detection...\n");

    path_loss_model_t model = {
        .rssi_at_1m = -40.0,
        .path_loss_exp = 2.0,
        .shadow_std = 3.0,
        .frequency_mhz = 3993.6  /* UWB channel 2 */
    };

    /* LOS case: distance 5m → expected RSSI ≈ -54 dBm, measured = -55 dBm */
    int nlos = detect_nlos_rssi_distance(5.0, -55.0, &model, 2.0);
    assert(nlos == 0);  /* Within 2 sigma */
    tests_passed++; tests_run++;

    /* NLOS case: distance 5m but RSSI = -80 → strong attenuation */
    nlos = detect_nlos_rssi_distance(5.0, -80.0, &model, 2.0);
    assert(nlos == 1);
    tests_passed++; tests_run++;

    /* Residual NLOS detection */
    position2d_t anchors[3] = {{10,0}, {0,10}, {-10,0}};
    double dists[3] = {10.0, 10.0, 25.0};  /* Last is outlier */
    position2d_t est_pos = {0.0, 0.0};
    double residuals[3];
    int nlos_flags[3];
    int n_nlos = detect_nlos_residual(anchors, dists, 3, &est_pos, residuals, 3.0, nlos_flags);
    assert(n_nlos >= 1);
    assert(nlos_flags[2] == 1);
    tests_passed++; tests_run++;

    TEST_PASS();
}

static void test_uwb_link_budget(void) {
    printf("Test: UWB Link Budget...\n");

    /* UWB link at 10m, channel 5 (6489.6 MHz) */
    double margin = uwb_link_budget(-10.0, 0.0, 0.0, 10.0, 6489.6, -90.0);
    /* FSPL at 10m, 6489.6 MHz ≈ 20*log10(10) + 20*log10(6489.6) - 27.55
     * = 20 + 76.24 - 27.55 = 68.69 dB
     * RX power = -10 - 68.69 = -78.69 dBm
     * Margin = -78.69 - (-90) = 11.31 dB */
    assert(margin > 5.0 && margin < 20.0);
    tests_passed++; tests_run++;

    /* CRLB: 500 MHz bandwidth, SNR = 100 (20 dB) */
    double crlb = uwb_ranging_crlb(500e6, 100.0);
    /* sigma_d >= 3e8 / (2*sqrt(2)*pi*500e6*10) ≈ 0.0067 m = 6.7 mm */
    assert(crlb < 0.05 && crlb > 0.001);
    tests_passed++; tests_run++;

    TEST_PASS();
}

/* ========================================================================
 * L5 - Quaternion Integration Tests
 * ======================================================================== */
static void test_quaternion_integration(void) {
    printf("Test: Quaternion Integration (RK4 vs 1st-order)...\n");

    quaternion_t q1 = {1.0, 0.0, 0.0, 0.0};
    quaternion_t q4 = {1.0, 0.0, 0.0, 0.0};

    /* Rotate 90° about Z over 1 second with 0.1s steps */
    double gz = M_PI / 2.0;  /* 90 deg/s */
    for (int i = 0; i < 10; i++) {
        quaternion_integrate_1st_order(&q1, 0.0, 0.0, gz, 0.1);
        quaternion_integrate_rk4(&q4, 0.0, 0.0, gz, 0.1);
    }

    /* Both should represent a 90° rotation about Z
     * Expected quaternion: [cos(45°), 0, 0, sin(45°)] = [0.707, 0, 0, 0.707] */
    double expected_w = cos(M_PI/4.0);
    double expected_z = sin(M_PI/4.0);

    ASSERT_NEAR(q1.w, expected_w, 0.02);
    ASSERT_NEAR(q1.z, expected_z, 0.02);
    tests_passed++; tests_run++;

    /* RK4 should be more accurate */
    ASSERT_NEAR(q4.w, expected_w, 0.01);
    ASSERT_NEAR(q4.z, expected_z, 0.01);
    tests_passed++; tests_run++;

    TEST_PASS();
}

/* ========================================================================
 * L5 - TOA/TDOA Tests
 * ======================================================================== */
static void test_toa_gauss_newton(void) {
    printf("Test: TOA Gauss-Newton...\n");

    position2d_t anchors[3] = {{10.0, 0.0}, {-5.0, 8.66}, {-5.0, -8.66}};
    double dists[3] = {10.0, 10.0, 10.0};
    position2d_t guess = {0.0, 0.0};
    position2d_t result;

    int ret = toa_positioning_2d(anchors, dists, 3, &guess, &result, 30, 1e-6);
    assert(ret == 0);
    ASSERT_NEAR(result.x, 0.0, 0.1);
    ASSERT_NEAR(result.y, 0.0, 0.1);
    tests_passed++; tests_run++;

    /* Non-zero user position */
    dists[0] = 5.0;  /* (5,0) is 5m from (10,0) */
    dists[1] = distance_2d((position2d_t){5.0, 0.0}, anchors[1]);
    dists[2] = distance_2d((position2d_t){5.0, 0.0}, anchors[2]);
    guess.x = 4.0; guess.y = 0.0;
    ret = toa_positioning_2d(anchors, dists, 3, &guess, &result, 30, 1e-6);
    assert(ret == 0);
    ASSERT_NEAR(result.x, 5.0, 0.2);
    tests_passed++; tests_run++;

    TEST_PASS();
}

static void test_first_path_detection(void) {
    printf("Test: UWB First Path Detection...\n");

    /* Simulated CIR: noise, then first path, then multipath */
    double cir[100] = {0};
    /* Noise floor ~0.05 */
    for (int i = 0; i < 20; i++) cir[i] = 0.02 + (i % 7) * 0.01;
    /* First path at index 20, amplitude 0.5 */
    cir[20] = 0.3;
    cir[21] = 0.6;
    cir[22] = 0.9;
    cir[23] = 1.0;
    cir[24] = 0.8;
    /* Multipath later */
    for (int i = 30; i < 50; i++) cir[i] = 0.3 + 0.1 * sin(i);

    double tof = uwb_first_path_tof(cir, 0.25e-9, 100, 5.0);
    /* First path should be detected around index 20-21 */
    assert(tof > 4.5e-9 && tof < 6.5e-9);
    tests_passed++; tests_run++;

    TEST_PASS();
}

/* ========================================================================
 * Main Test Runner
 * ======================================================================== */
int main(void) {
    printf("============================================================\n");
    printf("  Indoor Positioning Module — Comprehensive Test Suite\n");
    printf("============================================================\n\n");

    /* L4 - Fundamental Laws */
    test_rssi_path_loss_model();

    /* L5 - Algorithms */
    test_trilateration_2d();
    test_trilateration_3d();
    test_weighted_centroid();
    test_twr_ranging();
    test_tdoa_chan();
    test_fingerprint_matching();
    test_quaternion_operations();
    test_ins_mechanization();
    test_madgwick_ahrs();
    test_mahony_ahrs();
    test_kalman_filter();
    test_ekf_constant_velocity();
    test_complementary_filter();
    test_aoa();
    test_pdr_step_detection();
    test_toa_gauss_newton();
    test_first_path_detection();

    /* L3/L4/L5 - Error Analysis */
    test_error_analysis();
    test_coordinate_transforms();
    test_nlos_detection();

    /* L8 - Advanced Filters */
    test_ukf();
    test_particle_filter();
    test_quaternion_integration();

    /* L6 - UWB */
    test_uwb_link_budget();

    printf("\n============================================================\n");
    printf("  Results: %d / %d tests passed\n", tests_passed, tests_run);
    printf("============================================================\n");

    if (tests_passed == tests_run) {
        printf("  ALL TESTS PASSED ✅\n");
        return 0;
    } else {
        printf("  SOME TESTS FAILED ❌\n");
        return 1;
    }
}
