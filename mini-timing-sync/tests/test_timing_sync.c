/**
 * @file test_timing_sync.c
 * @brief Comprehensive tests for timing sync library
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "timing_sync.h"
#include "clock_model.h"
#include "phase_detector.h"
#include "ptp_engine.h"
#include "ntp_client.h"
#include "allan_variance.h"
#include "time_transfer.h"

static int tests_run = 0;
static int tests_passed = 0;

#define TEST(name) do { tests_run++; printf("  TEST: %s ... ", name); } while(0)
#define PASS() do { tests_passed++; printf("PASS\n"); } while(0)

static void test_timestamp_diff_ns(void) {
    TEST("timestamp_diff_ns");
    Timestamp ts1 = {100, 500000000};
    Timestamp ts2 = {100, 0};
    double diff = timing_timestamp_diff_ns(&ts1, &ts2);
    assert(fabs(diff - 5.0e8) < 1.0);
    PASS();
}

static void test_timestamp_cmp(void) {
    TEST("timestamp_cmp");
    Timestamp a = {100, 0};
    Timestamp b = {100, 100};
    assert(timing_timestamp_cmp(&a, &b) < 0);
    assert(timing_timestamp_cmp(&a, &a) == 0);
    assert(timing_timestamp_cmp(&b, &a) > 0);
    PASS();
}

static void test_timestamp_valid(void) {
    TEST("timestamp_valid");
    Timestamp valid_ts = {100, 500};
    Timestamp invalid_ns = {100, 1000000000};
    assert(timing_timestamp_valid(&valid_ts) == 1);
    assert(timing_timestamp_valid(&invalid_ns) == 0);
    assert(timing_timestamp_valid(NULL) == 0);
    PASS();
}

static void test_ptp_offset_delay(void) {
    TEST("ptp_offset_delay");
    PtpTimestamps ts;
    ts.t1.seconds = 1000; ts.t1.nanoseconds = 0;
    ts.t2 = ts.t1; timing_timestamp_add_ns(&ts.t2, 1100);
    ts.t3 = ts.t2; timing_timestamp_add_ns(&ts.t3, 1000000);
    ts.t4 = ts.t3; timing_timestamp_add_ns(&ts.t4, 900);
    double offset, delay;
    int ret = timing_compute_offset_delay(&ts, &offset, &delay);
    assert(ret == 0);
    assert(fabs(offset - 100.0) < 1.0);
    assert(fabs(delay - 1000.0) < 1.0);
    PASS();
}

static void test_ntp_offset_delay(void) {
    TEST("ntp_offset_delay");
    Timestamp T1 = {1000, 0};
    Timestamp T2, T3, T4;
    T2 = T1; timing_timestamp_add_ns(&T2, 2050);
    T3 = T2; timing_timestamp_add_ns(&T3, 1000000);
    T4 = T3; timing_timestamp_add_ns(&T4, 1950);
    double offset, delay;
    int ret = timing_ntp_offset_delay(&T1, &T2, &T3, &T4, &offset, &delay);
    assert(ret == 0);
    assert(fabs(offset - 50.0) < 1.0);
    PASS();
}

static void test_pi_servo(void) {
    TEST("pi_servo");
    PiServoConfig cfg;
    PiServoState state;
    pi_servo_init(&cfg, &state);
    double corr = pi_servo_update(&cfg, &state, 100.0);
    assert(fabs(corr) < 0.01);
    corr = pi_servo_update(&cfg, &state, 100.0);
    assert(corr > 50.0 && corr < 150.0);
    PASS();
}

static void test_kalman_clock(void) {
    TEST("kalman_clock");
    ClockState state;
    clock_state_init(&state, 100.0, 10.0, 0.0, 50.0, 5.0);
    assert(fabs(state.offset_ns - 100.0) < 0.01);
    double Q[9] = {1e-4, 0, 0, 0, 1e-6, 0, 0, 0, 1e-10};
    clock_kalman_predict(&state, 1.0, Q);
    assert(fabs(state.offset_ns - 110.0) < 1.0);
    clock_kalman_update(&state, 112.0, 25.0);
    assert(fabs(state.offset_ns - 110.0) < 5.0);
    PASS();
}

static void test_clock_regression(void) {
    TEST("clock_linear_regression");
    double t_ref[10], t_device[10];
    for (int i = 0; i < 10; i++) {
        t_ref[i] = (double)i;
        t_device[i] = 2.0 * t_ref[i] + 100.0;
    }
    double slope, intercept, r2;
    int ret = clock_linear_regression(t_ref, t_device, 10, &slope, &intercept, &r2);
    assert(ret == 0);
    assert(fabs(slope - 2.0) < 0.001);
    assert(r2 > 0.999);
    PASS();
}

static void test_bmca(void) {
    TEST("bmca_compare");
    BmcaDataset A = {1, {6, 0x21, 100}, 128, {{0,1,2,3,4,5,6,7}}, 0};
    BmcaDataset B = {2, {6, 0x21, 100}, 128, {{0,1,2,3,4,5,6,7}}, 0};
    BmcaResult r = bmca_compare_datasets(&A, &B);
    assert(r == BMCA_THIS_IS_BETTER);
    PASS();
}

static void test_allan_variance(void) {
    TEST("allan_variance");
    int N = 1000;
    double *yi = (double*)malloc((size_t)N * sizeof(double));
    for (int i = 0; i < N; i++) {
        yi[i] = 1e-10 * ((double)rand() / RAND_MAX - 0.5);
    }
    AllanResult result;
    int ret = allan_variance_compute(yi, N, 1.0, 10.0, &result);
    assert(ret == 0);
    assert(result.adev > 0.0);
    free(yi);
    PASS();
}

static void test_allan_noise_type(void) {
    TEST("allan_noise_type");
    double adev1 = 1e-10;
    double adev10 = adev1 / sqrt(10.0);
    NoiseType nt = allan_identify_noise_type(adev1, 1.0, adev10, 10.0);
    assert(nt == NOISE_WHITE_FM);
    PASS();
}

static void test_three_cornered_hat(void) {
    TEST("three_cornered_hat");
    double sA2, sB2, sC2;
    int ret = three_cornered_hat(5.0, 7.0, 8.0, &sA2, &sB2, &sC2);
    assert(ret == 0);
    assert(sA2 == 3.0);
    assert(sB2 == 2.0);
    assert(sC2 == 5.0);
    PASS();
}

static void test_sagnac(void) {
    TEST("sagnac_correction");
    double sagnac = sagnac_correction_ns(35.7, 139.8, 37.8, -122.4);
    assert(fabs(sagnac) > 10.0);
    assert(fabs(sagnac) < 500.0);
    PASS();
}

static void test_holdover(void) {
    TEST("holdover");
    HoldoverConfig cfg = {500.0, 100.0, 86400.0, 1000.0};
    assert(holdover_should_enter(200.0, &cfg) == 0);
    assert(holdover_should_enter(600.0, &cfg) == 1);
    assert(holdover_can_exit(50.0, &cfg, 100.0) == 1);
    double unc = holdover_estimate_uncertainty(3600.0, 1.0, 0.0, 10.0);
    assert(unc > 10.0);
    PASS();
}

static void test_clock_model(void) {
    TEST("clock_model");
    ClockParameters params = {0.0, 100.0, 1.0, 0.0, 50.0};
    double err = clock_model_time_error(&params, 86400.0);
    assert(err > 8000.0);
    double offset = clock_offset_evolution(100.0, 10.0, 100.0);
    assert(fabs(offset - 1100.0) < 1.0);
    PASS();
}

static void test_5g_fronthaul(void) {
    TEST("5g_fronthaul");
    assert(ptp_5g_fronthaul_check(50.0, 'A') == 1);
    assert(ptp_5g_fronthaul_check(600.0, 'A') == 0);
    assert(ptp_5g_fronthaul_check(10.0, 'C') == 1);
    PASS();
}

static void test_power_grid(void) {
    TEST("power_grid_timing");
    assert(ptp_power_grid_check(500.0) == 1);
    assert(ptp_power_grid_check(2000.0) == 0);
    PASS();
}

static void test_gpsdo(void) {
    TEST("gpsdo_stability");
    assert(gpsdo_stability_check(1e-12, 1e-13, 1e-14, 'U') == 1);
    PASS();
}

static void test_mifid(void) {
    TEST("mifid_compliance");
    assert(ntp_mifid_compliance(50.0, 100000.0) == 1);
    assert(ntp_mifid_compliance(200000.0, 100000.0) == 0);
    PASS();
}

static void test_ptp_slave(void) {
    TEST("ptp_slave");
    PtpSlaveState state;
    ptp_slave_init(&state, 0.0, 1.0);
    assert(state.sync_status == SYNC_FREE_RUNNING);
    PtpTimestamps ts;
    ts.t1.seconds = 1000; ts.t1.nanoseconds = 0;
    ts.t2 = ts.t1; timing_timestamp_add_ns(&ts.t2, 1100);
    ts.t3 = ts.t2; timing_timestamp_add_ns(&ts.t3, 1000000);
    ts.t4 = ts.t3; timing_timestamp_add_ns(&ts.t4, 900);
    double corr = ptp_slave_update(&state, &ts);
    (void)corr;
    assert(state.sync_status == SYNC_ACQUIRING);
    PASS();
}

static void test_ntp_client(void) {
    TEST("ntp_client");
    NtpClient client;
    ntp_client_init(&client);
    int idx = ntp_add_peer(&client, 0x0A000001, 2);
    assert(idx == 0);
    NtpSample sample;
    memset(&sample, 0, sizeof(sample));
    sample.orig.seconds = 1000;
    sample.recv.seconds = 1000; sample.recv.nanoseconds = 2050;
    sample.xmit.seconds = 1001;
    sample.dest.seconds = 1001; sample.dest.nanoseconds = 1950;
    double corr = ntp_client_update(&client, 0, &sample);
    (void)corr;
    assert(client.sync_status == SYNC_LOCKED);
    PASS();
}

static void test_dpll(void) {
    TEST("dpll");
    DpllState dpll;
    dpll_init(&dpll, 10e6, 100e6, 32);
    LoopFilterCoeffs coeffs = {0.5, 0.01, 0.0, 0.01};
    double freq_out = 0.0;
    for (int i = 0; i < 100; i++) {
        dpll_update(0.0, &coeffs, &dpll, &freq_out);
    }
    assert(freq_out > 0.0);
    /* Lock detection needs N_consecutive detections below threshold */
    int locked = 0;
    for (int i = 0; i < 20; i++) {
        locked = dpll_lock_detect(0.001, 0.0001, 10, &dpll);
    }
    assert(locked == 1);
    PASS();
}

int main(void) {
    printf("=== Timing Sync Library Test Suite ===\n\n");
    printf("L1: Timestamp Utilities\n");
    test_timestamp_diff_ns();
    test_timestamp_cmp();
    test_timestamp_valid();
    printf("\nL2: Core Concepts\n");
    test_pi_servo();
    test_dpll();
    printf("\nL3-L4: Kalman / Fundamental Laws\n");
    test_kalman_clock();
    test_ptp_offset_delay();
    test_ntp_offset_delay();
    printf("\nL5: Algorithms\n");
    test_clock_regression();
    test_bmca();
    test_allan_variance();
    test_allan_noise_type();
    test_three_cornered_hat();
    test_sagnac();
    printf("\nL6: Canonical Problems\n");
    test_holdover();
    test_clock_model();
    test_ptp_slave();
    test_ntp_client();
    printf("\nL7: Applications\n");
    test_5g_fronthaul();
    test_power_grid();
    test_gpsdo();
    test_mifid();
    printf("\n=== Results: %d/%d tests passed ===\n", tests_passed, tests_run);
    if (tests_passed == tests_run) {
        printf("ALL TESTS PASSED\n");
        return 0;
    }
    printf("SOME TESTS FAILED\n");
    return 1;
}
