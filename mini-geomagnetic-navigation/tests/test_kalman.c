/**
 * test_kalman.c -- Tests for Kalman filter
 */
#include "geomag_kalman.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

static int tests_run = 0, tests_passed = 0;
#define TEST(n) printf("  %s... ", n); tests_run++
#define PASS() do { printf("PASS\n"); tests_passed++; } while(0)
#define CHECK(c) do { if (!(c)) { printf("FAIL: %s line %d\n", __FILE__, __LINE__); return 1; } } while(0)

static int test_kalman_init_free(void) {
    TEST("kalman_init/free");
    KalmanFilter kf;
    CHECK(kalman_init(&kf, 4, 2) == 0);
    CHECK(kf.n_states == 4 && kf.n_meas == 2);
    CHECK(kf.x != NULL);
    CHECK(kf.P != NULL);
    CHECK(kf.F != NULL);
    CHECK(kf.K != NULL);
    kalman_free(&kf);
    PASS(); return 0;
}

static int test_kalman_predict(void) {
    TEST("kalman_predict");
    KalmanFilter kf;
    kalman_init(&kf, 3, 2);
    double x0[3] = {1.0, 2.0, 3.0};
    double P0[9] = {1,0,0,0,1,0,0,0,1};
    double F[9] = {1,0,0,0,1,0,0,0,1};
    kalman_set_initial(&kf, x0, P0);
    kalman_set_F(&kf, F);
    kalman_predict(&kf);
    CHECK(fabs(kf.x[0] - 1.0) < 0.01);
    CHECK(fabs(kf.x[1] - 2.0) < 0.01);
    CHECK(fabs(kf.x[2] - 3.0) < 0.01);
    kalman_free(&kf);
    PASS(); return 0;
}

static int test_ekf_init_free(void) {
    TEST("ekf_init/free");
    ExtendedKalmanFilter ekf;
    CHECK(ekf_init(&ekf, 6, 3) == 0);
    CHECK(ekf.x != NULL);
    CHECK(ekf.P != NULL);
    ekf_free(&ekf);
    PASS(); return 0;
}

static int test_position_accuracy(void) {
    TEST("ekf_position_accuracy");
    double P[225] = {0};
    P[0*15+0] = 100.0; P[1*15+1] = 25.0;
    double smaj, smin, orient;
    ekf_position_accuracy(P, 15, &smaj, &smin, &orient);
    CHECK(smaj > 0.0 && smin > 0.0);
    CHECK(smaj >= smin);
    PASS(); return 0;
}

int main(void) {
    printf("test_kalman\n");
    test_kalman_init_free();
    test_kalman_predict();
    test_ekf_init_free();
    test_position_accuracy();
    printf("  %d/%d passed\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
