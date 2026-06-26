/**
 * test_sensor.c -- Tests for magnetometer sensor models
 */
#include "geomag_sensor.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>

static int tests_run = 0, tests_passed = 0;
#define TEST(n) printf("  %s... ", n); tests_run++
#define PASS() do { printf("PASS\n"); tests_passed++; } while(0)
#define CHECK(c) do { if (!(c)) { printf("FAIL: %s line %d\n", __FILE__, __LINE__); return 1; } } while(0)

static int test_scalar_read(void) {
    TEST("scalar_magnetometer_read");
    MagVector B = {30000.0, 10000.0, 40000.0};
    double reading;
    scalar_magnetometer_read(&B, 5.0, 0.1, &reading);
    double expected = sqrt(30000.0*30000.0+10000.0*10000.0+40000.0*40000.0)+5.0;
    CHECK(fabs(reading - expected) < 0.01);
    PASS(); return 0;
}

static int test_triaxial_read(void) {
    TEST("triaxial_magnetometer_read");
    MagVector B = {20000.0, 5000.0, 35000.0};
    double bias[3] = {0,0,0}, scale[9] = {1,0,0,0,1,0,0,0,1}, si[9] = {1,0,0,0,1,0,0,0,1};
    double reading[3];
    triaxial_magnetometer_read(&B, bias, scale, si, 0.1, reading);
    CHECK(fabs(reading[0] - 20000.0) < 0.01);
    CHECK(fabs(reading[1] - 5000.0) < 0.01);
    CHECK(fabs(reading[2] - 35000.0) < 0.01);
    PASS(); return 0;
}

static int test_mad_threshold(void) {
    TEST("mad_threshold_detect");
    CHECK(mad_threshold_detect(50500.0, 50000.0, 100.0) == 1);
    CHECK(mad_threshold_detect(50050.0, 50000.0, 100.0) == 0);
    PASS(); return 0;
}

static int test_mad_snr(void) {
    TEST("mad_snr_db");
    CHECK(fabs(mad_snr_db(10.0, 1.0) - 20.0) < 0.01);
    PASS(); return 0;
}

int main(void) {
    printf("test_sensor\n");
    test_scalar_read();
    test_triaxial_read();
    test_mad_threshold();
    test_mad_snr();
    printf("  %d/%d passed\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
