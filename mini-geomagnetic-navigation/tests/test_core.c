/**
 * test_core.c -- Tests for geomag_core
 */
#include "geomag_core.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>

static int tests_run = 0, tests_passed = 0;
#define TEST(name) printf("  %s... ", name); tests_run++
#define PASS() do { printf("PASS\n"); tests_passed++; } while(0)
#define CHECK(cond) do { if (!(cond)) { printf("FAIL at line %d\n", __LINE__); return 1; } } while(0)

static int test_geodetic_to_ecef(void) {
    TEST("geodetic_to_ecef");
    GeodeticCoord geo = { 0.0, 0.0, 0.0 };
    ECEFCoord ecef;
    geodetic_to_ecef(&geo, &ecef);
    CHECK(fabs(ecef.x - WGS84_A) < 1.0);
    CHECK(fabs(ecef.y) < 1.0);
    CHECK(fabs(ecef.z) < 1.0);
    geo.lat = 90.0;
    geodetic_to_ecef(&geo, &ecef);
    CHECK(fabs(ecef.z - WGS84_B) < 1.0);
    PASS(); return 0;
}

static int test_ecef_to_geodetic(void) {
    TEST("ecef_to_geodetic");
    ECEFCoord ecef = { WGS84_A, 0.0, 0.0 };
    GeodeticCoord geo;
    ecef_to_geodetic(&ecef, &geo);
    CHECK(fabs(geo.lat) < 1e-6);
    CHECK(fabs(geo.lon) < 1e-6);
    PASS(); return 0;
}

static int test_mag_magnitude(void) {
    TEST("mag_magnitude");
    MagVector B = { 3.0, 4.0, 0.0 };
    CHECK(fabs(mag_magnitude(&B) - 5.0) < 1e-10);
    B.bz = 12.0;
    CHECK(fabs(mag_magnitude(&B) - 13.0) < 1e-10);
    PASS(); return 0;
}

static int test_magnetic_elements(void) {
    TEST("magnetic_elements");
    MagVector B = { 20000.0, 0.0, 30000.0 };
    MagneticElements elem;
    compute_magnetic_elements(&B, &elem);
    CHECK(fabs(elem.north_component - 20000.0) < 0.1);
    CHECK(fabs(elem.vertical - 30000.0) < 0.1);
    double H_expected = 20000.0;
    double F_expected = sqrt(20000.0*20000.0 + 30000.0*30000.0);
    CHECK(fabs(elem.horizontal - H_expected) < 0.1);
    CHECK(fabs(elem.total_intensity - F_expected) < 0.1);
    PASS(); return 0;
}

static int test_calibration(void) {
    TEST("magnetometer_calibrate");
    double raw[3] = { 101.0, 98.0, 103.0 };
    double bias[3] = { 1.0, -2.0, 3.0 };
    double scale_inv[9] = { 1,0,0, 0,1,0, 0,0,1 };
    double cal[3];
    magnetometer_calibrate(raw, bias, scale_inv, cal);
    CHECK(fabs(cal[0] - 100.0) < 0.01);
    CHECK(fabs(cal[1] - 100.0) < 0.01);
    CHECK(fabs(cal[2] - 100.0) < 0.01);
    PASS(); return 0;
}

int main(void) {
    printf("test_core\n");
    test_geodetic_to_ecef();
    test_ecef_to_geodetic();
    test_mag_magnitude();
    test_magnetic_elements();
    test_calibration();
    printf("  %d/%d passed\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
