/**
 * test_navigation.c -- Tests for magnetic navigation
 */
#include "geomag_navigation.h"
#include "geomag_math.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>

static int tests_run = 0, tests_passed = 0;
#define TEST(n) printf("  %s... ", n); tests_run++
#define PASS() do { printf("PASS\n"); tests_passed++; } while(0)
#define CHECK(c) do { if (!(c)) { printf("FAIL: %s line %d\n", __FILE__, __LINE__); return 1; } } while(0)

static int test_heading_basic(void) {
    TEST("mag_heading basic");
    MagVector B = { 1.0, 0.0, 0.0 };
    CHECK(fabs(mag_heading_from_triaxial(&B, 0.0)) < 0.1);
    B.bx = 0.0; B.by = -1.0;
    CHECK(fabs(mag_heading_from_triaxial(&B, 0.0) - 90.0) < 0.1);
    PASS(); return 0;
}

static int test_tilt_compensated_heading(void) {
    TEST("tilt_compensated_heading");
    MagVector B = { 20000.0, 0.0, 30000.0 };
    CHECK(fabs(tilt_compensated_heading(&B, 0.0, 0.0, 0.0)) < 0.1);
    PASS(); return 0;
}

static int test_compass_dr(void) {
    TEST("mag_compass_dr_update");
    NavSolution state;
    state.position.lat = 0.0; state.position.lon = 0.0; state.position.alt = 0.0;
    state.timestamp = 0.0;
    mag_compass_dr_update(&state, 90.0, 1.0, 3600.0);
    CHECK(state.position.lon > 0.0);
    CHECK(fabs(state.position.lat) < 0.1);
    PASS(); return 0;
}

static int test_great_circle(void) {
    TEST("great_circle_distance");
    GeodeticCoord p1 = {0.0, 0.0, 0.0}, p2 = {0.0, 1.0, 0.0};
    CHECK(fabs(great_circle_distance(&p1, &p2) - 111319.0) < 1000.0);
    PASS(); return 0;
}

int main(void) {
    printf("test_navigation\n");
    test_heading_basic();
    test_tilt_compensated_heading();
    test_compass_dr();
    test_great_circle();
    printf("  %d/%d passed\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
