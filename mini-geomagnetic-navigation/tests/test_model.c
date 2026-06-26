/**
 * test_model.c -- Minimal IGRF model tests
 */
#include "geomag_model.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

static int tests_run = 0, tests_passed = 0;
#define TEST(n) printf("  %s... ", n); tests_run++
#define PASS() do { printf("PASS\n"); tests_passed++; } while(0)
#define CHECK(c) do { if (!(c)) { printf("FAIL: %s line %d\n", __FILE__, __LINE__); return 1; } } while(0)

int main(void) {
    printf("test_model\n");

    TEST("igrf_init_model");
    IGRFModel model;
    CHECK(igrf_init_model(13, &model) == 0);
    CHECK(model.nmax == 13);
    CHECK(model.ncoeffs == 104);
    PASS();

    TEST("schmidt_index");
    CHECK(schmidt_index(1, 0) == 1);
    CHECK(schmidt_index(13, 13) == 104);
    PASS();

    TEST("alloc_legendre");
    LegendreState *ls = alloc_legendre_state(3);
    CHECK(ls != NULL);
    free_legendre_state(ls);
    PASS();

    TEST("field_compute");
    GeodeticCoord loc = {0.0, 0.0, 0.0};
    MagVector B;
    int rc = igrf_compute_field(&model, &loc, &B);
    CHECK(rc == 0);
    double F = mag_magnitude(&B);
    CHECK(F > 20000.0 && F < 70000.0);
    PASS();

    igrf_free_model(&model);
    printf("  %d/%d passed\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
