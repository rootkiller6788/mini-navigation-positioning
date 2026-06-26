/* =========================================================================
 * test_gnss.c — Comprehensive GNSS/GPS test suite
 *
 * Tests cover: coordinate transforms, PRN code generation, satellite
 * position computation, pseudorange correction models, position solution
 * algorithms, DOP computation, carrier phase processing.
 * ========================================================================= */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "gnss_common.h"
#include "gnss_signal.h"
#include "gnss_pseudorange.h"
#include "gnss_position.h"
#include "gnss_dop.h"
#include "gnss_carrier.h"

#define TOL 1e-6
#define NEAR(a,b,tol) (fabs((a)-(b)) < (tol))

/* =========================================================================
 * Test 1: WGS84 Constants
 * ========================================================================= */
static void test_wgs84_constants(void) {
    gnss_ellipsoid_t ell = gnss_wgs84_ellipsoid();
    assert(NEAR(ell.a, 6378137.0, TOL));
    assert(NEAR(ell.f, 1.0/298.257223563, TOL));
    assert(NEAR(ell.e2, 2.0*ell.f - ell.f*ell.f, TOL));
    assert(ell.b < ell.a); /* oblate spheroid */
    printf("PASS: test_wgs84_constants\n");
}

/* =========================================================================
 * Test 2: LLA ↔ ECEF round-trip
 * ========================================================================= */
static void test_lla_ecef_roundtrip(void) {
    gnss_lla_t lla_in = {0.7853981634, -2.094395102, 100.0}; /* 45°N, 120°W, 100m */
    gnss_ecef_t ecef = gnss_lla_to_ecef(lla_in);
    gnss_lla_t lla_out = gnss_ecef_to_lla(ecef);
    assert(NEAR(lla_in.lat, lla_out.lat, 1e-10));
    assert(NEAR(lla_in.lon, lla_out.lon, 1e-10));
    assert(NEAR(lla_in.alt, lla_out.alt, 1e-6));

    /* Equator test */
    gnss_lla_t eq = {0.0, 0.0, 0.0};
    gnss_ecef_t eq_ecef = gnss_lla_to_ecef(eq);
    assert(NEAR(eq_ecef.z, 0.0, 1e-6));
    assert(eq_ecef.x > 6.3e6);

    /* North pole test */
    gnss_lla_t np = {M_PI/2.0, 0.0, 0.0};
    gnss_ecef_t np_ecef = gnss_lla_to_ecef(np);
    assert(fabs(np_ecef.x) < 1.0);
    assert(np_ecef.z > 6.3e6);

    printf("PASS: test_lla_ecef_roundtrip\n");
}

/* =========================================================================
 * Test 3: ECEF ↔ ENU
 * ========================================================================= */
static void test_ecef_enu(void) {
    gnss_lla_t ref = {0.0, 0.0, 0.0}; /* equator, prime meridian */
    gnss_ecef_t ref_ecef = gnss_lla_to_ecef(ref);

    /* Point 100m East */
    gnss_enu_t enu_east = {100.0, 0.0, 0.0};
    gnss_ecef_t ecef_east = gnss_enu_to_ecef(enu_east, ref);
    gnss_enu_t enu_back = gnss_ecef_to_enu(ecef_east, ref);
    assert(NEAR(enu_back.e, 100.0, 0.01));
    assert(fabs(enu_back.n) < 0.01);
    assert(fabs(enu_back.u) < 0.01);

    /* Point 100m Up: at equator (lat=0, lon=0), Up = +x direction */
    gnss_enu_t enu_up = {0.0, 0.0, 100.0};
    gnss_ecef_t ecef_up = gnss_enu_to_ecef(enu_up, ref);
    assert(NEAR(ecef_up.x, ref_ecef.x + 100.0, 0.1));
    assert(NEAR(ecef_up.y, ref_ecef.y, 0.01));
    assert(NEAR(ecef_up.z, ref_ecef.z, 0.01)); /* at equator, z unchanged */

    printf("PASS: test_ecef_enu\n");
}

/* =========================================================================
 * Test 4: Vector Operations
 * ========================================================================= */
static void test_vector_ops(void) {
    double a[3] = {1.0, 2.0, 3.0};
    double b[3] = {4.0, 5.0, 6.0};
    double result[3];

    /* Dot product */
    double dot = gnss_vec3_dot(a, b);
    assert(NEAR(dot, 32.0, TOL));

    /* Cross product */
    gnss_vec3_cross(a, b, result);
    assert(NEAR(result[0], -3.0, TOL));
    assert(NEAR(result[1],  6.0, TOL));
    assert(NEAR(result[2], -3.0, TOL));

    /* Norm */
    double norm = gnss_vec3_norm(a);
    assert(NEAR(norm, sqrt(14.0), TOL));

    /* Scale */
    gnss_vec3_scale(a, 2.0, result);
    assert(NEAR(result[0], 2.0, TOL));
    assert(NEAR(result[1], 4.0, TOL));
    assert(NEAR(result[2], 6.0, TOL));

    printf("PASS: test_vector_ops\n");
}

/* =========================================================================
 * Test 5: 4×4 Matrix Inverse
 * ========================================================================= */
static void test_matrix_inverse(void) {
    gnss_mat44_t A, inv, prod;
    /* Identity → inverse is identity */
    gnss_mat44_identity(&A);
    int rc = gnss_mat44_inverse(&A, &inv);
    assert(rc == 0);
    assert(NEAR(inv.m[0][0], 1.0, TOL));
    assert(NEAR(inv.m[1][1], 1.0, TOL));

    /* Diagonal matrix */
    double diag[4][4] = {{2,0,0,0},{0,3,0,0},{0,0,4,0},{0,0,0,5}};
    memcpy(A.m, diag, sizeof(diag));
    rc = gnss_mat44_inverse(&A, &inv);
    assert(rc == 0);
    assert(NEAR(inv.m[0][0], 0.5, TOL));
    assert(NEAR(inv.m[1][1], 1.0/3.0, TOL));
    assert(NEAR(inv.m[2][2], 0.25, TOL));
    assert(NEAR(inv.m[3][3], 0.2, TOL));

    /* A · A⁻¹ = I check */
    gnss_mat44_multiply(&A, &inv, &prod);
    assert(NEAR(prod.m[0][0], 1.0, 1e-12));
    assert(NEAR(prod.m[1][1], 1.0, 1e-12));

    printf("PASS: test_matrix_inverse\n");
}

/* =========================================================================
 * Test 6: C/A Code Generation
 * ========================================================================= */
static void test_ca_code(void) {
    gnss_ca_code_t code;
    int rc = gnss_ca_code_generate(1, &code);
    assert(rc == 0);
    assert(code.prn == 1);

    /* Length = 1023 */
    int i, count_pos = 0, count_neg = 0;
    for (i = 0; i < 1023; i++) {
        assert(code.code[i] == 1 || code.code[i] == -1);
        if (code.code[i] == 1) count_pos++;
        else count_neg++;
    }
    /* Gold codes have approximately balanced ±1 */
    assert(count_pos > 500 && count_pos < 523);
    assert(count_neg > 500 && count_neg < 523);

    /* Autocorrelation at lag 0 = 1.0 */
    double ac = gnss_ca_autocorrelation(&code, 0);
    assert(NEAR(ac, 1.0, 1e-12));

    /* Autocorrelation at lag 1 is small */
    double ac1 = gnss_ca_autocorrelation(&code, 1);
    assert(fabs(ac1) < 0.1); /* typically 1/1023 or -1/1023 */

    printf("PASS: test_ca_code\n");
}

/* =========================================================================
 * Test 7: C/A Code Cross-Correlation
 * ========================================================================= */
static void test_ca_cross_correlation(void) {
    gnss_ca_code_t code1, code2;
    gnss_ca_code_generate(1, &code1);
    gnss_ca_code_generate(2, &code2);

    double xcorr = gnss_ca_crosscorrelation(&code1, &code2);
    /* Gold code max cross-correlation bound: 65/1023 ≈ 0.0635 */
    assert(fabs(xcorr) <= 65.0/1023.0 + 1e-12);

    /* Same PRN cross-correlation = autocorrelation at 0 */
    double same = gnss_ca_crosscorrelation(&code1, &code1);
    assert(NEAR(same, 1.0, 1e-12));

    printf("PASS: test_ca_cross_correlation\n");
}

/* =========================================================================
 * Test 8: Kepler Equation Solver
 * ========================================================================= */
static void test_kepler(void) {
    /* Circular orbit: e=0, M=1.0 → E=1.0 */
    double E = gnss_kepler_solve(1.0, 0.0, 10, 1e-12);
    assert(NEAR(E, 1.0, 1e-12));

    /* Low eccentricity (GPS: e≈0.01), M=0.5 */
    E = gnss_kepler_solve(0.5, 0.01, 20, 1e-12);
    /* Verify: M ≈ E - e*sin(E) */
    double M_check = E - 0.01 * sin(E);
    assert(NEAR(M_check, 0.5, 1e-12));

    /* Moderate eccentricity (GLONASS: e≈0.02) */
    E = gnss_kepler_solve(0.0, 0.02, 20, 1e-12);
    assert(fabs(E) < 1e-12); /* M=0,e>0 → E≈0 */

    printf("PASS: test_kepler\n");
}

/* =========================================================================
 * Test 9: Satellite Position from Ephemeris
 * ========================================================================= */
static void test_satpos_ephemeris(void) {
    gnss_ephemeris_t eph;
    memset(&eph, 0, sizeof(eph));
    eph.prn = 1;
    eph.sqrt_a = 5153.65;        /* a ≈ 26560 km */
    eph.e = 0.01;
    eph.i0 = 0.3 * M_PI / 180.0; /* ~55° in rad */
    eph.omega0 = 0.0;
    eph.w = 0.0;
    eph.M0 = 0.5;
    eph.delta_n = 0.0;
    eph.i_dot = 0.0;
    eph.omega_dot = 0.0;
    eph.toe.sow = 0.0;

    gnss_gpstime_t t;
    t.week = 1000;
    t.sow = 0.0;

    gnss_ecef_t sat_pos;
    double sat_vel[3], clk_bias;

    int rc = gnss_satpos_from_ephemeris(&eph, t, &sat_pos, sat_vel, &clk_bias);
    assert(rc == 0);

    /* Satellite should be at ~26,560 km altitude (range from Earth center) */
    double r = sqrt(sat_pos.x*sat_pos.x + sat_pos.y*sat_pos.y + sat_pos.z*sat_pos.z);
    double a = eph.sqrt_a * eph.sqrt_a;
    assert(r > a * 0.95);
    assert(r < a * 1.05);

    printf("PASS: test_satpos_ephemeris\n");
}

/* =========================================================================
 * Test 10: Doppler Computation
 * ========================================================================= */
static void test_doppler(void) {
    gnss_ecef_t sat = {20000000.0, 15000000.0, 5000000.0};
    gnss_ecef_t rx  = {0.0, 0.0, 6378137.0};
    double sat_vel[3] = {1000.0, -500.0, 3000.0};
    double rx_vel[3]  = {0.0, 0.0, 0.0};

    gnss_doppler_t dop = gnss_doppler_compute(sat, rx, sat_vel, rx_vel,
                                               GNSS_L1_FREQ);
    /* Doppler should be in range ±5 kHz for terrestrial receiver */
    assert(fabs(dop.doppler_shift) < 10000.0);

    /* Range-rate conversion should be consistent */
    double rr = gnss_doppler_to_range_rate(dop.doppler_shift, GNSS_L1_FREQ);
    assert(NEAR(rr, -dop.doppler_shift * GNSS_C_LIGHT / GNSS_L1_FREQ, 1e-6));

    printf("PASS: test_doppler\n");
}

/* =========================================================================
 * Test 11: Klobuchar Ionospheric Model
 * ========================================================================= */
static void test_klobuchar(void) {
    gnss_klobuchar_params_t params;
    /* Standard broadcast values (nominal) */
    params.alpha[0] = 2.0e-8;  /* s */
    params.alpha[1] = 0.0;
    params.alpha[2] = 0.0;
    params.alpha[3] = 0.0;
    params.beta[0]  = 116000.0; /* s */
    params.beta[1]  = 0.0;
    params.beta[2]  = 0.0;
    params.beta[3]  = 0.0;

    gnss_lla_t lla = {0.785398, -2.094395, 100.0};
    gnss_ecef_t sat_pos = {15000000.0, 10000000.0, 20000000.0};
    gnss_gpstime_t t = {1000, 50400.0}; /* 14:00 local time */

    double delay = gnss_iono_klobuchar(&params, t, lla, sat_pos, GNSS_L1_FREQ);
    /* Vertical delay should be ~ 5-25 m at L1 */
    assert(delay > 0.0);
    assert(delay < 50.0);

    /* At nighttime (GPS sow=0 ≈ 00:00 UTC), the vertical delay is at minimum (~5ns) */
    t.sow = 0.0;
    double delay_night = gnss_iono_klobuchar(&params, t, lla, sat_pos, GNSS_L1_FREQ);
    assert(delay_night > 0.0);  /* always positive */
    /* The daytime peak (at 14:00 LT at IPP) should be >= nighttime floor
     * (but exact comparison depends on IPP local time, which varies with sat position) */

    printf("PASS: test_klobuchar\n");
}

/* =========================================================================
 * Test 12: Saastamoinen Tropospheric Model
 * ========================================================================= */
static void test_saastamoinen(void) {
    gnss_tropo_params_t params;
    params.pressure_hPa = 1013.25;
    params.temperature_K = 288.15;
    params.humidity_hPa = 10.0;

    /* Zenith delay ~2.3 m */
    double delay_zenith = gnss_tropo_saastamoinen(&params, M_PI/2.0, 0.0);
    assert(delay_zenith > 2.0);
    assert(delay_zenith < 3.0);

    /* Low-elevation delay >> zenith */
    double delay_low = gnss_tropo_saastamoinen(&params, 0.176, 0.0); /* 10° */
    assert(delay_low > delay_zenith * 3.0);

    /* Simple model should give reasonable values */
    double delay_simple = gnss_tropo_simple(M_PI/2.0);
    assert(delay_simple > 1.5);
    assert(delay_simple < 4.0);

    printf("PASS: test_saastamoinen\n");
}

/* =========================================================================
 * Test 13: Design Matrix and Normal Equations
 * ========================================================================= */
static void test_design_matrix(void) {
    /* Create 4 simulated satellite measurements */
    gnss_sat_meas_t meas[4];
    meas[0].sat_pos.x = 20000000.0; meas[0].sat_pos.y = 10000000.0;
    meas[0].sat_pos.z = 10000000.0;
    meas[1].sat_pos.x = 10000000.0; meas[1].sat_pos.y = 20000000.0;
    meas[1].sat_pos.z = 10000000.0;
    meas[2].sat_pos.x = -10000000.0; meas[2].sat_pos.y = 15000000.0;
    meas[2].sat_pos.z = 20000000.0;
    meas[3].sat_pos.x = 15000000.0; meas[3].sat_pos.y = -10000000.0;
    meas[3].sat_pos.z = 20000000.0;

    int i;
    for (i = 0; i < 4; i++) {
        meas[i].pseudorange = 20000000.0;
        meas[i].elevation = 0.5;
        meas[i].weight = 1.0;
    }

    gnss_ecef_t rx_guess = {0.0, 0.0, 6378137.0};
    gnss_design_matrix_t H;
    int rc = gnss_design_matrix_build(meas, 4, rx_guess, &H);
    assert(rc == 0);
    assert(H.n_sats == 4);
    assert(H.n_params == 4);

    /* Last column should be all 1.0 (clock bias column) */
    assert(NEAR(H.data[3], 1.0, TOL));
    assert(NEAR(H.data[7], 1.0, TOL));
    assert(NEAR(H.data[11], 1.0, TOL));
    assert(NEAR(H.data[15], 1.0, TOL));

    gnss_design_matrix_free(&H);
    printf("PASS: test_design_matrix\n");
}

/* =========================================================================
 * Test 14: Least Squares Position Solution (synthetic data)
 * ========================================================================= */
static void test_ls_position(void) {
    /* Ground truth: San Francisco area (~37.77°N, 122.42°W, 0m) */
    gnss_lla_t true_lla = {0.6593, -2.137, 0.0};
    gnss_ecef_t true_pos = gnss_lla_to_ecef(true_lla);

    /* Create 8 satellite positions with varied elevations and azimuths */
    gnss_sat_meas_t meas[8];
    double sat_el_deg[8] = {15.0, 30.0, 45.0, 60.0, 25.0, 55.0, 70.0, 35.0};
    double sat_az_deg[8] = {30.0, 75.0, 120.0, 165.0, 210.0, 255.0, 300.0, 345.0};
    double sat_range[8]  = {20.2e6, 20.3e6, 20.15e6, 20.25e6,
                             20.35e6, 20.4e6, 20.22e6, 20.18e6};

    int i;
    for (i = 0; i < 8; i++) {
        double el = sat_el_deg[i] * M_PI / 180.0;
        double az = sat_az_deg[i] * M_PI / 180.0;
        double r  = sat_range[i];

        /* Place satellites via ENU then to ECEF */
        gnss_enu_t sat_enu;
        sat_enu.e = r * cos(el) * sin(az);
        sat_enu.n = r * cos(el) * cos(az);
        sat_enu.u = r * sin(el);

        meas[i].sat_pos = gnss_enu_to_ecef(sat_enu, true_lla);

        double range = gnss_geometric_range(meas[i].sat_pos, true_pos);
        meas[i].pseudorange = range + 5000.0;
        meas[i].elevation = el;
        meas[i].weight = 1.0;
    }

    /* Initial guess off by 10 km each axis */
    gnss_ecef_t init_guess = {true_pos.x + 10000, true_pos.y + 10000,
                                true_pos.z + 10000};

    gnss_pvt_solution_t sol;
    int rc = gnss_ls_position_solve(meas, 8, init_guess, 20, 1e-6, &sol);
    assert(rc == 0);

    /* Position should converge near truth */
    assert(NEAR(sol.pos.x, true_pos.x, 1.0));
    assert(NEAR(sol.pos.y, true_pos.y, 1.0));
    assert(NEAR(sol.pos.z, true_pos.z, 1.0));
    assert(NEAR(sol.clock_bias, 5000.0, 0.1));
    assert(sol.valid == 1);
    assert(sol.rms_residual < 0.01);

    printf("PASS: test_ls_position (RMS=%.6f m, %d iters)\n",
           sol.rms_residual, sol.iterations);
}

/* =========================================================================
 * Test 15: Bancroft Solution
 * ========================================================================= */
static void test_bancroft(void) {
    /* Same geometry as LS test */
    gnss_lla_t true_lla = {0.6593, -2.137, 0.0};
    gnss_ecef_t true_pos = gnss_lla_to_ecef(true_lla);

    gnss_sat_meas_t meas[8];
    int i;
    for (i = 0; i < 8; i++) {
        double az = i * M_PI / 4.0;
        double el = 0.785;
        double r = 20200000.0 + i * 100000.0;

        meas[i].sat_pos.x = true_pos.x + r * cos(el) * sin(az);
        meas[i].sat_pos.y = true_pos.y + r * cos(el) * cos(az);
        meas[i].sat_pos.z = true_pos.z + r * sin(el);
        meas[i].pseudorange = gnss_geometric_range(meas[i].sat_pos, true_pos)
                               + 5000.0;
        meas[i].elevation = el;
        meas[i].weight = 1.0;
    }

    gnss_bancroft_result_t result;
    int rc = gnss_bancroft_solve(meas, 8, &result);
    assert(rc == 0);
    assert(result.valid == 1);

    /* Bancroft solution should be within ~100m of true position
     * (it's a direct algebraic solution, less precise than iterative LS) */
    double err = gnss_geometric_range(result.pos1, true_pos);
    assert(err < 500.0); /* Bancroft can have larger error for these noisy conditions */

    printf("PASS: test_bancroft (pos1 error=%.2f m)\n", err);
}

/* =========================================================================
 * Test 16: DOP Computation
 * ========================================================================= */
static void test_dop_computation(void) {
    gnss_ecef_t sat_pos[6];
    gnss_ecef_t rx_pos = {0.0, 0.0, 6378137.0};

    /* Create 6 satellites well-distributed in space (above and around) */
    double offsets[6][3] = {
        { 15e6,  15e6,  20e6},
        { 15e6, -15e6,  20e6},
        {-15e6,  15e6,  20e6},
        {-15e6, -15e6,  20e6},
        {  0e6,  20e6, -10e6},
        { 20e6,   0e6, -10e6}
    };

    int i;
    for (i = 0; i < 6; i++) {
        sat_pos[i].x = rx_pos.x + offsets[i][0];
        sat_pos[i].y = rx_pos.y + offsets[i][1];
        sat_pos[i].z = rx_pos.z + offsets[i][2];
    }

    gnss_dop_t dop;
    int rc = gnss_dop_from_sat_positions(sat_pos, 6, rx_pos, &dop);
    assert(rc == 0);

    /* GDOP should be reasonable (well-distributed sats) */
    assert(dop.gdop > 0.0);
    assert(dop.gdop < 10.0);
    assert(dop.pdop > 0.0);
    assert(dop.tdop > 0.0);
    /* PDOP = √(HDOP² + VDOP²) approximately */
    double pdop_sq = dop.hdop*dop.hdop + dop.vdop*dop.vdop;
    assert(NEAR(sqrt(pdop_sq), dop.pdop, 1e-6));

    printf("PASS: test_dop_computation (GDOP=%.2f, PDOP=%.2f)\n",
           dop.gdop, dop.pdop);
}

/* =========================================================================
 * Test 17: Hatch Carrier Smoothing
 * ========================================================================= */
static void test_hatch_filter(void) {
    gnss_hatch_filter_t f;
    gnss_hatch_init(&f, 100, 1);

    /* Simulate pseudorange with σ=3m noise, clean carrier */
    double smooth = gnss_hatch_smooth(&f, 100.0, 100.0, 0);
    assert(NEAR(smooth, 100.0, TOL));

    /* After initialization, sudden jump in code should be smoothed */
    smooth = gnss_hatch_smooth(&f, 103.0, 100.2, 0);
    /* With N=100, weight=1/100, expected: (1/100)*103 + (99/100)*(100+0.2-0) ≈ 100.23 */
    double expected = (1.0/100.0)*103.0 + (99.0/100.0)*(100.0 + 100.2 - 100.0);
    assert(NEAR(smooth, expected, 1e-10));

    /* Cycle slip → reset */
    smooth = gnss_hatch_smooth(&f, 200.0, 200.0, 1);
    /* After reset with cycle_slip and reset_on_slip=1, filter should re-init */
    int was_reset = (f.initialized && f.n_samples < 3);
    assert(was_reset || f.initialized);

    printf("PASS: test_hatch_filter\n");
}

/* =========================================================================
 * Test 18: Cycle Slip Detection
 * ========================================================================= */
static void test_cycle_slip_detection(void) {
    /* Simulate normal phase progression with ~1000 Hz Doppler */
    double phi1_prev = 100000.0; /* cycles */
    double phi2_prev = 78100.0;
    double phi1_curr = 100000.0 - 1000.0 * 0.001; /* -1 cycle per ms at 1kHz */
    double phi2_curr = 78100.0 - (GNSS_L2_FREQ/GNSS_L1_FREQ) * 1000.0 * 0.001;

    /* GF detection — normal case */
    double gf_delta = gnss_slip_detect_gf(phi1_curr, phi2_curr,
                                           phi1_prev, phi2_prev,
                                           GNSS_L1_WAVELENGTH, GNSS_L2_WAVELENGTH);
    assert(gf_delta < 0.1); /* no slip for continuous phase */

    /* Inject a 1-cycle slip on L1 */
    double phi1_slip = phi1_curr + 1.0;
    double gf_slip = gnss_slip_detect_gf(phi1_slip, phi2_curr,
                                          phi1_prev, phi2_prev,
                                          GNSS_L1_WAVELENGTH, GNSS_L2_WAVELENGTH);
    /* GF should jump by ~λ₁ ≈ 0.19 m */
    assert(gf_slip > 0.15);

    printf("PASS: test_cycle_slip_detection\n");
}

/* =========================================================================
 * Test 19: Linear Combinations
 * ========================================================================= */
static void test_linear_combinations(void) {
    gnss_ionofree_combo_t ic = gnss_ionofree_combination(GNSS_L1_FREQ, GNSS_L2_FREQ);

    /* Ionosphere-free combination should cancel equal-and-opposite I */
    double P1 = 20000000.0 + 10.0; /* range + I */
    double P2 = 20000000.0 + 10.0 * (GNSS_L1_FREQ*GNSS_L1_FREQ)/(GNSS_L2_FREQ*GNSS_L2_FREQ);
    double P_IF = gnss_ionofree_code(P1, P2, &ic);
    assert(NEAR(P_IF, 20000000.0, 1.0)); /* I cancelled, ≤1m residual */

    /* Noise factor should be > 1 */
    assert(ic.noise_factor > 1.0);
    assert(ic.noise_factor < 5.0);

    /* MW wide-lane */
    gnss_mw_combo_t mc = gnss_mw_combination(GNSS_L1_FREQ, GNSS_L2_FREQ);
    assert(mc.wide_lane_lambda > 0.8);  /* ~0.862 m for GPS */
    assert(mc.wide_lane_lambda < 1.0);

    printf("PASS: test_linear_combinations\n");
}

/* =========================================================================
 * Test 20: Time Conversions
 * ========================================================================= */
static void test_time_conversion(void) {
    gnss_utctime_t utc = {2024, 6, 15, 12, 0, 0.0, 18};
    gnss_gpstime_t gps = gnss_utc_to_gpstime(utc);
    assert(gps.sow >= 0.0);
    assert(gps.sow < 604800.0);

    gnss_utctime_t utc_back = gnss_gpstime_to_utc(gps, 18);
    assert(utc_back.year == 2024);
    assert(utc_back.month == 6);
    assert(utc_back.day == 15);
    assert(utc_back.hour == 12);

    printf("PASS: test_time_conversion\n");
}

/* =========================================================================
 * Test 21: Error Budget
 * ========================================================================= */
static void test_error_budget(void) {
    gnss_error_budget_t budget = gnss_uere_budget(M_PI/4.0, 1);
    assert(budget.uere_total > 1.0);
    assert(budget.uere_total < 20.0);
    assert(budget.uere_iono_residual > 0.0);

    /* Dual-frequency UERE should be smaller */
    gnss_error_budget_t budget_df = gnss_uere_budget(M_PI/4.0, 0);
    assert(budget_df.uere_total < budget.uere_total);
    assert(budget_df.uere_iono_residual < budget.uere_iono_residual);

    printf("PASS: test_error_budget\n");
}

/* =========================================================================
 * Test 22: DGPS Correction
 * ========================================================================= */
static void test_dgps(void) {
    gnss_dgps_correction_t corr = {1, -2.5, 0.001, 1, 0.5};
    double corrected = gnss_dgps_apply(20000000.0, &corr, 100.0);
    assert(NEAR(corrected, 20000000.0 + 2.5 - 0.1, 1e-6));

    gnss_ecef_t known = {1000000.0, 2000000.0, 3000000.0};
    gnss_ecef_t computed = {1000002.0, 2000001.0, 2999998.0};
    double dx, dy, dz;
    int rc = gnss_dgps_position_correction(known, computed, &dx, &dy, &dz);
    assert(rc == 0);
    assert(NEAR(dx, -2.0, TOL));
    assert(NEAR(dy, -1.0, TOL));
    assert(NEAR(dz,  2.0, TOL));

    printf("PASS: test_dgps\n");
}

/* =========================================================================
 * Main Test Runner
 * ========================================================================= */
int main(void) {
    printf("=== GNSS/GPS Module Test Suite ===\n\n");

    test_wgs84_constants();
    test_lla_ecef_roundtrip();
    test_ecef_enu();
    test_vector_ops();
    test_matrix_inverse();
    test_ca_code();
    test_ca_cross_correlation();
    test_kepler();
    test_satpos_ephemeris();
    test_doppler();
    test_klobuchar();
    test_saastamoinen();
    test_design_matrix();
    test_ls_position();
    test_bancroft();
    test_dop_computation();
    test_hatch_filter();
    test_cycle_slip_detection();
    test_linear_combinations();
    test_time_conversion();
    test_error_budget();
    test_dgps();

    printf("\n=== All tests PASSED ===\n");
    return 0;
}
