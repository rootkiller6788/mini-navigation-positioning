/* =========================================================================
 * example_position_solve.c — End-to-end GPS single-point positioning
 *
 * Demonstrates: L4 (pseudorange equations), L5 (iterative least squares),
 * L6 (SPP canonical problem), L7 (urban canyon GPS receiver simulation).
 *
 * Simulates a GPS receiver at San Francisco Airport (SFO) with 8 visible
 * satellites. Applies atmospheric corrections and solves for position.
 * ========================================================================= */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "gnss_common.h"
#include "gnss_signal.h"
#include "gnss_pseudorange.h"
#include "gnss_position.h"
#include "gnss_dop.h"

int main(void) {
    printf("=== GPS Single-Point Positioning Example ===\n\n");

    /* Ground truth: SFO (37.6189°N, 122.3750°W, 4m) */
    gnss_lla_t true_lla = {0.65669, -2.1359, 4.0};
    gnss_ecef_t true_pos = gnss_lla_to_ecef(true_lla);

    printf("True position (WGS84):\n");
    printf("  Lat = %.6f deg, Lon = %.6f deg, Alt = %.1f m\n",
           true_lla.lat * 180.0/M_PI, true_lla.lon * 180.0/M_PI, true_lla.alt);
    printf("  ECEF: (%.1f, %.1f, %.1f) m\n\n",
           true_pos.x, true_pos.y, true_pos.z);

    /* Simulate 8 visible satellites */
    printf("Simulating 8 visible GPS satellites...\n\n");

    gnss_sat_meas_t meas[8];
    int n_sats = 8;
    int i;

    double sat_az_deg[8]  = {30, 75, 120, 165, 210, 255, 300, 345};
    double sat_el_deg[8]  = {15, 45, 60, 30, 25, 55, 70, 40};
    double sat_orbit_r[8] = {26.0e6, 26.3e6, 25.8e6, 26.1e6,
                              26.2e6, 25.9e6, 26.4e6, 26.0e6};
    double sat_clock_err_ns[8] = {1.5, -0.8, 2.3, -1.2, 0.9, -0.4, 1.8, -0.6}; /* nanoseconds */

    for (i = 0; i < n_sats; i++) {
        double az = sat_az_deg[i] * M_PI / 180.0;
        double el = sat_el_deg[i] * M_PI / 180.0;
        double r  = sat_orbit_r[i];

        gnss_enu_t sat_enu;
        sat_enu.e = r * cos(el) * sin(az);
        sat_enu.n = r * cos(el) * cos(az);
        sat_enu.u = r * sin(el);

        meas[i].sat_pos = gnss_enu_to_ecef(sat_enu, true_lla);

        double true_range = gnss_geometric_range(meas[i].sat_pos, true_pos);
        double sat_clk_m = sat_clock_err_ns[i] * 1e-9 * GNSS_C_LIGHT;

        double el_rad = el;
        double iono = 5.0 / sin(el_rad + 0.1);
        double tropo = 2.3 / (sin(el_rad) + 0.15);
        double multipath = 0.3 * sin(i * 1.7);
        double noise = 0.5 * ((double)((i * 12345) % 10000) / 10000.0 - 0.5);
        double rx_clock = 10000.0;

        meas[i].pseudorange = true_range + rx_clock + sat_clk_m
                              + iono + tropo + multipath + noise;
        meas[i].elevation = el_rad;
        meas[i].azimuth = az;
        meas[i].weight = 1.0;
        meas[i].used = 1;

        printf("PRN %2d: el=%5.1f deg az=%6.1f deg range=%10.3f m, PR=%10.3f m\n",
               i+1, sat_el_deg[i], sat_az_deg[i],
               true_range, meas[i].pseudorange);
    }

    printf("\n--- Position Solution ---\n");

    gnss_ecef_t init_guess = {0.0, 0.0, 0.0};

    gnss_bancroft_result_t bancroft;
    int rc = gnss_bancroft_solve(meas, n_sats, &bancroft);
    if (rc == 0 && bancroft.valid) {
        printf("Bancroft initial: (%.1f, %.1f, %.1f), bias=%.1f m\n",
               bancroft.pos1.x, bancroft.pos1.y, bancroft.pos1.z, bancroft.bias1);
        init_guess = bancroft.pos1;
    }

    gnss_pvt_solution_t sol;
    rc = gnss_ls_position_solve(meas, n_sats, init_guess, 15, 1e-4, &sol);

    if (rc == 0 && sol.valid) {
        gnss_lla_t sol_lla = gnss_ecef_to_lla(sol.pos);

        double horz_err = sqrt(pow((sol_lla.lat - true_lla.lat) *
                                    gnss_rad_curvature_meridian(true_lla.lat), 2.0)
                              + pow((sol_lla.lon - true_lla.lon) *
                                     gnss_rad_curvature_prime(true_lla.lat)
                                     * cos(true_lla.lat), 2.0));
        double vert_err = sol_lla.alt - true_lla.alt;

        printf("\n=== Solution Converged ===\n");
        printf("Position (ECEF): (%.3f, %.3f, %.3f) m\n",
               sol.pos.x, sol.pos.y, sol.pos.z);
        printf("Position (LLA):  Lat=%.8f deg, Lon=%.8f deg, Alt=%.3f m\n",
               sol_lla.lat * 180.0/M_PI, sol_lla.lon * 180.0/M_PI, sol_lla.alt);
        printf("Clock bias:     %.6f m (%.3f us)\n",
               sol.clock_bias, sol.clock_bias / GNSS_C_LIGHT * 1e6);
        printf("Iterations:     %d\n", sol.iterations);
        printf("RMS residual:  %.4f m\n", sol.rms_residual);
        printf("\nDOP Values:\n");
        printf("  GDOP = %.2f   PDOP = %.2f\n", sol.gdop, sol.pdop);
        printf("  HDOP = %.2f   VDOP = %.2f   TDOP = %.2f\n",
               sol.hdop, sol.vdop, sol.tdop);
        printf("\nPosition Errors:\n");
        printf("  Horizontal: %.3f m\n", horz_err);
        printf("  Vertical:   %.3f m\n", vert_err);

        gnss_dop_rating_t rating = gnss_dop_rate_pdop(sol.pdop);
        printf("\nDOP Rating: %s\n", gnss_dop_rating_string(rating));
    } else {
        printf("Position solution FAILED (rc=%d, valid=%d)\n", rc, sol.valid);
    }

    printf("\n=== Example Complete ===\n");
    return 0;
}
