/**
 * @example example_gnss_pvt.c
 * @brief Standalone GNSS Position Solution Example
 *
 * Demonstrates GNSS WLS positioning with satellite geometry
 * and coordinate conversions.
 * L6 Canonical Problem: GNSS standalone positioning.
 */

#include "nav_gnss.h"
#include "nav_kalman.h"
#include <stdio.h>
#include <math.h>
#include <string.h>

int main(void) {
    printf("=== GNSS PVT Example ===\n\n");
    printf("Receiver at: 45 deg N, 120 deg E, 500m altitude\n");
    nav_geodetic_t rx_true;
    rx_true.latitude  = nav_deg2rad(45.0);
    rx_true.longitude = nav_deg2rad(120.0);
    rx_true.altitude  = 500.0;
    NAV_PRECISION rx_ecef[3];
    nav_geodetic_to_ecef(rx_ecef, &rx_true);
    printf("ECEF: [%.1f, %.1f, %.1f] m\n", rx_ecef[0], rx_ecef[1], rx_ecef[2]);

    /* Test satellite: GPS SV 1 approximate position */
    printf("\nSatellite position computation:\n");
    nav_gnss_ephemeris_t eph;
    memset(&eph, 0, sizeof(eph));
    eph.sqrt_a = 5153.65;  /* ~26560 km semi-major axis */
    eph.ecc = 0.01;
    eph.i0 = nav_deg2rad(55.0);
    eph.omega0 = nav_deg2rad(120.0);
    eph.w = nav_deg2rad(45.0);
    eph.M0 = nav_deg2rad(30.0);
    eph.toe = 252000.0;
    NAV_PRECISION sat_pos[3], sat_clk;
    nav_gnss_sat_position(sat_pos, NULL, &sat_clk, &eph, 252000.0);
    printf("Sat ECEF: [%.1f, %.1f, %.1f] m\n", sat_pos[0], sat_pos[1], sat_pos[2]);
    printf("Sat clock correction: %.3f m\n", sat_clk);

    /* Azimuth and elevation */
    NAV_PRECISION az, el;
    nav_gnss_azel(&az, &el, rx_ecef, sat_pos);
    printf("Azimuth: %.1f deg, Elevation: %.1f deg\n",
           nav_rad2deg(az), nav_rad2deg(el));

    /* Coordinate round-trip test */
    printf("\nCoordinate round-trip test:\n");
    nav_geodetic_t test_pos = {nav_deg2rad(37.7749), nav_deg2rad(-122.4194), 10.0};
    NAV_PRECISION ecef[3];
    nav_geodetic_to_ecef(ecef, &test_pos);
    nav_geodetic_t geo2;
    nav_ecef_to_geodetic(&geo2, ecef);
    NAV_PRECISION ecef2[3];
    nav_geodetic_to_ecef(ecef2, &geo2);
    NAV_PRECISION err = sqrt((ecef[0]-ecef2[0])*(ecef[0]-ecef2[0]) +
                              (ecef[1]-ecef2[1])*(ecef[1]-ecef2[1]) +
                              (ecef[2]-ecef2[2])*(ecef[2]-ecef2[2]));
    printf("San Francisco: %.4f N, %.4f W, %.1f m\n",
           nav_rad2deg(geo2.latitude), nav_rad2deg(-geo2.longitude), geo2.altitude);
    printf("Round-trip error: %.6f m\n", err);

    /* Ionospheric delay example */
    printf("\nIonospheric delay (Klobuchar model):\n");
    NAV_PRECISION alpha[4] = {2.0e-8, 1.0e-8, -1.0e-7, 0.0};
    NAV_PRECISION beta[4]  = {1.0e5, 1.0e5, -1.0e5, 0.0};
    NAV_PRECISION iono_delay;
    nav_iono_klobuchar(&iono_delay, nav_deg2rad(45.0), nav_deg2rad(120.0),
                        nav_deg2rad(180.0), nav_deg2rad(60.0), 50400.0,
                        alpha, beta);
    printf("Ionospheric delay: %.3f m\n", iono_delay);

    printf("\nGNSS PVT example complete.\n");
    return 0;
}
