/* =========================================================================
 * example_dop_calc.c — DOP analysis and constellation selection
 *
 * Demonstrates: L1 (DOP definitions), L2 (satellite geometry impact),
 * L6 (constellation analysis), L7 (DOP prediction for mission planning).
 * ========================================================================= */

#include <stdio.h>
#include <math.h>
#include "gnss_common.h"
#include "gnss_dop.h"

int main(void) {
    printf("=== DOP Analysis & Constellation Selection Example ===\n\n");

    gnss_lla_t ref_lla = {0.7119, -1.2909, 30.0};
    gnss_ecef_t ref_pos = gnss_lla_to_ecef(ref_lla);

    printf("Reference: NYC Central Park (%.4fN, %.4fW, %.0fm)\n\n",
           ref_lla.lat*180.0/M_PI, -ref_lla.lon*180.0/M_PI, ref_lla.alt);

    /* Scenario 1: Good geometry — 8 satellites evenly distributed */
    printf("=== Scenario 1: Good Geometry (8 sats, even distribution) ===\n");

    gnss_ecef_t sats_good[8];
    double az_step = 2.0 * M_PI / 8.0;
    int i;
    for (i = 0; i < 8; i++) {
        double az = i * az_step;
        double el = 0.785398;
        double r = 26.0e6;
        gnss_enu_t enu;
        enu.e = r * cos(el) * sin(az);
        enu.n = r * cos(el) * cos(az);
        enu.u = r * sin(el);
        sats_good[i] = gnss_enu_to_ecef(enu, ref_lla);
    }

    gnss_dop_t dop_good;
    gnss_dop_from_sat_positions(sats_good, 8, ref_pos, &dop_good);
    printf("  GDOP=%.2f, PDOP=%.2f, HDOP=%.2f, VDOP=%.2f, TDOP=%.2f\n",
           dop_good.gdop, dop_good.pdop, dop_good.hdop, dop_good.vdop, dop_good.tdop);
    printf("  PDOP rating: %s\n", gnss_dop_rating_string(gnss_dop_rate_pdop(dop_good.pdop)));

    /* Scenario 2: Poor geometry — all satellites clustered */
    printf("\n=== Scenario 2: Poor Geometry (8 sats clustered in NW quadrant) ===\n");

    gnss_ecef_t sats_poor[8];
    for (i = 0; i < 8; i++) {
        double az = 0.5 + i * 0.1;
        double el = 0.3 + i * 0.02;
        double r = 26.0e6;
        gnss_enu_t enu;
        enu.e = r * cos(el) * sin(az);
        enu.n = r * cos(el) * cos(az);
        enu.u = r * sin(el);
        sats_poor[i] = gnss_enu_to_ecef(enu, ref_lla);
    }

    gnss_dop_t dop_poor;
    gnss_dop_from_sat_positions(sats_poor, 8, ref_pos, &dop_poor);
    printf("  GDOP=%.2f, PDOP=%.2f, HDOP=%.2f, VDOP=%.2f, TDOP=%.2f\n",
           dop_poor.gdop, dop_poor.pdop, dop_poor.hdop, dop_poor.vdop, dop_poor.tdop);
    printf("  PDOP rating: %s\n", gnss_dop_rating_string(gnss_dop_rate_pdop(dop_poor.pdop)));

    double uere = 4.0;
    printf("\n  Position accuracy comparison:\n");
    printf("    Good geometry: sigma_3D = %.1f m, sigma_H = %.1f m\n",
           dop_good.pdop * uere, dop_good.hdop * uere);
    printf("    Poor geometry: sigma_3D = %.1f m, sigma_H = %.1f m\n",
           dop_poor.pdop * uere, dop_poor.hdop * uere);

    /* Scenario 3: Elevation mask effect */
    printf("\n=== Scenario 3: Elevation Mask Effect ===\n");
    double masks[4] = {0.0, 5.0, 15.0, 30.0};
    int m;
    for (m = 0; m < 4; m++) {
        double mask_rad = masks[m] * M_PI / 180.0;
        gnss_geometry_stats_t stats;
        gnss_skyplot_entry_t entries[12];
        gnss_ecef_t sats_mix[12];
        for (i = 0; i < 12; i++) {
            double az = i * 2.0 * M_PI / 12.0;
            double el = 0.2 + (i % 4) * 0.2;
            gnss_enu_t enu;
            enu.e = 26e6 * cos(el) * sin(az);
            enu.n = 26e6 * cos(el) * cos(az);
            enu.u = 26e6 * sin(el);
            sats_mix[i] = gnss_enu_to_ecef(enu, ref_lla);
        }
        gnss_geometry_analyze(sats_mix, 12, ref_pos, mask_rad, entries, &stats);
        printf("  Mask %3.0f deg: %2d sats above, PDOP=%.2f\n",
               masks[m], stats.n_sats_above_mask,
               stats.dop.pdop);
    }

    printf("\n=== Example Complete ===\n");
    return 0;
}
