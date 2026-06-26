/**
 * @file example_rssi_trilateration.c
 * @brief End-to-end example: RSSI-based trilateration in an office floor
 *
 * L6 Canonical Problem: WiFi RSSI positioning using trilateration.
 * Simulates a 20m × 15m office with 4 WiFi access points.
 * Converts RSSI measurements to distances using a path loss model,
 * then solves for user position using linear least-squares trilateration.
 *
 * Reference: Bahl & Padmanabhan (2000) — RADAR system.
 */

#include <stdio.h>
#include <math.h>
#include "../include/indoor_positioning.h"

int main(void) {
    printf("╔══════════════════════════════════════════════════════════╗\n");
    printf("║  RSSI-Based Indoor Positioning — Office Floorplan       ║\n");
    printf("╚══════════════════════════════════════════════════════════╝\n\n");

    /* Office floorplan: 20m × 15m
     * APs at corners: (0,0), (20,0), (0,15), (20,15) */
    position2d_t aps[4] = {
        {0.0, 0.0},    /* AP1: Southwest corner */
        {20.0, 0.0},   /* AP2: Southeast corner */
        {0.0, 15.0},   /* AP3: Northwest corner */
        {20.0, 15.0}   /* AP4: Northeast corner */
    };

    /* Path loss model for office environment (2.4 GHz WiFi)
     * n=3.0 typical for office with partitions */
    path_loss_model_t model = {
        .rssi_at_1m = -30.0,    /* RSSI at 1m reference */
        .path_loss_exp = 3.0,   /* Office path loss exponent */
        .shadow_std = 4.0,      /* Shadow fading std dev */
        .frequency_mhz = 2400.0 /* 2.4 GHz WiFi */
    };

    printf("Office Layout: %.0fm × %.0fm\n", 20.0, 15.0);
    printf("WiFi APs at: (0,0), (20,0), (0,15), (20,15)\n");
    printf("Path Loss Model: n=%.1f, RSSI_1m=%.0f dBm\n\n",
           model.path_loss_exp, model.rssi_at_1m);

    /* Simulate user at position (7.0, 4.0) — near cubicle area */
    position2d_t true_pos = {7.0, 4.0};
    printf("True User Position: (%.1f, %.1f)\n\n", true_pos.x, true_pos.y);

    /* Compute true distances and expected RSSI at each AP */
    printf("AP   | Distance (m) | RSSI (dBm)\n");
    printf("-----+-------------+------------\n");

    double true_dists[4];
    double measured_rssi[4];

    for (int i = 0; i < 4; i++) {
        true_dists[i] = distance_2d(true_pos, aps[i]);
        double rssi_clean = distance_to_rssi(true_dists[i], &model);
        /* Add shadow fading noise */
        double noise = ((i * 7 + 3) % 17 - 8.0) * model.shadow_std / 4.0;
        measured_rssi[i] = rssi_clean + noise;
        printf("AP%1d  | %11.2f | %10.1f\n", i+1, true_dists[i], measured_rssi[i]);
    }

    /* Convert RSSI to estimated distances */
    printf("\n--- Converting RSSI → Distance Estimates ---\n");
    printf("AP   | RSSI (dBm) | Est Dist (m) | True Dist (m) | Error (m)\n");
    printf("-----+------------+--------------+---------------+----------\n");

    double est_dists[4];
    for (int i = 0; i < 4; i++) {
        est_dists[i] = rssi_to_distance(measured_rssi[i], &model);
        double error = fabs(est_dists[i] - true_dists[i]);
        printf("AP%1d  | %10.1f | %12.2f | %13.2f | %9.2f\n",
               i+1, measured_rssi[i], est_dists[i], true_dists[i], error);
    }

    /* Solve trilateration */
    printf("\n--- Solving Trilateration ---\n");
    position2d_t est_pos;
    int ret = trilateration_2d(aps, est_dists, 4, &est_pos);
    double error_2d = 99.0;

    if (ret == 0) {
        error_2d = distance_2d(true_pos, est_pos);
        printf("Estimated Position: (%.2f, %.2f)\n", est_pos.x, est_pos.y);
        printf("Positioning Error: %.2f meters\n", error_2d);
        printf("\nResult: Trilateration succeeded with %.2fm accuracy.\n", error_2d);
    } else {
        printf("Trilateration failed (degenerate geometry or bad RSSI).\n");
    }

    /* Compare with weighted centroid */
    printf("\n--- Weighted Centroid (Baseline) ---\n");
    position2d_t centroid;
    weighted_centroid_2d(aps, est_dists, 4, &centroid);
    double centroid_error = distance_2d(true_pos, centroid);
    printf("Centroid Position: (%.2f, %.2f)\n", centroid.x, centroid.y);
    printf("Centroid Error: %.2f meters\n", centroid_error);
    if (ret == 0) {
        printf("Improvement over centroid: %.1f%%\n",
               (1.0 - error_2d / centroid_error) * 100.0);
    }

    printf("\n════════════════════════════════════════════════\n");
    printf("  L6: WiFi RSSI Positioning Complete\n");
    printf("  Trilateration accuracy: %.2f m in %.0fm×%.0fm office\n",
           error_2d, 20.0, 15.0);
    printf("════════════════════════════════════════════════\n");

    return 0;
}
