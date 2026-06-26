/**
 * @file example_fingerprint_navigation.c
 * @brief End-to-end example: WiFi fingerprint positioning in a museum
 *
 * L7 Application: Indoor navigation for museum visitors using smartphone
 * WiFi RSSI fingerprinting. Demonstrates radio map construction, floor
 * estimation, and k-NN positioning.
 *
 * Simulated museum: 3 floors with gallery rooms, 6 WiFi APs.
 * Surveys 24 reference points, performs online positioning.
 *
 * Reference: Youssef & Agrawala (2005) — Horus system
 *            Real-world: Google Indoor Maps, Apple Indoor Survey
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../include/indoor_positioning.h"
#include "../include/fingerprint_positioning.h"
#include "../include/positioning_error.h"

/* ================================================================
 * Helper: create a MAC address from an ID
 * ================================================================ */
static mac_address_t make_bssid(int id) {
    mac_address_t mac = {0};
    mac.octets[5] = (uint8_t)(id & 0xFF);
    mac.octets[4] = (uint8_t)((id >> 8) & 0xFF);
    return mac;
}

/* ================================================================
 * Helper: add a survey point to the radio map
 * ================================================================ */
static void add_point(radio_map_t *map, double x, double y, double z,
                      int floor, const char *label,
                      const double *rssi_vals, const mac_address_t *bssids,
                      int n_aps) {
    survey_point_t sp;
    memset(&sp, 0, sizeof(sp));
    sp.position.x = x;
    sp.position.y = y;
    sp.position.z = z;
    sp.floor_id = floor;
    strncpy(sp.label, label, 31);
    sp.label[31] = '\0';
    for (int i = 0; i < n_aps && i < FP_MAX_APS; i++) {
        sp.readings[i].bssid = bssids[i];
        sp.readings[i].rssi_mean = rssi_vals[i];
        sp.readings[i].rssi_std = 3.0;  /* Survey std */
        sp.readings[i].n_samples = 30;
    }
    sp.n_aps = n_aps;
    radio_map_add_point(map, &sp);
}

int main(void) {
    printf("╔══════════════════════════════════════════════════════════╗\n");
    printf("║  Museum Indoor Navigation — WiFi Fingerprint System      ║\n");
    printf("╚══════════════════════════════════════════════════════════╝\n\n");

    /* Museum: 40m × 30m, 3 floors */
    printf("Museum Layout: 40m × 30m, 3 floors\n");
    printf("WiFi APs: 6 (2 per floor)\n\n");

    /* Define 6 APs */
    mac_address_t ap_macs[6];
    for (int i = 0; i < 6; i++) ap_macs[i] = make_bssid(i + 1);

    /* Build radio map */
    radio_map_t rmap;
    radio_map_init(&rmap);

    /* Register APs */
    for (int i = 0; i < 6; i++) {
        access_point_t ap;
        memset(&ap, 0, sizeof(ap));
        ap.bssid = ap_macs[i];
        snprintf(ap.ssid, FP_MAX_SSID_LEN, "Museum-WiFi-%d", i+1);
        ap.tx_power_dbm = 20.0;
        ap.frequency_mhz = 2400.0;
        radio_map_register_ap(&rmap, &ap);
    }

    /* Survey points: 3 floors × 8 points = 24 total
     * Each with RSSI readings from visible APs (simulated) */
    printf("--- Building Radio Map (Offline Survey Phase) ---\n");
    printf("Floor | Room              | Visible APs\n");
    printf("------+-------------------+------------\n");

    /* Floor 1: Lobby, Gift Shop, Egyptian Gallery, Greek Gallery */
    double fp1_rssi[6][6] = {
        { -40, -60, -75, -80, -85, -90 },  /* Near AP1 */
        { -55, -45, -70, -78, -82, -88 },  /* Near AP2 */
        { -70, -72, -45, -60, -80, -85 },  /* Near AP3 */
        { -68, -70, -58, -42, -78, -82 },  /* Near AP4 */
        { -75, -78, -48, -65, -82, -87 },  /* Mid */
        { -62, -65, -68, -55, -80, -85 },  /* Mid */
        { -80, -82, -52, -50, -78, -83 },  /* Corridor */
        { -72, -75, -45, -62, -82, -88 },  /* Corridor */
    };

    const char *labels_f1[] = {
        "F1-Lobby", "F1-GiftShop", "F1-Egyptian", "F1-Greek",
        "F1-Corridor-W", "F1-Corridor-E", "F1-Atrium-N", "F1-Atrium-S"
    };

    for (int i = 0; i < 8; i++) {
        double x = (i % 4) * 10.0 + 5.0;
        double y = (i < 4) ? 20.0 : 5.0;
        /* Visible APs: AP1, AP2 (floor 1) + AP3, AP4 (floor 2, weaker) */
        mac_address_t visible[4] = {ap_macs[0], ap_macs[1], ap_macs[2], ap_macs[3]};
        double rssi[4] = {fp1_rssi[i][0], fp1_rssi[i][1], fp1_rssi[i][2], fp1_rssi[i][3]};
        add_point(&rmap, x, y, 0.0, 1, labels_f1[i], rssi, visible, 4);
        printf("  F1  | %-18s| AP1-AP4\n", labels_f1[i]);
    }

    /* Floor 2: Renaissance, Modern, Impressionist, Sculpture */
    double fp2_rssi[8][4] = {
        { -72, -75, -42, -58 }, { -78, -80, -38, -55 },
        { -70, -73, -55, -40 }, { -68, -72, -50, -45 },
        { -65, -70, -48, -52 }, { -75, -78, -42, -55 },
        { -80, -82, -38, -58 }, { -73, -76, -45, -50 },
    };

    const char *labels_f2[] = {
        "F2-Renaissance", "F2-Modern", "F2-Impressionist", "F2-Sculpture",
        "F2-Corridor-N", "F2-Corridor-S", "F2-Balcony-E", "F2-Balcony-W"
    };

    for (int i = 0; i < 8; i++) {
        double x = (i % 4) * 10.0 + 5.0;
        double y = (i < 4) ? 20.0 : 5.0;
        /* Visible APs: AP3, AP4 (floor 2) + AP5, AP6 (floor 3, weaker) */
        mac_address_t visible[4] = {ap_macs[2], ap_macs[3], ap_macs[4], ap_macs[5]};
        double rssi[4] = {fp2_rssi[i][0], fp2_rssi[i][1], fp2_rssi[i][2], fp2_rssi[i][3]};
        add_point(&rmap, x, y, 5.0, 2, labels_f2[i], rssi, visible, 4);
        printf("  F2  | %-18s| AP3-AP6\n", labels_f2[i]);
    }

    printf("\nRadio Map: %d survey points, %d known APs\n",
           rmap.n_points, rmap.n_known_aps);

    /* ================================================================
     * Online Positioning Phase
     * ================================================================ */
    printf("\n--- Online Positioning Phase ---\n");
    printf("Scenario: Visitor at F2-Modern Gallery (true: 15, 20, F2)\n\n");

    /* Simulate observed RSSI at visitor's location (near Modern gallery) */
    double obs_rssi[4] = {-77.0, -79.0, -40.0, -56.0};
    mac_address_t obs_bssids[4];
    for (int i = 0; i < 4; i++) obs_bssids[i] = ap_macs[i+2];

    /* Step 1: Estimate floor */
    int estimated_floor = estimate_floor_from_rssi(obs_rssi, obs_bssids, 4, &rmap);
    printf("Step 1 - Floor Estimation: Floor %d\n", estimated_floor);

    /* Step 2: k-NN positioning (k=4) */
    position3d_t pos_nn, pos_knn, pos_wknn;

    fingerprint_match_nn(obs_rssi, obs_bssids, 4, &rmap, &pos_nn);
    printf("Step 2a - 1-NN:  (%.1f, %.1f, Floor %d)\n",
           pos_nn.x, pos_nn.y, estimated_floor);

    fingerprint_match_knn(obs_rssi, obs_bssids, 4, &rmap, 4, &pos_knn);
    printf("Step 2b - 4-NN:  (%.1f, %.1f, Floor %d)\n",
           pos_knn.x, pos_knn.y, estimated_floor);

    fingerprint_match_wknn(obs_rssi, obs_bssids, 4, &rmap, 4, 0.01, &pos_wknn);
    printf("Step 2c - WKNN:  (%.1f, %.1f, Floor %d)\n",
           pos_wknn.x, pos_wknn.y, estimated_floor);

    /* Accuracy assessment */
    double true_x = 15.0, true_y = 20.0;
    double err_nn = sqrt((pos_nn.x-true_x)*(pos_nn.x-true_x)
                       + (pos_nn.y-true_y)*(pos_nn.y-true_y));
    double err_knn = sqrt((pos_knn.x-true_x)*(pos_knn.x-true_x)
                         + (pos_knn.y-true_y)*(pos_knn.y-true_y));
    double err_wknn = sqrt((pos_wknn.x-true_x)*(pos_wknn.x-true_x)
                          + (pos_wknn.y-true_y)*(pos_wknn.y-true_y));

    printf("\nAccuracy (vs true position %.0f, %.0f):\n", true_x, true_y);
    printf("  1-NN error:  %.2f m\n", err_nn);
    printf("  4-NN error:  %.2f m\n", err_knn);
    printf("  WKNN error:  %.2f m\n", err_wknn);

    /* Display nearby rooms */
    printf("\nNearby Points:\n");
    printf("  Room              | Distance (sig) | Position\n");
    printf("  -------------------+---------------+----------\n");

    /* Top 3 nearest survey points */
    for (int p = 0; p < 3; p++) {
        /* Find p-th nearest */
        int best_idx = -1;
        double best_dist = 1e9;
        static int visited[FP_MAX_SURVEY_PTS] = {0};

        for (int i = 0; i < rmap.n_points; i++) {
            if (visited[i]) continue;
            /* Simple distance in signal space to first AP */
            double d = 0.0;
            int matches = 0;
            for (int j = 0; j < rmap.points[i].n_aps && j < 4; j++) {
                if (memcmp(&rmap.points[i].readings[j].bssid, &obs_bssids[0],
                           sizeof(mac_address_t)) == 0) {
                    d += fabs(obs_rssi[0] - rmap.points[i].readings[j].rssi_mean);
                    matches++;
                    break;
                }
            }
            if (matches > 0 && d < best_dist) {
                best_dist = d;
                best_idx = i;
            }
        }
        if (best_idx >= 0) {
            printf("  %-18s | %13.1f | (%.0f,%.0f, F%d)\n",
                   rmap.points[best_idx].label, best_dist,
                   rmap.points[best_idx].position.x,
                   rmap.points[best_idx].position.y,
                   rmap.points[best_idx].floor_id);
            visited[best_idx] = 1;
        }
    }

    printf("\n══════════════════════════════════════════════════════════\n");
    printf("  L7: Museum Indoor Navigation Complete\n");
    printf("  WiFi fingerprinting: %.1fm accuracy (WKNN)\n", err_wknn);
    printf("  Application: Museum visitor guidance system\n");
    printf("══════════════════════════════════════════════════════════\n");

    return 0;
}
