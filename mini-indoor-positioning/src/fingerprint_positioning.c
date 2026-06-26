/**
 * @file fingerprint_positioning.c
 * @brief Fingerprint-based indoor positioning algorithms
 *
 * Implements: radio map management, NN/k-NN/WKNN matching,
 * probabilistic fingerprinting, signal distance metrics,
 * magnetic field matching, floor estimation.
 */

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "../include/fingerprint_positioning.h"

/* ============================================================================
 * Signal Distance Metrics
 * ============================================================================ */

double signal_distance_euclidean(const double *rssi_a, const double *rssi_b, int n) {
    if (!rssi_a || !rssi_b || n <= 0) return 1e9;
    double sum_sq = 0.0;
    int valid = 0;
    for (int i = 0; i < n; i++) {
        /* Skip invalid readings */
        if (rssi_a[i] <= FP_RSSI_INVALID || rssi_b[i] <= FP_RSSI_INVALID) continue;
        double diff = rssi_a[i] - rssi_b[i];
        sum_sq += diff * diff;
        valid++;
    }
    if (valid == 0) return 1e9;
    return sqrt(sum_sq / valid);
}

double signal_distance_manhattan(const double *rssi_a, const double *rssi_b, int n) {
    if (!rssi_a || !rssi_b || n <= 0) return 1e9;
    double sum_abs = 0.0;
    int valid = 0;
    for (int i = 0; i < n; i++) {
        if (rssi_a[i] <= FP_RSSI_INVALID || rssi_b[i] <= FP_RSSI_INVALID) continue;
        sum_abs += fabs(rssi_a[i] - rssi_b[i]);
        valid++;
    }
    if (valid == 0) return 1e9;
    return sum_abs / valid;
}

double signal_similarity_cosine(const double *rssi_a, const double *rssi_b, int n) {
    if (!rssi_a || !rssi_b || n <= 0) return 0.0;
    double dot = 0.0, norm_a = 0.0, norm_b = 0.0;
    for (int i = 0; i < n; i++) {
        if (rssi_a[i] <= FP_RSSI_INVALID || rssi_b[i] <= FP_RSSI_INVALID) continue;
        dot += rssi_a[i] * rssi_b[i];
        norm_a += rssi_a[i] * rssi_a[i];
        norm_b += rssi_b[i] * rssi_b[i];
    }
    double denom = sqrt(norm_a) * sqrt(norm_b);
    if (denom < 1e-12) return 0.0;
    return dot / denom;
}

/* ============================================================================
 * Radio Map Management
 * ============================================================================ */

void radio_map_init(radio_map_t *map) {
    if (!map) return;
    memset(map, 0, sizeof(radio_map_t));
}

int radio_map_add_point(radio_map_t *map, const survey_point_t *point) {
    if (!map || !point) return -1;
    if (map->n_points >= FP_MAX_SURVEY_PTS) return -1;
    memcpy(&map->points[map->n_points], point, sizeof(survey_point_t));
    map->n_points++;

    /* Update bounding box */
    if (map->n_points == 1) {
        map->survey_area_x_min = point->position.x;
        map->survey_area_x_max = point->position.x;
        map->survey_area_y_min = point->position.y;
        map->survey_area_y_max = point->position.y;
    } else {
        if (point->position.x < map->survey_area_x_min) map->survey_area_x_min = point->position.x;
        if (point->position.x > map->survey_area_x_max) map->survey_area_x_max = point->position.x;
        if (point->position.y < map->survey_area_y_min) map->survey_area_y_min = point->position.y;
        if (point->position.y > map->survey_area_y_max) map->survey_area_y_max = point->position.y;
    }
    return 0;
}

int radio_map_register_ap(radio_map_t *map, const access_point_t *ap) {
    if (!map || !ap) return -1;
    if (map->n_known_aps >= FP_MAX_APS) return -1;
    memcpy(&map->known_aps[map->n_known_aps], ap, sizeof(access_point_t));
    map->n_known_aps++;
    return 0;
}

int radio_map_get_floor_points(const radio_map_t *map, int floor_id,
                               int *indices, int max_indices) {
    if (!map || !indices || max_indices <= 0) return 0;
    int count = 0;
    for (int i = 0; i < map->n_points && count < max_indices; i++) {
        if (map->points[i].floor_id == floor_id) {
            indices[count++] = i;
        }
    }
    return count;
}

int radio_map_spatial_distances(const radio_map_t *map, double *dist_matrix) {
    if (!map || !dist_matrix || map->n_points <= 0) return -1;
    int N = map->n_points;
    for (int i = 0; i < N; i++) {
        dist_matrix[i * N + i] = 0.0;
        for (int j = i + 1; j < N; j++) {
            double dx = map->points[i].position.x - map->points[j].position.x;
            double dy = map->points[i].position.y - map->points[j].position.y;
            double dz = map->points[i].position.z - map->points[j].position.z;
            double d = sqrt(dx*dx + dy*dy + dz*dz);
            dist_matrix[i * N + j] = d;
            dist_matrix[j * N + i] = d;
        }
    }
    return 0;
}

/* ============================================================================
 * Helper: compute RSSI vector distance to a survey point
 * ============================================================================ */

/**
 * @brief Compute Euclidean signal distance between observed RSSI vector
 *        and a single survey point, matching by BSSID.
 *
 * @param observed_rssi Array of observed RSSI values
 * @param observed_bssids Corresponding BSSID MACs
 * @param n_observed Number of observed APs
 * @param sp Survey point to compare against
 * @return Signal distance (lower = better match). Returns large number on no overlap.
 */
static double compute_signal_distance_to_point(const double *observed_rssi,
                                                const mac_address_t *observed_bssids,
                                                int n_observed,
                                                const survey_point_t *sp) {
    double sum_sq = 0.0;
    int matches = 0;

    for (int i = 0; i < n_observed; i++) {
        /* Find this AP in the survey point's readings */
        for (int j = 0; j < sp->n_aps; j++) {
            if (memcmp(&observed_bssids[i], &sp->readings[j].bssid,
                       sizeof(mac_address_t)) == 0) {
                double diff = observed_rssi[i] - sp->readings[j].rssi_mean;
                sum_sq += diff * diff;
                matches++;
                break;
            }
        }
    }
    if (matches == 0) return 1e9;
    return sqrt(sum_sq / matches);
}

/* ============================================================================
 * L5 - Nearest-Neighbor Fingerprint Matching
 * ============================================================================ */

int fingerprint_match_nn(const double *observed_rssi,
                         const mac_address_t *observed_bssids,
                         int n_observed,
                         const radio_map_t *radio_map,
                         position3d_t *result) {
    if (!observed_rssi || !observed_bssids || !radio_map || !result) return -1;
    if (n_observed <= 0) return -1;
    if (radio_map->n_points <= 0) return -1;

    double best_dist = 1e15;
    int best_idx = -1;

    for (int i = 0; i < radio_map->n_points; i++) {
        double d = compute_signal_distance_to_point(observed_rssi, observed_bssids,
                                                    n_observed, &radio_map->points[i]);
        if (d < best_dist) {
            best_dist = d;
            best_idx = i;
        }
    }

    if (best_idx < 0) return -1;
    *result = radio_map->points[best_idx].position;
    return 0;
}

/* ============================================================================
 * L5 - k-NN Fingerprint Matching
 * ============================================================================ */

/** Structure for sorting (distance, index) pairs */
typedef struct {
    double distance;
    int    index;
} dist_idx_pair_t;

static int compare_dist_idx(const void *a, const void *b) {
    double da = ((const dist_idx_pair_t *)a)->distance;
    double db = ((const dist_idx_pair_t *)b)->distance;
    if (da < db) return -1;
    if (da > db) return 1;
    return 0;
}

int fingerprint_match_knn(const double *observed_rssi,
                          const mac_address_t *observed_bssids,
                          int n_observed,
                          const radio_map_t *radio_map,
                          int k,
                          position3d_t *result) {
    if (!observed_rssi || !observed_bssids || !radio_map || !result) return -1;
    if (n_observed <= 0 || k <= 0) return -1;
    if (radio_map->n_points <= 0) return -1;

    /* Cap k to available points */
    if (k > radio_map->n_points) k = radio_map->n_points;

    dist_idx_pair_t *pairs = (dist_idx_pair_t *)malloc(
        radio_map->n_points * sizeof(dist_idx_pair_t));
    if (!pairs) return -1;

    for (int i = 0; i < radio_map->n_points; i++) {
        pairs[i].distance = compute_signal_distance_to_point(
            observed_rssi, observed_bssids, n_observed, &radio_map->points[i]);
        pairs[i].index = i;
    }

    qsort(pairs, radio_map->n_points, sizeof(dist_idx_pair_t), compare_dist_idx);

    /* Average the k nearest positions */
    double sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;
    for (int i = 0; i < k; i++) {
        int idx = pairs[i].index;
        sum_x += radio_map->points[idx].position.x;
        sum_y += radio_map->points[idx].position.y;
        sum_z += radio_map->points[idx].position.z;
    }

    result->x = sum_x / k;
    result->y = sum_y / k;
    result->z = sum_z / k;

    free(pairs);
    return 0;
}

/* ============================================================================
 * L5 - Weighted k-NN Fingerprint Matching
 * ============================================================================ */

int fingerprint_match_wknn(const double *observed_rssi,
                           const mac_address_t *observed_bssids,
                           int n_observed,
                           const radio_map_t *radio_map,
                           int k,
                           double epsilon,
                           position3d_t *result) {
    if (!observed_rssi || !observed_bssids || !radio_map || !result) return -1;
    if (n_observed <= 0 || k <= 0) return -1;
    if (radio_map->n_points <= 0) return -1;
    if (epsilon < 1e-12) epsilon = 1e-6;

    if (k > radio_map->n_points) k = radio_map->n_points;

    dist_idx_pair_t *pairs = (dist_idx_pair_t *)malloc(
        radio_map->n_points * sizeof(dist_idx_pair_t));
    if (!pairs) return -1;

    for (int i = 0; i < radio_map->n_points; i++) {
        pairs[i].distance = compute_signal_distance_to_point(
            observed_rssi, observed_bssids, n_observed, &radio_map->points[i]);
        pairs[i].index = i;
    }

    qsort(pairs, radio_map->n_points, sizeof(dist_idx_pair_t), compare_dist_idx);

    /* Weighted average: w_i = 1/(d_i + epsilon) */
    double sum_wx = 0.0, sum_wy = 0.0, sum_wz = 0.0, sum_w = 0.0;
    for (int i = 0; i < k; i++) {
        int idx = pairs[i].index;
        double w = 1.0 / (pairs[i].distance + epsilon);
        sum_wx += w * radio_map->points[idx].position.x;
        sum_wy += w * radio_map->points[idx].position.y;
        sum_wz += w * radio_map->points[idx].position.z;
        sum_w += w;
    }

    if (sum_w > 0.0) {
        result->x = sum_wx / sum_w;
        result->y = sum_wy / sum_w;
        result->z = sum_wz / sum_w;
    } else {
        /* Fall back to unweighted */
        double sx = 0, sy = 0, sz = 0;
        for (int i = 0; i < k; i++) {
            int idx = pairs[i].index;
            sx += radio_map->points[idx].position.x;
            sy += radio_map->points[idx].position.y;
            sz += radio_map->points[idx].position.z;
        }
        result->x = sx / k;
        result->y = sy / k;
        result->z = sz / k;
    }

    free(pairs);
    return 0;
}

/* ============================================================================
 * L5 - Probabilistic Fingerprint Matching (Horus-like)
 * ============================================================================ */

int fingerprint_match_probabilistic(const double *observed_rssi,
                                    const mac_address_t *observed_bssids,
                                    int n_observed,
                                    const radio_map_t *radio_map,
                                    position3d_t *result) {
    if (!observed_rssi || !observed_bssids || !radio_map || !result) return -1;
    if (n_observed <= 0 || radio_map->n_points <= 0) return -1;

    double best_log_likelihood = -1e15;
    int best_idx = -1;

    for (int i = 0; i < radio_map->n_points; i++) {
        const survey_point_t *sp = &radio_map->points[i];
        double log_likelihood = 0.0;
        int matches = 0;

        for (int j = 0; j < n_observed; j++) {
            /* Find matching AP in survey point */
            for (int a = 0; a < sp->n_aps; a++) {
                if (memcmp(&observed_bssids[j], &sp->readings[a].bssid,
                           sizeof(mac_address_t)) == 0) {
                    double diff = observed_rssi[j] - sp->readings[a].rssi_mean;
                    double sigma = sp->readings[a].rssi_std;
                    if (sigma < 1.0) sigma = 1.0;  /* Minimum std dev */
                    /* Log of Gaussian PDF: -0.5*log(2*pi*sigma^2) - 0.5*(diff/sigma)^2 */
                    log_likelihood += -0.5 * log(2.0 * M_PI * sigma * sigma)
                                    - 0.5 * (diff * diff) / (sigma * sigma);
                    matches++;
                    break;
                }
            }
            /* Missing APs: penalize with a uniform background likelihood */
            /* This is a simplified approach; Horus uses discrete distributions */
        }

        /* Penalize points with very few matching APs */
        if (matches < (n_observed / 4)) {
            log_likelihood -= 10.0;  /* Penalty for poor coverage */
        }

        if (log_likelihood > best_log_likelihood) {
            best_log_likelihood = log_likelihood;
            best_idx = i;
        }
    }

    if (best_idx < 0) return -1;
    *result = radio_map->points[best_idx].position;
    return 0;
}

/* ============================================================================
 * Floor Estimation from RSSI
 * ============================================================================ */

int estimate_floor_from_rssi(const double *observed_rssi,
                             const mac_address_t *observed_bssids,
                             int n_observed,
                             const radio_map_t *map) {
    if (!observed_rssi || !observed_bssids || !map) return -1;
    if (n_observed <= 0 || map->n_points <= 0) return -1;

    /* Count how many APs from each floor match the observed APs */
    /* Simple voting: for each survey point, compare signal distance */
    double best_dist = 1e15;
    int best_floor = -1;

    for (int i = 0; i < map->n_points; i++) {
        double d = compute_signal_distance_to_point(observed_rssi, observed_bssids,
                                                    n_observed, &map->points[i]);
        if (d < best_dist) {
            best_dist = d;
            best_floor = map->points[i].floor_id;
        }
    }

    return best_floor;
}

/* ============================================================================
 * L6 - Magnetic Field Matching (MAGCOM)
 * ============================================================================ */

int magnetic_field_match(const magnetic_vector_t *observed,
                         const mag_survey_point_t *database,
                         int n_db,
                         position3d_t *result) {
    if (!observed || !database || !result || n_db <= 0) return -1;

    double best_score = 1e15;
    int best_idx = -1;

    for (int i = 0; i < n_db; i++) {
        /* Use combined metric of magnitude difference and angular difference */
        double mag_diff = fabs(observed->magnitude - database[i].mag.magnitude);

        /* Angular difference via dot product of unit vectors */
        double obs_norm = observed->magnitude;
        double db_norm = database[i].mag.magnitude;
        if (obs_norm < 1e-6 || db_norm < 1e-6) continue;

        double dot = observed->mag_x * database[i].mag.mag_x
                   + observed->mag_y * database[i].mag.mag_y
                   + observed->mag_z * database[i].mag.mag_z;
        double cos_angle = dot / (obs_norm * db_norm);
        if (cos_angle > 1.0) cos_angle = 1.0;
        if (cos_angle < -1.0) cos_angle = -1.0;
        double angle_diff = acos(cos_angle);  /* radians */

        /* Combined score: weighted sum */
        double score = mag_diff * 0.5 + angle_diff * 50.0;  /* ~2.5m equivalent per rad */

        if (score < best_score) {
            best_score = score;
            best_idx = i;
        }
    }

    if (best_idx < 0) return -1;
    *result = database[best_idx].position;
    return 0;
}
