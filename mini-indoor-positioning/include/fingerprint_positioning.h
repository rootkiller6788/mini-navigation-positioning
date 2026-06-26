/**
 * @file fingerprint_positioning.h
 * @brief WiFi/BLE/Magnetic fingerprint-based indoor positioning
 *
 * Knowledge Coverage:
 *   L1 - Definitions: fingerprint, radio map, AP/BSSID, RSSI vector,
 *        magnetic field vector, survey point
 *   L2 - Core Concepts: nearest-neighbor matching, k-NN, weighted k-NN,
 *        probabilistic fingerprinting (Bayesian), clustering
 *   L5 - Algorithms: Euclidean distance matching, Manhattan distance,
 *        cosine similarity, histogram matching, k-NN classification,
 *        Gaussian likelihood with covariance
 *   L6 - Canonical Problems: WiFi RSSI fingerprint positioning in
 *        an office building (RADAR system, Microsoft Research)
 *
 * Reference: Bahl & Padmanabhan, "RADAR: An In-Building RF-based User
 * Location and Tracking System," IEEE INFOCOM, 2000.
 *            Youssef & Agrawala, "The Horus WLAN location determination
 * system," ACM MobiSys, 2005.
 *
 * Course Alignment:
 *   - MIT 6.450 (Digital Communications)
 *   - Berkeley EE123 (DSP) — pattern matching
 *   - Stanford EE359 (Wireless)
 *   - 清华 信号与系统 — 模式识别基础
 */

#ifndef FINGERPRINT_POSITIONING_H
#define FINGERPRINT_POSITIONING_H

#include "indoor_positioning.h"

/* ============================================================================
 * L1 - Definitions: Fingerprint Data Structures
 * ============================================================================ */

#define FP_MAX_APS          32   /**< Maximum access points per fingerprint */
#define FP_MAX_SURVEY_PTS   256  /**< Maximum survey points in radio map */
#define FP_MAX_BSSID_LEN    20   /**< BSSID string length */
#define FP_MAX_SSID_LEN     32   /**< SSID string length */

/** MAC address (BSSID) representation for WiFi access points */
typedef struct {
    uint8_t octets[6];
} mac_address_t;

/**
 * @brief Access point (AP) information
 *
 * L1 Definition: A WiFi access point is a fixed infrastructure node
 * that transmits beacon frames at regular intervals. Each AP is
 * uniquely identified by its BSSID (MAC address).
 */
typedef struct {
    mac_address_t bssid;         /**< AP MAC address */
    char          ssid[FP_MAX_SSID_LEN]; /**< Network name */
    int           channel;       /**< WiFi channel (1-14 for 2.4GHz, 36-165 for 5GHz) */
    double        tx_power_dbm;  /**< Nominal transmit power in dBm */
    double        frequency_mhz; /**< Center frequency in MHz */
    position3d_t  position;      /**< Known or estimated AP position */
    int           is_anchor;     /**< 1 if position is surveyed, 0 if unknown */
} access_point_t;

/**
 * @brief RSSI fingerprint: RSSI readings from multiple APs at one location
 *
 * L1 Definition: A fingerprint is a vector of RSSI values
 * F = [RSSI_1, RSSI_2, ..., RSSI_N] measured at a known position.
 * The set of all fingerprints forms the radio map.
 *
 * Each entry records the RSSI from one access point. Missing APs
 * are marked with FP_RSSI_INVALID.
 */
#define FP_RSSI_INVALID (-200.0) /**< Sentinel for missing AP */

typedef struct {
    mac_address_t bssid;       /**< AP MAC address */
    double        rssi_mean;   /**< Mean RSSI in dBm over survey window */
    double        rssi_std;    /**< RSSI standard deviation in dB */
    int           n_samples;   /**< Number of survey samples collected */
} rssi_reading_t;

/**
 * @brief A single survey point (fingerprint at known location)
 *
 * L1 Definition: A survey point is one reference position in the
 * radio map, containing the RSSI vector from all visible APs at
 * that location, along with metadata.
 */
typedef struct {
    position3d_t  position;    /**< Surveyed position (ground truth) */
    rssi_reading_t readings[FP_MAX_APS]; /**< RSSI readings per AP */
    int           n_aps;       /**< Number of APs seen at this point */
    uint64_t      timestamp;   /**< Survey timestamp in microseconds */
    int           floor_id;    /**< Floor identifier */
    char          label[32];   /**< Human-readable label (e.g., "Room 301") */
} survey_point_t;

/**
 * @brief Magnetic field fingerprint
 *
 * L1 Definition: The Earth's magnetic field is distorted by
 * steel structures in buildings, creating unique vector patterns
 * at each location. This enables magnetic-field-based positioning.
 */
typedef struct {
    double mag_x;     /**< Magnetic field X component in microtesla */
    double mag_y;     /**< Magnetic field Y component in microtesla */
    double mag_z;     /**< Magnetic field Z component in microtesla */
    double magnitude; /**< Total magnetic field magnitude in uT */
    double declination_deg; /**< Declination angle in degrees */
    double inclination_deg; /**< Inclination angle in degrees */
} magnetic_vector_t;

/**
 * @brief Radio map: complete fingerprint database for indoor positioning
 *
 * L2 Concept: The radio map is the training dataset for fingerprint-based
 * positioning. It is created during the offline survey phase and used
 * during the online positioning phase.
 */
typedef struct {
    survey_point_t  points[FP_MAX_SURVEY_PTS]; /**< All survey points */
    int             n_points;                    /**< Number of survey points */
    access_point_t  known_aps[FP_MAX_APS];       /**< Known AP metadata */
    int             n_known_aps;                  /**< Number of known APs */
    double          survey_area_x_min;            /**< Bounding box of survey area */
    double          survey_area_x_max;
    double          survey_area_y_min;
    double          survey_area_y_max;
    int             n_floors;                     /**< Number of floors surveyed */
    uint64_t        survey_date;                  /**< Survey completion timestamp */
} radio_map_t;

/* ============================================================================
 * L5 - Algorithms: Fingerprint Matching
 * ============================================================================ */

/**
 * @brief Match an observed RSSI vector against the radio map using
 * nearest-neighbor (1-NN) in signal space.
 *
 * The "signal distance" between two RSSI vectors:
 *   d(F_obs, F_ref) = sqrt( sum_i (RSSI_obs,i - RSSI_ref,i)^2 )
 *
 * Only APs present in both vectors are compared.
 *
 * @param observed_rssi Array of observed RSSI values
 * @param observed_bssids Corresponding BSSID MAC addresses
 * @param n_observed Number of observed APs
 * @param radio_map The trained radio map
 * @param[out] result Estimated position
 * @return 0 on success, -1 if radio map is empty or no match
 *
 * L5: Nearest-neighbor classification in signal space.
 * Reference: Bahl & Padmanabhan (2000) — RADAR system.
 *
 * Complexity: O(M * K) where M = survey points, K = APs per point
 */
int fingerprint_match_nn(const double *observed_rssi,
                         const mac_address_t *observed_bssids,
                         int n_observed,
                         const radio_map_t *radio_map,
                         position3d_t *result);

/**
 * @brief k-Nearest-Neighbors (k-NN) fingerprint matching
 *
 * Finds the k closest survey points in signal space and returns
 * the centroid of their positions as the estimate.
 *
 * @param observed_rssi Observed RSSI array
 * @param observed_bssids Corresponding BSSID addresses
 * @param n_observed Number of observed APs
 * @param radio_map The radio map
 * @param k Number of neighbors (typical: 3-5)
 * @param[out] result Estimated position (average of k nearest)
 * @return 0 on success, -1 on error
 *
 * L5: k-NN classification in signal space.
 * Complexity: O(M * K + M log M) = O(M * (K + log M))
 */
int fingerprint_match_knn(const double *observed_rssi,
                          const mac_address_t *observed_bssids,
                          int n_observed,
                          const radio_map_t *radio_map,
                          int k,
                          position3d_t *result);

/**
 * @brief Weighted k-NN fingerprint matching
 *
 * Weights each neighbor inversely proportional to its signal distance:
 *   w_i = 1 / (d_i + epsilon)
 *   pos_est = sum(w_i * pos_i) / sum(w_i)
 *
 * This gives more influence to closer matches in signal space,
 * improving accuracy over plain k-NN.
 *
 * @param observed_rssi Observed RSSI array
 * @param observed_bssids Corresponding BSSID addresses
 * @param n_observed Number of observed APs
 * @param radio_map The radio map
 * @param k Number of neighbors
 * @param epsilon Small constant to avoid division by zero
 * @param[out] result Weighted position estimate
 * @return 0 on success
 *
 * L5: Weighted k-NN with inverse-distance weighting.
 * Complexity: O(M * (K + log M))
 */
int fingerprint_match_wknn(const double *observed_rssi,
                           const mac_address_t *observed_bssids,
                           int n_observed,
                           const radio_map_t *radio_map,
                           int k,
                           double epsilon,
                           position3d_t *result);

/**
 * @brief Probabilistic fingerprinting: Maximum Likelihood estimation
 * with Gaussian RSSI distribution per AP per survey point.
 *
 * For each survey point, computes the likelihood of the observed
 * RSSI vector assuming independent Gaussian distributions:
 *   P(F_obs | pos_i) = prod_j N(RSSI_obs,j; mu_ij, sigma_ij^2)
 *   pos_est = argmax_i P(F_obs | pos_i)
 *
 * @param observed_rssi Observed RSSI array
 * @param observed_bssids Corresponding BSSID addresses
 * @param n_observed Number of observed APs
 * @param radio_map Radio map with mean/std per AP per point
 * @param[out] result Maximum likelihood position
 * @return 0 on success
 *
 * L5: Probabilistic fingerprinting (Horus system).
 * Reference: Youssef & Agrawala (2005).
 *
 * Complexity: O(M * K)
 */
int fingerprint_match_probabilistic(const double *observed_rssi,
                                    const mac_address_t *observed_bssids,
                                    int n_observed,
                                    const radio_map_t *radio_map,
                                    position3d_t *result);

/**
 * @brief Compute signal distance between two RSSI vectors (Euclidean in RSSI space)
 *
 * @param rssi_a First RSSI vector
 * @param rssi_b Second RSSI vector
 * @param n Number of APs common to both
 * @return Signal distance (RMS of RSSI differences)
 *
 * Complexity: O(n)
 */
double signal_distance_euclidean(const double *rssi_a, const double *rssi_b, int n);

/**
 * @brief Compute signal distance using Manhattan metric
 *
 * d = sum_i |RSSI_a,i - RSSI_b,i|
 *
 * @param rssi_a First RSSI vector
 * @param rssi_b Second RSSI vector
 * @param n Number of APs common to both
 * @return Manhattan signal distance
 *
 * Complexity: O(n)
 */
double signal_distance_manhattan(const double *rssi_a, const double *rssi_b, int n);

/**
 * @brief Compute cosine similarity between two RSSI vectors
 *
 * sim = (A · B) / (||A|| * ||B||)
 *
 * Returns a value in [-1, 1], where 1 means identical direction.
 * For positioning, higher similarity = closer match.
 *
 * @param rssi_a First RSSI vector
 * @param rssi_b Second RSSI vector
 * @param n Number of components
 * @return Cosine similarity [-1, 1]
 *
 * Complexity: O(n)
 */
double signal_similarity_cosine(const double *rssi_a, const double *rssi_b, int n);

/* ============================================================================
 * L2 - Core Concepts: Radio Map Construction
 * ============================================================================ */

/**
 * @brief Initialize an empty radio map
 *
 * @param[out] map Radio map to initialize
 */
void radio_map_init(radio_map_t *map);

/**
 * @brief Add a survey point to the radio map
 *
 * @param map Radio map to modify
 * @param point Survey point to add
 * @return 0 on success, -1 if map is full
 */
int radio_map_add_point(radio_map_t *map, const survey_point_t *point);

/**
 * @brief Register a known access point in the radio map
 *
 * @param map Radio map
 * @param ap Access point information
 * @return 0 on success, -1 if AP table is full
 */
int radio_map_register_ap(radio_map_t *map, const access_point_t *ap);

/**
 * @brief Get survey points on a specific floor
 *
 * @param map Radio map
 * @param floor_id Target floor
 * @param[out] indices Array of indices of matching survey points
 * @param max_indices Maximum number of indices to return
 * @return Number of survey points found on this floor
 */
int radio_map_get_floor_points(const radio_map_t *map, int floor_id,
                               int *indices, int max_indices);

/**
 * @brief Compute the inter-point distance matrix for a radio map
 *
 * Used for analyzing spatial coverage and clustering survey points.
 *
 * @param map Radio map
 * @param[out] dist_matrix Pre-allocated NxN matrix (N = map->n_points)
 * @return 0 on success
 */
int radio_map_spatial_distances(const radio_map_t *map, double *dist_matrix);

/**
 * @brief Estimate floor from observed RSSI using vertical AP discrimination
 *
 * Different floors have different AP visibility patterns. This function
 * estimates which floor the user is on before 2D positioning.
 *
 * @param observed_rssi Observed RSSI values
 * @param observed_bssids Corresponding BSSID addresses
 * @param n_observed Number of observed APs
 * @param map Radio map
 * @return Estimated floor ID, or -1 on failure
 */
int estimate_floor_from_rssi(const double *observed_rssi,
                             const mac_address_t *observed_bssids,
                             int n_observed,
                             const radio_map_t *map);

/* ============================================================================
 * L6 - Canonical Problem: Magnetic Field Positioning
 * ============================================================================ */

/**
 * @brief Magnetic field fingerprint data point
 */
typedef struct {
    position3d_t     position;
    magnetic_vector_t mag;
    uint64_t         timestamp;
} mag_survey_point_t;

/**
 * @brief Estimate position using magnetic field matching (MAGCOM approach)
 *
 * Compares observed magnetic vector against database using
 * vector magnitude and angular differences.
 *
 * @param observed Measured magnetic vector
 * @param database Magnetic survey database
 * @param n_db Number of database entries
 * @param[out] result Best-matching position
 * @return 0 on success
 *
 * Reference: Li et al., "How feasible is the use of magnetic field
 * alone for indoor positioning?" IEEE IPIN, 2012.
 */
int magnetic_field_match(const magnetic_vector_t *observed,
                         const mag_survey_point_t *database,
                         int n_db,
                         position3d_t *result);

#endif /* FINGERPRINT_POSITIONING_H */
