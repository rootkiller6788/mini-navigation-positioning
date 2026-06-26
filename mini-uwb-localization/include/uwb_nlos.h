/**
 * mini-uwb-localization: NLOS Detection and Mitigation
 *
 * Algorithmic identification of Non-Line-of-Sight (NLOS) conditions in
 * UWB ranging using Channel Impulse Response (CIR) statistical features,
 * and mitigation strategies to improve localization accuracy.
 *
 * Reference: Marano et al. (2010) "NLOS Identification and Mitigation
 *            for UWB Localization"
 * Reference: Guvenc et al. (2007) "NLOS Identification for UWB"
 * Reference: Wymeersch et al. (2012) "A Machine Learning Approach to
 *            Ranging Error Mitigation for UWB Localization"
 *
 * Knowledge Coverage: L5 Algorithms (NLOS Detection, Mitigation)
 *                      L8 Advanced (Machine Learning NLOS Classification)
 */

#ifndef UWB_NLOS_H
#define UWB_NLOS_H

#include "uwb_types.h"

/* CIR feature vector dimension */
#define UWB_NLOS_FEATURE_DIM 14

/* NLOS classifier types */
typedef enum {
    NLOS_CLASSIFIER_JOINT_LRT    = 0,
    NLOS_CLASSIFIER_SVM         = 1,
    NLOS_CLASSIFIER_DECISION_TREE = 2,
    NLOS_CLASSIFIER_NEURAL_NET  = 3,
    NLOS_CLASSIFIER_LOGISTIC    = 4
} nlos_classifier_type_t;

/* CIR statistical features for NLOS detection (L2: Core Concept) */
typedef struct {
    double kurtosis;
    double skewness;
    double rms_delay_spread_ns;
    double mean_excess_delay_ns;
    double max_excess_delay_ns;
    double total_energy;
    double first_path_energy_ratio;
    double peak_to_energy_ratio;
    double fp_to_total_energy_ratio;
    double cir_growth_rate;
    double rise_time_ps;
    double fp_amplitude_decay;
    double peak_before_fp_db;
    double noise_std_energy;
    double path_loss_db;
} nlos_features_t;

/* NLOS detection result */
typedef struct {
    int is_nlos;
    double nlos_probability;
    double confidence;
    nlos_classifier_type_t classifier_used;
    nlos_features_t features;
} nlos_detection_result_t;

/* NLOS mitigation strategy */
typedef enum {
    NLOS_MITIGATE_REJECT    = 0,
    NLOS_MITIGATE_WEIGHT    = 1,
    NLOS_MITIGATE_CORRECT   = 2,
    NLOS_MITIGATE_SMOOTH    = 3
} nlos_mitigation_t;

/* NLOS mitigation configuration */
typedef struct {
    nlos_mitigation_t strategy;
    double nlos_prob_threshold;
    double nlos_range_bias_m;
    double nlos_variance_scale;
    double smoothing_factor;
} nlos_mitigation_config_t;

/* Feature extraction from CIR */
void nlos_extract_features(const uwb_cir_t *cir, nlos_features_t *features);

/*
 * Joint likelihood ratio test for LOS/NLOS classification.
 * Compares P(CIR | NLOS) / P(CIR | LOS) against a threshold.
 * Reference: Marano et al. (2010) IEEE JSAC
 */
void nlos_lrt_detect(const nlos_features_t *features,
                     nlos_detection_result_t *result);

/*
 * Decision-tree NLOS classifier using CIR features.
 * Uses kurtosis, skewness, and rms delay spread as primary splits.
 */
void nlos_decision_tree_detect(const nlos_features_t *features,
                               nlos_detection_result_t *result);

/*
 * Logistic regression NLOS classifier.
 * P(NLOS) = 1 / (1 + exp(-(w^T * f + b)))
 * with pre-trained weights w and bias b.
 */
void nlos_logistic_detect(const nlos_features_t *features,
                          const double *weights, double bias,
                          nlos_detection_result_t *result);

/*
 * Apply NLOS mitigation to a range measurement.
 * Strategies: reject (set quality=REJECT), weight (inflate variance),
 *             correct (subtract estimated bias), smooth (apply temporal filter)
 */
void nlos_mitigate_range(uwb_ranging_meas_t *meas,
                         const nlos_detection_result_t *nlos_result,
                         const nlos_mitigation_config_t *config);

/*
 * Estimate NLOS range bias from CIR features.
 * Uses a polynomial model: bias = a0 + a1*kurtosis + a2*skewness + ...
 * Reference: Al-Jazzar & Caffery (2002)
 */
double nlos_estimate_range_bias(const nlos_features_t *features);

/*
 * Compute RMS delay spread from CIR.
 * tau_rms = sqrt(sum(P_i * tau_i^2) / sum(P_i) - (mean_tau)^2)
 */
double cir_compute_rms_delay_spread(const uwb_cir_t *cir);

/*
 * Compute kurtosis of CIR envelope.
 * K = E[(|h| - mu)^4] / sigma^4
 * K=3 for Rayleigh (rich multipath), K>3 for sparse/impulsive (NLOS indicator)
 */
double cir_compute_kurtosis(const uwb_cir_t *cir);

/*
 * Compute skewness of CIR envelope.
 * S = E[(|h| - mu)^3] / sigma^3
 * Positive skewness indicates leading-edge NLOS obstruction.
 */
double cir_compute_skewness(const uwb_cir_t *cir);

/*
 * Detect leading edge using energy-based threshold crossing.
 * Searches CIR for first sample exceeding noise_floor * threshold_factor.
 */
int cir_detect_leading_edge(const uwb_cir_t *cir, double threshold_factor,
                            uint16_t *leading_edge_index);

#endif /* UWB_NLOS_H */
