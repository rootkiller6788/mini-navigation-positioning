/**
 * mini-uwb-localization: UWB Channel Models
 *
 * Implements IEEE 802.15.4a UWB channel models for indoor localization
 * simulation: path loss, multipath profiles (LOS/NLOS residential,
 * office, industrial), Saleh-Valenzuela model, and small-scale fading.
 *
 * Reference: Molisch et al. (2004) "IEEE 802.15.4a Channel Model"
 * Reference: Foerster et al. (2003) "Channel Modeling Sub-committee
 *            Report Final" IEEE P802.15-02/490r1-SG3a
 *
 * Knowledge Coverage: L3 Mathematical Structures (SV model, fading statistics)
 *                      L4 Fundamental Laws (Friis for UWB)
 *                      L6 Canonical Problems (Indoor channel modeling)
 */

#ifndef UWB_CHANNEL_H
#define UWB_CHANNEL_H

#include "uwb_types.h"

#define UWB_CHANNEL_MAX_PATHS  200
#define UWB_CHANNEL_MAX_CLUSTERS 10

/* IEEE 802.15.4a channel environment */
typedef enum {
    UWB_ENV_RESIDENTIAL_LOS   = 0,
    UWB_ENV_RESIDENTIAL_NLOS  = 1,
    UWB_ENV_OFFICE_LOS        = 2,
    UWB_ENV_OFFICE_NLOS       = 3,
    UWB_ENV_INDUSTRIAL_LOS    = 4,
    UWB_ENV_INDUSTRIAL_NLOS   = 5,
    UWB_ENV_OPEN_OUTDOOR_LOS  = 6,
    UWB_ENV_OPEN_OUTDOOR_NLOS = 7
} uwb_channel_environment_t;

/* Single multipath component in SV model */
typedef struct {
    double amplitude;
    double delay_ns;
    double phase_rad;
    int cluster_index;
} uwb_mpc_t;

/* Saleh-Valenzuela channel impulse response model */
typedef struct {
    uwb_mpc_t paths[UWB_CHANNEL_MAX_PATHS];
    int num_paths;
    int num_clusters;
    double cluster_decay_ns;
    double ray_decay_ns;
    double cluster_arrival_rate_1_per_ns;
    double ray_arrival_rate_1_per_ns;
    double path_loss_db;
    double shadowing_db;
    double distance_m;
    uwb_channel_environment_t environment;
    double center_frequency_hz;
    double bandwidth_hz;
} uwb_sv_channel_t;

/* Path loss model parameters (L4: Friis extension for UWB) */
typedef struct {
    double pl0_db;
    double path_loss_exponent;
    double shadowing_stddev_db;
    double frequency_dependence;
    double breakpoint_distance_m;
} uwb_path_loss_model_t;

/* Channel statistics */
typedef struct {
    double mean_excess_delay_ns;
    double rms_delay_spread_ns;
    double coherence_bandwidth_hz;
    double max_excess_delay_ns;
    int num_significant_paths_10db;
    int num_significant_paths_20db;
    double total_energy;
    double nakagami_m;
} uwb_channel_statistics_t;

/*
 * Initialize SV channel model with environment-specific parameters.
 * IEEE 802.15.4a parameters per environment.
 */
void uwb_channel_sv_init(uwb_sv_channel_t *channel,
                         uwb_channel_environment_t env,
                         double distance_m,
                         double fc_hz, double bw_hz);

/*
 * Generate multipath components using Saleh-Valenzuela model.
 * Clusters and rays have Poisson-distributed arrival times,
 * amplitudes decay exponentially with cluster and ray decay constants.
 */
void uwb_channel_generate_paths(uwb_sv_channel_t *channel, uint32_t seed);

/*
 * Compute path loss for given distance and environment.
 * PL(d) = PL0 + 10*n*log10(d/d0) + X_sigma
 * where X_sigma ~ N(0, sigma_shadowing^2)
 */
double uwb_channel_path_loss_db(double distance_m,
                                const uwb_path_loss_model_t *model);

/*
 * Compute path loss model parameters for a given environment.
 * Table-driven parameters per IEEE 802.15.4a.
 */
void uwb_channel_get_path_loss_model(uwb_channel_environment_t env,
                                     uwb_path_loss_model_t *model);

/*
 * Compute wideband channel statistics from generated paths.
 * RMS delay spread: tau_rms = sqrt(mean(tau^2) - mean(tau)^2)
 * Coherence BW: B_c = 1/(5 * tau_rms) (50% correlation)
 *              or 1/(50 * tau_rms) (90% correlation)
 */
void uwb_channel_compute_statistics(const uwb_sv_channel_t *channel,
                                    uwb_channel_statistics_t *stats);

/*
 * Apply log-normal shadowing to a channel realization.
 * Adds spatially correlated shadowing using Gudmundson model.
 */
void uwb_channel_apply_shadowing(uwb_sv_channel_t *channel,
                                 double shadowing_std_db, uint32_t seed);

/*
 * Compute Nakagami-m fading parameter from environment.
 * m = 1 is Rayleigh (worst), m > 1 is less severe, m -> inf is no fading.
 */
double uwb_channel_nakagami_m(uwb_channel_environment_t env);

/*
 * Generate a discrete-time CIR from the SV model paths
 * for use in receiver simulation.
 */
void uwb_channel_to_cir(const uwb_sv_channel_t *sv_channel,
                        uwb_cir_t *cir,
                        double sampling_period_ps,
                        int num_samples);

/*
 * Compute received signal power using Friis transmission equation
 * extended for UWB bandwidth.
 * P_rx = P_tx + G_tx + G_rx - PL(d)
 */
double uwb_friis_rx_power_dbm(double tx_power_dbm,
                               double tx_gain_dbi,
                               double rx_gain_dbi,
                               double distance_m,
                               const uwb_path_loss_model_t *model);

/*
 * Time-varying channel: apply Doppler effect to CIR.
 * Each path phase rotates by 2*pi*f_d*t where f_d = v*fc/c.
 */
void uwb_channel_apply_doppler(uwb_sv_channel_t *channel,
                               double velocity_ms,
                               double time_s);

/*
 * Compute signal bandwidth needed for given range accuracy.
 * From CRLB: sigma_d = c / (2*pi*sqrt(2*SNR)*B_eff)
 * Rearranged: B_eff = c / (2*pi*sqrt(2*SNR)*sigma_d)
 */
double uwb_bandwidth_for_accuracy(double required_accuracy_m, double snr_linear);

#endif /* UWB_CHANNEL_H */
