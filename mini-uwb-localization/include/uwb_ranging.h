/**
 * mini-uwb-localization: UWB Ranging Algorithms
 *
 * Implements Two-Way Ranging (TWR) protocols for precision distance
 * measurement between UWB nodes. Includes Single-Sided TWR (SS-TWR),
 * Double-Sided TWR (DS-TWR), and clock drift compensation.
 *
 * Reference: IEEE 802.15.4z-2020, Section 15.4 (UWB PHY ranging)
 * Reference: Decawave DW1000 User Manual (APS013: TWR)
 *
 * Knowledge Coverage: L2 Core Concepts, L4 Error Budgets, L5 Algorithms
 */

#ifndef UWB_RANGING_H
#define UWB_RANGING_H

#include "uwb_types.h"

#define UWB_TWR_DEFAULT_REPLY_DELAY_US  5000
#define UWB_TWR_MIN_REPLY_DELAY_US      100
#define UWB_TWR_MAX_REPLY_DELAY_US      50000
#define UWB_TWR_CLOCK_TOLERANCE_PPM     20.0
#define UWB_TWR_DEFAULT_BURST_SIZE      8
#define UWB_TWR_TRIMMED_MEAN_RATIO      0.25

typedef enum {
    TWR_STATE_IDLE        = 0,
    TWR_STATE_POLL_SENT   = 1,
    TWR_STATE_RESP_WAIT   = 2,
    TWR_STATE_FINAL_SENT  = 3,
    TWR_STATE_COMPLETE    = 4,
    TWR_STATE_TIMEOUT     = 5,
    TWR_STATE_ERROR       = 6
} twr_state_t;

typedef struct {
    uint64_t t_poll_tx;
    uint64_t t_poll_rx;
    uint64_t t_resp_tx;
    uint64_t t_resp_rx;
    double clock_offset_ppm;
} twr_ss_timestamps_t;

typedef struct {
    uint64_t t_poll_tx;
    uint64_t t_poll_rx;
    uint64_t t_resp_tx;
    uint64_t t_resp_rx;
    uint64_t t_final_tx;
    uint64_t t_final_rx;
    double clock_offset_ppm;
} twr_ds_timestamps_t;

typedef struct {
    double distance_m;
    double tof_seconds;
    double distance_stddev_m;
    double clock_error_est_m;
    double rssi_dbm;
    twr_state_t final_state;
    uint32_t round_trip_time_us;
    uint8_t num_retries;
} twr_result_t;

typedef struct {
    twr_state_t state;
    uint16_t session_id;
    uint16_t initiator_id;
    uint16_t responder_id;
    twr_ss_timestamps_t ss_ts;
    twr_ds_timestamps_t ds_ts;
    twr_result_t last_result;
    double antenna_delay_initiator_ps;
    double antenna_delay_responder_ps;
    uint32_t reply_delay_us;
    uint32_t timeout_us;
    uint8_t burst_size;
} twr_session_t;

typedef struct {
    double distances[UWB_TWR_DEFAULT_BURST_SIZE];
    int count;
    double mean;
    double stddev;
    double trimmed_mean;
    double median;
    double min_val;
    double max_val;
} ranging_burst_stats_t;

typedef struct {
    double offset_ppm;
    double drift_rate_ppb;
    double allan_variance;
    double temperature_coeff_ppm_per_c;
    uint64_t last_calib_ts;
} clock_drift_model_t;

/* TWR Protocol Functions */
void twr_session_init(twr_session_t *session, uint16_t initiator_id,
                      uint16_t responder_id);

/*
 * SS-TWR ToF:
 * T_tof = (T_round - T_reply) / 2
 * where T_round = t_resp_rx - t_poll_tx, T_reply = t_resp_tx - t_poll_rx
 * Reference: Decawave APS013, Section 2.1
 */
double twr_ss_compute_tof(const twr_ss_timestamps_t *ts, double tick_period_ps);

/*
 * SS-TWR distance with antenna delay compensation:
 * d = c * (T_tof - T_ant_init - T_ant_resp) / 2
 */
double twr_ss_compute_distance(const twr_ss_timestamps_t *ts,
                               double tick_period_ps,
                               double ant_delay_init_ps,
                               double ant_delay_resp_ps);

/*
 * DS-TWR ToF (symmetric double-sided):
 * T_tof = (T_round1 * T_round2 - T_reply1 * T_reply2) /
 *         (T_round1 + T_round2 + T_reply1 + T_reply2)
 * Eliminates first-order clock drift error.
 * Reference: Neirynck et al. (2016)
 */
double twr_ds_compute_tof(const twr_ds_timestamps_t *ts, double tick_period_ps);

double twr_ds_compute_distance(const twr_ds_timestamps_t *ts,
                               double tick_period_ps,
                               double ant_delay_init_ps,
                               double ant_delay_resp_ps);

/*
 * Estimate clock frequency offset from DS-TWR:
 * k = f_initiator / f_responder
 * k ~= (T_round1 + T_reply2) / (T_reply1 + T_round2)
 * clock_offset_ppm = (k - 1) * 1e6
 */
double twr_ds_estimate_clock_offset(twr_ds_timestamps_t *ts, double tick_period_ps);

double twr_compensate_clock_offset(double distance_m, double clock_offset_ppm);
double twr_estimate_clock_error_m(double reply_time_us, double clock_offset_ppm);

/* Ranging Burst Processing */
void ranging_burst_compute_stats(const double *distances, int count,
                                 ranging_burst_stats_t *stats);
int ranging_burst_median_filter(double *distances, int count, int window_size);

/*
 * Skewness-and-kurtosis NLOS detection on a ranging burst.
 * NLOS indicators: positive skewness, high kurtosis (>3).
 * Reference: Guvenc et al. (2007) "NLOS Identification for UWB"
 */
double ranging_burst_nlos_score(const double *distances, int count);
double ranging_ewma_filter(const double *history, int count, double alpha);

/*
 * Total ranging error budget (RSS of contributors):
 * sigma_total = sqrt(sigma_noise^2 + sigma_clock^2 + sigma_mp^2 + sigma_ant^2)
 * sigma_noise = c / (2*pi*B_eff*sqrt(2*SNR))
 * sigma_clock = T_reply * c * ppm_error / 2
 */
double twr_error_budget(double snr_linear, double bw_effective,
                        double reply_time_us, double clock_ppm_error,
                        double multipath_error_m, double ant_cal_error_m);

#endif /* UWB_RANGING_H */
