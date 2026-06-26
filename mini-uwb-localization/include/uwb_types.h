/**
 * mini-uwb-localization: UWB Localization Core Types
 *
 * Defines fundamental data structures for Ultra-Wideband (UWB) precision
 * indoor localization systems operating per IEEE 802.15.4z.
 *
 * Reference: Tsui (2005) "Fundamentals of Global Positioning System Receivers"
 * Reference: Sahinoglu, Gezici, Guvenc (2008) "Ultra-wideband Positioning Systems"
 * Reference: IEEE Std 802.15.4z-2020
 *
 * Knowledge Coverage: L1 Definitions
 */

#ifndef UWB_TYPES_H
#define UWB_TYPES_H

#include <stdint.h>
#include <stddef.h>
#include <float.h>

/*
 * IEEE 802.15.4z UWB PHY Parameters
 * L1: Core physical layer definitions
 */

typedef enum {
    UWB_CHANNEL_1  = 1,
    UWB_CHANNEL_2  = 2,
    UWB_CHANNEL_3  = 3,
    UWB_CHANNEL_4  = 4,
    UWB_CHANNEL_5  = 5,
    UWB_CHANNEL_7  = 7,
    UWB_CHANNEL_9  = 9,
    UWB_CHANNEL_11 = 11,
    UWB_CHANNEL_MAX = 11
} uwb_channel_t;

typedef enum {
    UWB_PRF_16MHZ  = 0,
    UWB_PRF_64MHZ  = 1
} uwb_prf_t;

typedef enum {
    UWB_PREAMBLE_64   = 64,
    UWB_PREAMBLE_128  = 128,
    UWB_PREAMBLE_256  = 256,
    UWB_PREAMBLE_512  = 512,
    UWB_PREAMBLE_1024 = 1024,
    UWB_PREAMBLE_1536 = 1536,
    UWB_PREAMBLE_2048 = 2048,
    UWB_PREAMBLE_4096 = 4096
} uwb_preamble_len_t;

typedef enum {
    UWB_DATARATE_110K  = 0,
    UWB_DATARATE_850K  = 1,
    UWB_DATARATE_6M8   = 2,
    UWB_DATARATE_27M   = 3,
    UWB_DATARATE_31M2  = 4
} uwb_datarate_t;

/* Speed of light in m/s (exact, defined by SI) */
#define UWB_C            299792458.0
#define UWB_BW_DEFAULT   499200000.0
#define UWB_BW_EFFECTIVE 450000000.0
#define UWB_COARSE_RESOLUTION_NS  2.003
#define UWB_FINE_RESOLUTION_NS    (1.0 / UWB_BW_EFFECTIVE * 1e9)
#define UWB_DISTANCE_RESOLUTION_M (UWB_C / (2.0 * UWB_BW_EFFECTIVE))
#define UWB_MAX_RANGE_16MHZ_M     (UWB_C / (2.0 * 16e6))
#define UWB_MAX_RANGE_64MHZ_M     (UWB_C / (2.0 * 64e6))

typedef struct {
    double x;
    double y;
} uwb_pos2d_t;

typedef struct {
    double x;
    double y;
    double z;
} uwb_pos3d_t;

typedef struct {
    double var_x;
    double var_y;
    double var_z;
    double cov_xy;
    double cov_xz;
    double cov_yz;
} uwb_covariance_t;

/* L1: Anchor (fixed reference) node in UWB localization system */
typedef struct {
    uint16_t id;
    uwb_pos3d_t position;
    uwb_channel_t channel;
    uwb_prf_t prf;
    double tx_power_dbm;
    double antenna_gain_dbi;
    double clock_offset_ppm;
    int is_active;
    int is_synchronized;
    uint64_t last_ranging_ts;
} uwb_anchor_t;

/* L1: Tag (mobile) node in UWB localization system */
typedef struct {
    uint16_t id;
    uwb_pos3d_t position;
    uwb_pos3d_t velocity;
    uwb_covariance_t pos_cov;
    double clock_offset_ppm;
    double battery_voltage;
    int motion_state;
    uint64_t last_update_ts;
} uwb_tag_t;

/* L1: Ranging protocol type */
typedef enum {
    UWB_RANGING_SS_TWR  = 0,
    UWB_RANGING_DS_TWR  = 1,
    UWB_RANGING_TDOA    = 2,
    UWB_RANGING_TOA     = 3,
    UWB_RANGING_PDOA    = 4
} uwb_ranging_type_t;

/* L1: Quality indicator for a ranging measurement */
typedef enum {
    UWB_RANGE_QUALITY_EXCELLENT = 0,
    UWB_RANGE_QUALITY_GOOD      = 1,
    UWB_RANGE_QUALITY_FAIR      = 2,
    UWB_RANGE_QUALITY_POOR      = 3,
    UWB_RANGE_QUALITY_REJECT    = 4
} uwb_range_quality_t;

/* L1: Single UWB ranging measurement between tag and anchor */
typedef struct {
    uint16_t anchor_id;
    uint16_t tag_id;
    uwb_ranging_type_t type;
    double distance_m;
    double distance_variance;
    double tof_seconds;
    double rssi_dbm;
    double fp_power_dbm;
    double fp_amplitude;
    double noise_stddev;
    double cir_growth_rate;
    double cir_peak_to_leading;
    uwb_range_quality_t quality;
    uint16_t sequence_number;
    uint64_t timestamp;
} uwb_ranging_meas_t;

/* L1: Channel Impulse Response - raw CIR data from UWB receiver */
#define UWB_CIR_MAX_SAMPLES 1016

typedef struct {
    int32_t i;
    int32_t q;
} uwb_cir_sample_t;

typedef struct {
    uwb_cir_sample_t samples[UWB_CIR_MAX_SAMPLES];
    uint16_t num_samples;
    uint16_t first_path_index;
    uint16_t peak_path_index;
    double sampling_period_ps;
    double noise_floor;
    double cir_power;
    double preamble_accumulation_count;
} uwb_cir_t;

/* L1: UWB system configuration */
typedef struct {
    uwb_channel_t channel;
    uwb_prf_t prf;
    uwb_preamble_len_t preamble_len;
    uwb_datarate_t datarate;
    uint8_t preamble_code;
    uint8_t sfd_id;
    double tx_power_dbm;
    double antenna_delay_tx_ps;
    double antenna_delay_rx_ps;
    int smart_power_enabled;
    int leading_edge_detection;
    uint32_t ranging_interval_ms;
    uint16_t slot_duration_us;
    uint16_t response_timeout_us;
} uwb_config_t;

/* L1: Localization error statistics for performance evaluation */
typedef struct {
    double mean_error_2d;
    double mean_error_3d;
    double rmse_2d;
    double rmse_3d;
    double cep50;
    double cep90;
    double sep50;
    double sep90;
    double max_error;
    double min_error;
    double stddev_error;
    uint32_t num_samples;
    uint32_t num_outliers;
    double crlb;
    double gdop_avg;
} uwb_error_metrics_t;

/* Initialization and Utility Functions */
void uwb_anchor_init(uwb_anchor_t *anchor, uint16_t id, double x, double y, double z);
void uwb_tag_init(uwb_tag_t *tag, uint16_t id);
void uwb_config_init_default(uwb_config_t *config);
void uwb_ranging_meas_init(uwb_ranging_meas_t *meas);
void uwb_error_metrics_init(uwb_error_metrics_t *metrics);

double uwb_distance_2d(const uwb_pos2d_t *a, const uwb_pos2d_t *b);
double uwb_distance_3d(const uwb_pos3d_t *a, const uwb_pos3d_t *b);
double uwb_channel_frequency_hz(uwb_channel_t channel);
double uwb_channel_bandwidth_hz(uwb_channel_t channel);
double uwb_crlb_distance(double snr_linear, double bw_effective);
double uwb_compute_gdop(const uwb_pos3d_t *anchors, int num_anchors,
                        const uwb_pos3d_t *tag_pos);

static inline double uwb_tof_to_distance(double tof_seconds) {
    return UWB_C * tof_seconds / 2.0;
}

static inline double uwb_distance_to_tof(double distance_m) {
    return 2.0 * distance_m / UWB_C;
}

static inline double uwb_dbm_to_mw(double dbm) {
    double x = dbm / 10.0, xln10 = x * 2.302585092994046;
    double result = 1.0, term = 1.0;
    for (int i = 1; i < 20; i++) {
        term *= xln10 / i;
        result += term;
        if (term < 1e-15) break;
    }
    return result;
}

static inline double uwb_mw_to_dbm(double mw) {
    if (mw <= 0.0) return -200.0;
    double x = mw;
    int exp = 0;
    while (x >= 10.0) { x /= 10.0; exp++; }
    while (x < 1.0)   { x *= 10.0; exp--; }
    double y = (x - 1.0) / (x + 1.0), y2 = y * y, ln = 0.0, yp = y;
    for (int i = 1; i < 30; i += 2) {
        ln += yp / i;
        yp *= y2;
        if (yp / i < 1e-15) break;
    }
    ln *= 2.0;
    return 10.0 * (ln / 2.302585092994046 + exp);
}

static inline double uwb_compute_snr_db(double fp_amplitude, double noise_stddev) {
    if (noise_stddev <= 0.0) return 100.0;
    double snr_linear = (fp_amplitude * fp_amplitude) / (noise_stddev * noise_stddev);
    return uwb_mw_to_dbm(snr_linear);
}

#endif /* UWB_TYPES_H */
