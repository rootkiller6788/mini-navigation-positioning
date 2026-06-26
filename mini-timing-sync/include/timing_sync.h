/**
 * @file timing_sync.h
 * @brief Core timing synchronization data structures and API
 *
 * Knowledge Coverage:
 *   L1 - Core definitions: TimeOffset, ClockSkew, ClockDrift, ClockAging,
 *        SyncState, HoldoverState, Timestamp
 *   L2 - Core concepts: Two-way time transfer, clock discipline, holdover
 *   L3 - Mathematical structures: Clock state vector, Kalman state
 *   L4 - Fundamental laws: Two-way time transfer equation, offset-delay relation
 *
 * Reference: IEEE 1588-2019, ITU-T G.8271, IETF RFC 5905 (NTPv4)
 * Textbook: Mills, D.L. "Computer Network Time Synchronization" (2011)
 *           Tsui, J. "Fundamentals of GPS Receivers" (2005) Ch.10
 *
 * Course Mapping:
 *   Stanford EE359 - Wireless Communications (timing recovery)
 *   MIT 6.450 - Digital Communications (symbol timing sync)
 *   Berkeley EE123 - DSP (clock recovery)
 *   CMU 18-345 - Intro to Telecom Networks (clock sync)
 *   ETH 227-0436 - Communications (synchronization)
 */

#ifndef TIMING_SYNC_H
#define TIMING_SYNC_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * L1: Core Definitions - Data Structures
 * ============================================================================
 */

/** Timestamp in nanoseconds since epoch (PTP/TAI epoch: 1970-01-01) */
typedef struct {
    int64_t seconds;      /**< Integer seconds since epoch */
    int32_t nanoseconds;  /**< Nanosecond fraction [0, 999999999] */
} Timestamp;

/** Time offset: difference between two clocks (slave - master) */
typedef struct {
    double offset_ns;     /**< Time offset in nanoseconds */
    double variance_ns2;  /**< Variance estimate of offset measurement */
} TimeOffset;

/** Clock skew: first derivative of time offset (frequency offset) */
typedef struct {
    double skew_ppb;      /**< Frequency offset in parts-per-billion */
    double variance_ppb2; /**< Variance estimate of skew estimate */
} ClockSkew;

/** Clock drift: second derivative of time offset (aging rate) */
typedef struct {
    double drift_ppb_per_day;  /**< Aging rate in ppb per day */
} ClockDrift;

/** Full clock state vector for Kalman tracking */
typedef struct {
    double offset_ns;         /**< Time offset state [ns] */
    double freq_offset_ppb;   /**< Frequency offset state [ppb] */
    double drift_ppb_per_day; /**< Drift/aging state [ppb/day] */
    double P[3][3];           /**< 3x3 state covariance matrix */
    int64_t last_update_time; /**< Time of last state update */
} ClockState;

/** PTP delay mechanism types (IEEE 1588 Section 8.3) */
typedef enum {
    DELAY_MECH_E2E = 0,    /**< End-to-End delay mechanism */
    DELAY_MECH_P2P = 1,    /**< Peer-to-Peer delay mechanism */
    DELAY_MECH_DISABLED = 2
} DelayMechanism;

/** PTP port state (IEEE 1588 Section 9.2.6.11) */
typedef enum {
    PTP_INITIALIZING  = 0,
    PTP_FAULTY        = 1,
    PTP_DISABLED      = 2,
    PTP_LISTENING     = 3,
    PTP_PRE_MASTER    = 4,
    PTP_MASTER        = 5,
    PTP_PASSIVE       = 6,
    PTP_UNCALIBRATED  = 7,
    PTP_SLAVE         = 8
} PtpPortState;

/** Clock synchronization status */
typedef enum {
    SYNC_FREE_RUNNING = 0,   /**< No external reference */
    SYNC_ACQUIRING    = 1,   /**< Locking to reference */
    SYNC_LOCKED       = 2,   /**< Phase-locked to reference */
    SYNC_HOLDOVER     = 3,   /**< Reference lost, using internal model */
    SYNC_LOS          = 4    /**< Loss of sync, no valid time */
} SyncStatus;

/** Holdover specification parameters */
typedef struct {
    double entry_threshold_ns;  /**< Max offset to enter holdover [ns] */
    double exit_threshold_ns;   /**< Max offset to exit holdover [ns] */
    double max_holdover_duration_s; /**< Maximum holdover duration [s] */
    double holdover_accuracy_ns;    /**< Expected holdover accuracy [ns] */
} HoldoverConfig;

/** NTP stratum levels (RFC 5905 Section 3) */
typedef enum {
    NTP_STRATUM_UNSPECIFIED = 0,
    NTP_STRATUM_PRIMARY     = 1,   /**< Directly connected to reference clock */
    NTP_STRATUM_SECONDARY   = 2,   /**< Synced to stratum 1 via NTP */
    NTP_STRATUM_MAX         = 15,  /**< Maximum valid stratum */
    NTP_STRATUM_UNSYNC      = 16   /**< Unsynchronized */
} NtpStratum;

/** PTP clock quality (IEEE 1588 Section 5.3.7) */
typedef struct {
    uint8_t clock_class;
    uint8_t clock_accuracy;
    uint16_t offset_scaled_log_variance;
} ClockQuality;

/** PTP timestamp record for two-way exchange */
typedef struct {
    Timestamp t1;  /**< Sync departure time at master */
    Timestamp t2;  /**< Sync arrival time at slave */
    Timestamp t3;  /**< Delay_req departure time at slave */
    Timestamp t4;  /**< Delay_req arrival time at master */
} PtpTimestamps;

/** NTP packet header (RFC 5905 Section 7.3) */
typedef struct {
    uint8_t  li_vn_mode;
    uint8_t  stratum;
    int8_t   poll_interval;
    int8_t   precision;
    int32_t  root_delay;
    int32_t  root_dispersion;
    uint32_t ref_id;
    Timestamp ref_timestamp;
    Timestamp orig_timestamp;
    Timestamp recv_timestamp;
    Timestamp xmit_timestamp;
} NtpPacket;

/* Core API */

int timing_compute_offset_delay(const PtpTimestamps *ts,
                                double *offset, double *delay);

int timing_ntp_offset_delay(const Timestamp *T1, const Timestamp *T2,
                            const Timestamp *T3, const Timestamp *T4,
                            double *offset, double *delay);

double timing_timestamp_diff_ns(const Timestamp *ts1, const Timestamp *ts2);
void timing_timestamp_add_ns(Timestamp *ts, int64_t ns_to_add);
int timing_timestamp_cmp(const Timestamp *a, const Timestamp *b);
int timing_timestamp_valid(const Timestamp *ts);
double timing_timestamp_to_double(const Timestamp *ts);
Timestamp timing_double_to_timestamp(double ns);

/* PI Servo */

typedef struct {
    double Kp;
    double Ki;
    double integral_limit_ns;
    double max_correction_ns;
} PiServoConfig;

typedef struct {
    double integral_ns;
    int first_update;
} PiServoState;

void pi_servo_init(PiServoConfig *cfg, PiServoState *state);
double pi_servo_update(PiServoConfig *cfg, PiServoState *state,
                       double measured_offset_ns);
void pi_servo_reset(PiServoState *state);

/* Kalman Filter Clock Tracking */

void clock_state_init(ClockState *state, double initial_offset_ns,
                      double initial_freq_ppb, double initial_drift_ppb_per_day,
                      double offset_uncertainty_ns,
                      double freq_uncertainty_ppb);
void clock_kalman_predict(ClockState *state, double dt, const double Q[9]);
void clock_kalman_update(ClockState *state, double measured_offset_ns,
                         double measurement_noise_var_ns2);
double clock_predict_offset(const ClockState *state, double dt_s);
double clock_get_freq_ppb(const ClockState *state);

/* Holdover */

int holdover_should_enter(double offset_ns, const HoldoverConfig *cfg);
int holdover_can_exit(double offset_ns, const HoldoverConfig *cfg,
                      double holdover_duration_s);
double holdover_estimate_uncertainty(double holdover_duration_s,
                                     double drift_ppb_per_day,
                                     double aging_ppb_per_day2,
                                     double initial_error_ns);

/* GPS Timing */

typedef struct {
    double antenna_cable_delay_ns;
    double receiver_processing_delay_ns;
    int    sawtooth_correction_enable;
} GpsTimingConfig;

int gps_1pps_compute_offset(const Timestamp *local_timestamp,
                            const GpsTimingConfig *cfg, double *offset);
double gps_sawtooth_correction(double raw_offset_ns,
                               double sawtooth_amplitude_ns,
                               double phase_radians);

#ifdef __cplusplus
}
#endif

#endif /* TIMING_SYNC_H */
