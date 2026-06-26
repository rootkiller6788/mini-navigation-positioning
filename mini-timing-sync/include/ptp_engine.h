/**
 * @file ptp_engine.h
 * @brief Precision Time Protocol (IEEE 1588-2019) engine
 *
 * L1: PTP message types, port states, clock identity
 * L2: BMCA (Best Master Clock Algorithm), delay measurement
 * L3: PTP clock servo control theory
 * L4: IEEE 1588 two-way time transfer equations
 * L5: BMCA implementation, PTP message processing
 * L6: Master-slave synchronization loop
 * L7: 5G fronthaul timing (IEEE 802.1CM profile), power grid sync
 */

#ifndef PTP_ENGINE_H
#define PTP_ENGINE_H
#include <stdint.h>
#include "timing_sync.h"
#ifdef __cplusplus
extern "C" {
#endif

/* L1: PTP clock identity (EUI-64) */
typedef struct {
    uint8_t octets[8];
} PtpClockIdentity;

/* L1: PTP message types (IEEE 1588 Table 21) */
typedef enum {
    PTP_MSG_SYNC              = 0,
    PTP_MSG_DELAY_REQ         = 1,
    PTP_MSG_PDELAY_REQ        = 2,
    PTP_MSG_PDELAY_RESP       = 3,
    PTP_MSG_FOLLOW_UP         = 8,
    PTP_MSG_DELAY_RESP        = 9,
    PTP_MSG_PDELAY_RESP_FU    = 10,
    PTP_MSG_ANNOUNCE          = 11,
    PTP_MSG_SIGNALING         = 12,
    PTP_MSG_MANAGEMENT        = 13
} PtpMessageType;

/* L1: PTP common message header (IEEE 1588 Section 13.3) */
typedef struct {
    uint8_t  transport_specific;
    uint8_t  message_type;
    uint8_t  reserved0;
    uint8_t  version_ptp;
    uint16_t message_length;
    uint8_t  domain_number;
    uint8_t  reserved1;
    uint16_t flags;
    int64_t  correction_field;
    uint32_t reserved2;
    PtpClockIdentity source_port_identity;
    uint16_t sequence_id;
    uint8_t  control_field;
    int8_t   log_message_interval;
} PtpHeader;

/* L1: PTP announce message body */
typedef struct {
    Timestamp      origin_timestamp;
    uint16_t       current_utc_offset;
    uint8_t        reserved;
    uint8_t        grandmaster_priority1;
    ClockQuality   grandmaster_clock_quality;
    uint8_t        grandmaster_priority2;
    PtpClockIdentity grandmaster_identity;
    uint16_t       steps_removed;
    uint8_t        time_source;
} PtpAnnounceBody;

/* L1: PTP sync/follow-up message body */
typedef struct {
    Timestamp origin_timestamp;
} PtpSyncBody;

/* L2: BMCA dataset comparison result */
typedef enum {
    BMCA_THIS_IS_BETTER    = -1,
    BMCA_EQUAL             = 0,
    BMCA_OTHER_IS_BETTER   = 1
} BmcaResult;

/* L2: BMCA input dataset for comparison (IEEE 1588 Section 10.3) */
typedef struct {
    uint8_t        priority1;
    ClockQuality   clock_quality;
    uint8_t        priority2;
    PtpClockIdentity clock_identity;
    uint16_t       steps_removed;
} BmcaDataset;

/* L2: PTP slave clock tracking state */
typedef struct {
    PtpTimestamps   last_timestamps;  /* T1, T2, T3, T4 */
    double          offset_from_master_ns;
    double          mean_path_delay_ns;
    double          offset_std_dev_ns;
    double          delay_std_dev_ns;
    PiServoConfig   servo_config;
    PiServoState    servo_state;
    ClockState      kalman_clock;
    SyncStatus      sync_status;
    int             sync_count;
    int             valid_sync_count;
} PtpSlaveState;

/* L2: PTP master clock state */
typedef struct {
    PtpClockIdentity clock_id;
    ClockQuality     clock_quality;
    uint8_t          priority1;
    uint8_t          priority2;
    uint8_t          domain_number;
    double           local_time_offset_ns;
    uint16_t         current_utc_offset;
    uint8_t          time_source;
} PtpMasterState;

/**
 * L4: Process PTP two-way exchange to get offset and delay.
 * Equation: offset = ((t2-t1) - (t4-t3)) / 2
 *           delay  = ((t2-t1) + (t4-t3)) / 2
 */
int ptp_process_timestamps(PtpSlaveState *state);

/**
 * L5: BMCA dataset comparison algorithm (IEEE 1588 Section 10.3.8).
 *
 * Comparison order:
 * 1. priority1 (lower is better)
 * 2. clock_class (lower is better)
 * 3. clock_accuracy (lower is better)
 * 4. offset_scaled_log_variance (lower is better)
 * 5. priority2 (lower is better)
 * 6. clock_identity (lower is better, tiebreaker)
 */
BmcaResult bmca_compare_datasets(const BmcaDataset *A, const BmcaDataset *B);

/**
 * L5: Compute PTP announce interval from log_message_interval.
 * interval = 2^log_message_interval seconds
 */
double ptp_announce_interval(int8_t log_message_interval);

/**
 * L5: Initialize PTP slave state for synchronization.
 *
 * @param state            State to initialize
 * @param initial_offset_ns Estimated initial offset
 * @param update_interval_s Servo update interval [seconds]
 */
void ptp_slave_init(PtpSlaveState *state, double initial_offset_ns,
                    double update_interval_s);

/**
 * L5: Update PTP slave with new timestamps and compute correction.
 *
 * Full processing chain:
 * 1. Validate timestamps (t1 < t2, t3 < t4, logical consistency)
 * 2. Compute offset and path delay
 * 3. Update Kalman filter with measured offset
 * 4. Run PI servo to compute clock correction
 *
 * @param state  Slave state (updated in place)
 * @param ts     New PTP timestamps
 * @return Clock correction value in nanoseconds
 */
double ptp_slave_update(PtpSlaveState *state, const PtpTimestamps *ts);

/**
 * L6: Compute PTP sync accuracy estimate.
 *
 * Uncertainty budget:
 *   U_total^2 = U_timestamp^2 + U_delay_asymmetry^2 + U_oscillator^2
 *
 * @param state  Current slave state
 * @return Estimated sync accuracy in nanoseconds (1-sigma)
 */
double ptp_sync_accuracy(const PtpSlaveState *state);

/**
 * L6: Check if PTP slave is within specified accuracy target.
 *
 * @param state           Current slave state
 * @param accuracy_ns     Required accuracy [ns]
 * @param confidence_k    Confidence factor (k=2 for 95%, k=3 for 99.7%)
 * @return 1 if within spec, 0 otherwise
 */
int ptp_within_accuracy(const PtpSlaveState *state, double accuracy_ns,
                        double confidence_k);

/**
 * L5: Validate PTP timestamp consistency.
 *
 * Checks:
 * - t1, t2, t3, t4 have valid format
 * - t1 < t2 (causality: sync arrives after departure)
 * - t3 < t4 (causality: delay_req arrives after departure)
 * - t1 < t3 (sync before delay_req) for normal sequence
 *
 * @param ts  Timestamps to validate
 * @return 1 if valid, 0 if inconsistent
 */
int ptp_validate_timestamps(const PtpTimestamps *ts);

/**
 * L7: 5G fronthaul timing profile check (ITU-T G.8271.1).
 *
 * Class A: |offset| < 1100 ns (peak-to-peak)
 * Class B: |offset| < 260 ns (peak-to-peak)
 * Class C: |offset| < 30 ns (peak-to-peak)
 *
 * @param offset_ns  Measured offset from master
 * @param class_level 'A', 'B', or 'C'
 * @return 1 if within spec, 0 if not
 */
int ptp_5g_fronthaul_check(double offset_ns, char class_level);

/**
 * L7: IEC 61850 power grid timing check (sampled values).
 *
 * IEC 61850-9-2 requires timing accuracy < 1 us for
 * sampled value messages in digital substations.
 *
 * @param offset_ns  Measured offset from master
 * @return 1 if within spec (< 1000 ns), 0 otherwise
 */
int ptp_power_grid_check(double offset_ns);

/**
 * L6: Compute next sync transmission time based on log interval.
 *
 * @param log_sync_interval  Log2 of sync interval in seconds
 * @param current_time       Current master time
 * @param next_time          [out] Computed next sync time
 */
void ptp_next_sync_time(int8_t log_sync_interval,
                        const Timestamp *current_time,
                        Timestamp *next_time);

#ifdef __cplusplus
}
#endif
#endif
