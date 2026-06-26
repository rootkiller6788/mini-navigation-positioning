/**
 * @file ptp_engine.c
 * @brief Precision Time Protocol (IEEE 1588) engine implementation
 *
 * Implements: BMCA comparison, PTP timestamp processing,
 * slave clock servo loop, sync accuracy estimation,
 * 5G fronthaul and power grid compliance checks.
 */

#include "ptp_engine.h"
#include <math.h>
#include <string.h>

/* ===================================================================
 * L4: PTP Two-way Exchange Processing
 * =================================================================== */

int ptp_process_timestamps(PtpSlaveState *state)
{
    if (!state) return -1;

    double offset_ns, delay_ns;

    /* Delegate to core two-way computation */
    if (timing_compute_offset_delay(&state->last_timestamps,
                                    &offset_ns, &delay_ns) != 0) {
        return -1;
    }

    /* Update state with new measurements */
    state->offset_from_master_ns = offset_ns;
    state->mean_path_delay_ns = delay_ns;

    /* Running statistics using exponential moving average of
     * absolute deviation from running mean. This approach avoids
     * static variables and is state-object safe.
     *
     * MAD (Mean Absolute Deviation) as robust estimator:
     *   mad[n] = (1-alpha)*mad[n-1] + alpha*|x[n] - mean[n]|
     *   mean[n] = (1-alpha)*mean[n-1] + alpha*x[n]
     */
    double alpha = 0.1;

    /* Use the offset_std_dev fields to store running MAD estimate */
    double old_mad_offset = state->offset_std_dev_ns;
    double old_mad_delay  = state->delay_std_dev_ns;

    /* Update MAD: MAD_new = (1-alpha)*MAD_old + alpha*|new_value - old_mean|
     *
     * For the first few updates, use raw absolute value as initial estimate.
     * The running mean is approximated by the stored offset/delay values.
     */
    if (state->valid_sync_count == 0) {
        state->offset_std_dev_ns = fabs(offset_ns) * 0.1;
        state->delay_std_dev_ns = 0.1; /* Minimum estimate */
    } else {
        double offset_mean = state->offset_from_master_ns;
        double delay_mean  = state->mean_path_delay_ns;

        state->offset_std_dev_ns = (1.0 - alpha) * old_mad_offset
                                 + alpha * fabs(offset_ns - offset_mean);
        state->delay_std_dev_ns = (1.0 - alpha) * old_mad_delay
                                + alpha * fabs(delay_ns - delay_mean);
    }

    return 0;
}

/* ===================================================================
 * L5: BMCA Dataset Comparison (IEEE 1588 Section 10.3.8)
 * =================================================================== */

BmcaResult bmca_compare_datasets(const BmcaDataset *A, const BmcaDataset *B)
{
    if (!A || !B) return BMCA_EQUAL;

    /* Comparison order per IEEE 1588 Table 34 */

    /* 1. Priority1: lower is better */
    if (A->priority1 < B->priority1) return BMCA_THIS_IS_BETTER;
    if (A->priority1 > B->priority1) return BMCA_OTHER_IS_BETTER;

    /* 2. Clock class: lower is better */
    if (A->clock_quality.clock_class < B->clock_quality.clock_class)
        return BMCA_THIS_IS_BETTER;
    if (A->clock_quality.clock_class > B->clock_quality.clock_class)
        return BMCA_OTHER_IS_BETTER;

    /* 3. Clock accuracy: lower is better */
    if (A->clock_quality.clock_accuracy < B->clock_quality.clock_accuracy)
        return BMCA_THIS_IS_BETTER;
    if (A->clock_quality.clock_accuracy > B->clock_quality.clock_accuracy)
        return BMCA_OTHER_IS_BETTER;

    /* 4. Offset scaled log variance: lower is better */
    if (A->clock_quality.offset_scaled_log_variance <
        B->clock_quality.offset_scaled_log_variance)
        return BMCA_THIS_IS_BETTER;
    if (A->clock_quality.offset_scaled_log_variance >
        B->clock_quality.offset_scaled_log_variance)
        return BMCA_OTHER_IS_BETTER;

    /* 5. Priority2: lower is better */
    if (A->priority2 < B->priority2) return BMCA_THIS_IS_BETTER;
    if (A->priority2 > B->priority2) return BMCA_OTHER_IS_BETTER;

    /* 6. Clock identity: lower is better (tiebreaker) */
    for (int i = 0; i < 8; i++) {
        if (A->clock_identity.octets[i] < B->clock_identity.octets[i])
            return BMCA_THIS_IS_BETTER;
        if (A->clock_identity.octets[i] > B->clock_identity.octets[i])
            return BMCA_OTHER_IS_BETTER;
    }

    return BMCA_EQUAL;
}

double ptp_announce_interval(int8_t log_message_interval)
{
    /* IEEE 1588: announce interval = 2^log_message_interval seconds.
     * int8_t range [-128, 127] covers all valid PTP intervals. */
    (void)(log_message_interval + 128); /* Suppress unused when debug */
    return pow(2.0, (double)log_message_interval);
}

/* ===================================================================
 * L5: PTP Slave Initialization and Update
 * =================================================================== */

void ptp_slave_init(PtpSlaveState *state, double initial_offset_ns,
                    double update_interval_s)
{
    if (!state) return;
    (void)update_interval_s; /* Reserved for future servo tuning */

    memset(state, 0, sizeof(PtpSlaveState));

    /* Initialize PI servo */
    pi_servo_init(&state->servo_config, &state->servo_state);

    /* Initialize Kalman clock state */
    clock_state_init(&state->kalman_clock, initial_offset_ns,
                     0.0, 0.0,  /* Unknown freq and drift initially */
                     1000.0,     /* 1 us initial offset uncertainty */
                     10.0);      /* 10 ppb initial freq uncertainty */

    state->sync_status = SYNC_FREE_RUNNING;
    state->sync_count = 0;
    state->valid_sync_count = 0;
}

double ptp_slave_update(PtpSlaveState *state, const PtpTimestamps *ts)
{
    if (!state || !ts) return 0.0;

    /* 1. Validate timestamps */
    if (!ptp_validate_timestamps(ts)) {
        state->sync_status = SYNC_LOS;
        return 0.0;
    }

    /* 2. Store and process timestamps */
    state->last_timestamps = *ts;
    if (ptp_process_timestamps(state) != 0) {
        state->sync_status = SYNC_LOS;
        return 0.0;
    }

    state->sync_count++;

    /* 3. Update Kalman filter with measured offset */
    double meas_noise_var = state->offset_std_dev_ns * state->offset_std_dev_ns;
    if (meas_noise_var < 1.0) meas_noise_var = 1.0; /* Minimum 1 ns^2 */

    if (state->valid_sync_count == 0) {
        /* First valid measurement: initialize Kalman state */
        clock_state_init(&state->kalman_clock,
                         state->offset_from_master_ns,
                         0.0, 0.0,
                         meas_noise_var, 10.0);
    } else {
        clock_kalman_update(&state->kalman_clock,
                           state->offset_from_master_ns,
                           meas_noise_var);
    }

    state->valid_sync_count++;

    /* 4. Run PI servo to compute clock correction */
    double correction = pi_servo_update(&state->servo_config,
                                        &state->servo_state,
                                        state->offset_from_master_ns);

    /* 5. Update sync status */
    if (state->valid_sync_count < 3) {
        state->sync_status = SYNC_ACQUIRING;
    } else if (fabs(state->offset_from_master_ns) < 100.0) {
        state->sync_status = SYNC_LOCKED;
    } else {
        state->sync_status = SYNC_ACQUIRING;
    }

    return correction;
}

/* ===================================================================
 * L5: PTP Timestamp Validation
 * =================================================================== */

int ptp_validate_timestamps(const PtpTimestamps *ts)
{
    if (!ts) return 0;

    /* All timestamps must be valid */
    if (!timing_timestamp_valid(&ts->t1) ||
        !timing_timestamp_valid(&ts->t2) ||
        !timing_timestamp_valid(&ts->t3) ||
        !timing_timestamp_valid(&ts->t4)) {
        return 0;
    }

    /* Causality checks */
    /* t1 must precede t2 (sync message travels forward in time) */
    if (timing_timestamp_cmp(&ts->t1, &ts->t2) > 0) return 0;

    /* t3 must precede t4 (delay_req travels forward in time) */
    if (timing_timestamp_cmp(&ts->t3, &ts->t4) > 0) return 0;

    /* t1 should precede t3 (sync before delay_req in normal sequence) */
    if (timing_timestamp_cmp(&ts->t1, &ts->t3) > 0) return 0;

    /* t2 should precede t3 (slave receives sync before sending delay_req) */
    if (timing_timestamp_cmp(&ts->t2, &ts->t3) > 0) return 0;

    /* Check for excessive delays (more than 10 seconds round trip) */
    double delay_ns = timing_timestamp_diff_ns(&ts->t4, &ts->t1);
    if (delay_ns > 10.0e9) return 0;

    /* Check for negative apparent one-way times (indicates invalid timestamps) */
    double t2_m_t1 = timing_timestamp_diff_ns(&ts->t2, &ts->t1);
    double t4_m_t3 = timing_timestamp_diff_ns(&ts->t4, &ts->t3);
    if (t2_m_t1 < -1.0 || t4_m_t3 < -1.0) return 0;

    return 1;
}

/* ===================================================================
 * L6: PTP Sync Accuracy
 * =================================================================== */

double ptp_sync_accuracy(const PtpSlaveState *state)
{
    if (!state) return 0.0;

    /* Uncertainty budget (1-sigma):
     *   U_total^2 = U_timestamp^2 + U_delay_var^2 + U_osc_stab^2
     *
     * U_timestamp: from std dev of offset measurements
     * U_delay_var: from std dev of delay measurements
     * U_osc_stab: oscillator stability contribution
     */

    double U_ts = state->offset_std_dev_ns;
    double U_delay = state->delay_std_dev_ns;

    /* Oscillator stability: assume 1e-10 fractional frequency stability
     * over the sync interval -> time error ~ 1 ns */
    double U_osc = 1.0;

    return sqrt(U_ts * U_ts + U_delay * U_delay + U_osc * U_osc);
}

int ptp_within_accuracy(const PtpSlaveState *state, double accuracy_ns,
                        double confidence_k)
{
    if (!state) return 0;

    double U_total = ptp_sync_accuracy(state);
    return (U_total * confidence_k <= accuracy_ns) ? 1 : 0;
}

/* ===================================================================
 * L6: Next Sync Time
 * =================================================================== */

void ptp_next_sync_time(int8_t log_sync_interval,
                        const Timestamp *current_time,
                        Timestamp *next_time)
{
    if (!current_time || !next_time) return;

    double interval_s = pow(2.0, (double)log_sync_interval);
    double interval_ns = interval_s * 1.0e9;

    *next_time = *current_time;
    timing_timestamp_add_ns(next_time, (int64_t)interval_ns);
}

/* ===================================================================
 * L7: 5G Fronthaul Timing (ITU-T G.8271.1)
 * =================================================================== */

int ptp_5g_fronthaul_check(double offset_ns, char class_level)
{
    double abs_offset = fabs(offset_ns);

    switch (class_level) {
    case 'A':
    case 'a':
        /* Class A: +/- 1100 ns peak-to-peak -> |offset| < 550 ns */
        return (abs_offset <= 550.0) ? 1 : 0;

    case 'B':
    case 'b':
        /* Class B: +/- 260 ns peak-to-peak -> |offset| < 130 ns */
        return (abs_offset <= 130.0) ? 1 : 0;

    case 'C':
    case 'c':
        /* Class C: +/- 30 ns peak-to-peak -> |offset| < 15 ns */
        return (abs_offset <= 15.0) ? 1 : 0;

    default:
        return 0;
    }
}

/* ===================================================================
 * L7: IEC 61850 Power Grid Timing
 * =================================================================== */

int ptp_power_grid_check(double offset_ns)
{
    /* IEC 61850-9-2 Sampled Values:
     * Timing error must be < 1 us for protection class P1/P2/P3.
     *
     * IEC 61850-90-4 Synchrophasor:
     * Phase angle error < 0.1 degrees = 5.5 us at 50 Hz.
     * Timing error corresponds to phase error:
     *   phase_err = 360 * f * time_err
     *   0.1 deg = 360 * 50 * time_err -> time_err = 5.56 us
     */

    double abs_offset = fabs(offset_ns);

    /* Sampled Values requirement: < 1000 ns */
    if (abs_offset > 1000.0) return 0;

    /* Synchrophasor requirement: < 4000 ns (safety margin on 5.5 us) */
    if (abs_offset > 4000.0) return 0;

    return 1;
}
