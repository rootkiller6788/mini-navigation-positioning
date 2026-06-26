/**
 * mini-uwb-localization: UWB Ranging Implementation
 *
 * Complete implementations of Single-Sided TWR (SS-TWR) and
 * Double-Sided TWR (DS-TWR) ranging protocols, including
 * clock drift estimation and compensation, statistical burst
 * processing, and error budget analysis.
 *
 * Knowledge Coverage:
 *   L2 Core Concepts: Time-of-Flight, TWR protocol states
 *   L4 Fundamental Laws: CRLB for time delay, Friis equation
 *   L5 Algorithms: SS-TWR, DS-TWR, skewness-kurtosis NLOS test,
 *                  EWMA filtering, trimmed mean
 *   L6 Canonical Problems: Precision ranging with clock correction
 */

#include "uwb_ranging.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* =========================================================================
 * Helper: comparison function for qsort (double ascending)
 * ========================================================================= */

static int compare_double(const void *a, const void *b)
{
    double da = *(const double *)a;
    double db = *(const double *)b;
    if (da < db) return -1;
    if (da > db) return 1;
    return 0;
}

/* =========================================================================
 * TWR Session Management
 * ========================================================================= */

void twr_session_init(twr_session_t *session, uint16_t initiator_id,
                      uint16_t responder_id)
{
    if (!session) return;
    memset(session, 0, sizeof(*session));
    session->state = TWR_STATE_IDLE;
    session->initiator_id = initiator_id;
    session->responder_id = responder_id;
    session->reply_delay_us = UWB_TWR_DEFAULT_REPLY_DELAY_US;
    session->timeout_us = 20000;
    session->burst_size = UWB_TWR_DEFAULT_BURST_SIZE;
    session->antenna_delay_initiator_ps = 16450.0;
    session->antenna_delay_responder_ps = 16450.0;
}

/* =========================================================================
 * Single-Sided Two-Way Ranging (SS-TWR)
 *
 * Protocol (L2 Core Concept):
 *   1. Initiator sends Poll
 *   2. Responder receives Poll, waits T_reply, sends Response
 *   3. Initiator receives Response
 *
 * SS-TWR ToF Equation (L4):
 *   T_tof = (T_round - T_reply) / 2
 *
 * where:
 *   T_round = t_resp_rx - t_poll_tx  (measured on initiator)
 *   T_reply = t_resp_tx - t_poll_rx  (measured on responder)
 *
 * Error Analysis (L4 Fundamental Law):
 *   The primary SS-TWR error source is clock frequency offset.
 *   If initiator clock runs at f_A and responder at f_B:
 *     measured T_round = T_round_true * (1 + e_A)
 *     measured T_reply = T_reply_true * (1 + e_B)
 *     estimated T_tof = T_tof_true + T_reply*(e_A - e_B)/2 + ...
 *
 *   For a 5 ms reply time and 20 ppm offset difference:
 *     error = 5e-3 * 20e-6 * 3e8 / 2 = 15 mm
 *
 *   For large reply times, the error can exceed 1m ！ motivating DS-TWR.
 *
 * Reference: Decawave APS013 Application Note, Section 2.1
 */
double twr_ss_compute_tof(const twr_ss_timestamps_t *ts, double tick_period_ps)
{
    uint64_t T_round, T_reply;
    double tof;

    if (!ts || tick_period_ps <= 0.0) return -1.0;

    /* Guard against timestamp wraparound (40-bit timers on DW1000) */
    if (ts->t_resp_rx < ts->t_poll_tx || ts->t_resp_tx < ts->t_poll_rx) {
        return -1.0; /* Invalid timestamp ordering */
    }

    T_round = ts->t_resp_rx - ts->t_poll_tx;
    T_reply = ts->t_resp_tx - ts->t_poll_rx;

    /* T_tof = (T_round - T_reply) / 2, converted to seconds */
    tof = (double)(T_round - T_reply) * tick_period_ps * 1e-12 / 2.0;

    /* Sanity check: ToF must be positive and reasonable (< 1us) */
    if (tof < 0.0 || tof > 1e-6) {
        return -1.0;
    }

    return tof;
}

double twr_ss_compute_distance(const twr_ss_timestamps_t *ts,
                               double tick_period_ps,
                               double ant_delay_init_ps,
                               double ant_delay_resp_ps)
{
    double tof, corrected_tof;

    tof = twr_ss_compute_tof(ts, tick_period_ps);
    if (tof < 0.0) return -1.0;

    /* Compensate for antenna delays: subtract from measured ToF */
    corrected_tof = tof - (ant_delay_init_ps + ant_delay_resp_ps) * 1e-12 / 2.0;
    if (corrected_tof < 0.0) return 0.0;

    return UWB_C * corrected_tof / 2.0;
}

/* =========================================================================
 * Double-Sided Two-Way Ranging (DS-TWR)
 *
 * Protocol (L2 Core Concept):
 *   1. Initiator sends Poll
 *   2. Responder sends Response (after T_reply1)
 *   3. Initiator sends Final (after T_reply2)
 *   4. Responder receives Final
 *
 * DS-TWR ToF Equation (L4 Fundamental Law):
 *   T_tof = (T_round1 * T_round2 - T_reply1 * T_reply2) /
 *           (T_round1 + T_round2 + T_reply1 + T_reply2)
 *
 * where:
 *   T_round1 = t_resp_rx - t_poll_tx   (initiator side, round 1)
 *   T_reply1 = t_resp_tx - t_poll_rx   (responder side, reply 1)
 *   T_round2 = t_final_rx - t_resp_tx  (responder side, round 2)
 *   T_reply2 = t_final_tx - t_resp_rx  (initiator side, reply 2)
 *
 * Why DS-TWR eliminates first-order clock error (L4):
 *   Let e_A, e_B be clock frequency offsets (ppm/1e6).
 *   measured times: T' = T * (1 + e)
 *
 *   SS-TWR: T_tof_est = (T_rnd*(1+e_A) - T_rpy*(1+e_B))/2
 *           = T_tof_true + (T_rnd*e_A - T_rpy*e_B)/2
 *           error ~ T_reply * (e_A - e_B) / 2  [first-order]
 *
 *   DS-TWR: the symmetric product form cancels the (e_A+e_B) terms,
 *           leaving only second-order error ~ (e_A - e_B)^2
 *
 *   For 20 ppm offset: SS-TWR error ~ 150 mm, DS-TWR error ~ 0.03 mm
 *
 * Reference: Neirynck et al. (2016) "An Alternative Double-Sided
 *            Two-Way Ranging Method", MDPI Sensors
 */
double twr_ds_compute_tof(const twr_ds_timestamps_t *ts, double tick_period_ps)
{
    uint64_t T_round1, T_reply1, T_round2, T_reply2;
    double T1, T2, T3, T4;
    double numerator, denominator;
    double tof;

    if (!ts || tick_period_ps <= 0.0) return -1.0;

    /* Validate timestamp ordering */
    if (ts->t_resp_rx < ts->t_poll_tx || ts->t_resp_tx < ts->t_poll_rx ||
        ts->t_final_rx < ts->t_resp_tx || ts->t_final_tx < ts->t_resp_rx) {
        return -1.0;
    }

    T_round1 = ts->t_resp_rx - ts->t_poll_tx;
    T_reply1 = ts->t_resp_tx - ts->t_poll_rx;
    T_round2 = ts->t_final_rx - ts->t_resp_tx;
    T_reply2 = ts->t_final_tx - ts->t_resp_rx;

    /* Convert to seconds */
    T1 = (double)T_round1 * tick_period_ps * 1e-12;
    T2 = (double)T_reply1 * tick_period_ps * 1e-12;
    T3 = (double)T_round2 * tick_period_ps * 1e-12;
    T4 = (double)T_reply2 * tick_period_ps * 1e-12;

    /* DS-TWR symmetric formula */
    numerator   = T1 * T3 - T2 * T4;
    denominator = T1 + T2 + T3 + T4;

    if (denominator <= 0.0) return -1.0;

    tof = numerator / denominator;

    if (tof < 0.0 || tof > 1e-6) return -1.0;

    return tof;
}

double twr_ds_compute_distance(const twr_ds_timestamps_t *ts,
                               double tick_period_ps,
                               double ant_delay_init_ps,
                               double ant_delay_resp_ps)
{
    double tof, corrected_tof;

    tof = twr_ds_compute_tof(ts, tick_period_ps);
    if (tof < 0.0) return -1.0;

    corrected_tof = tof - (ant_delay_init_ps + ant_delay_resp_ps) * 1e-12 / 2.0;
    if (corrected_tof < 0.0) return 0.0;

    return UWB_C * corrected_tof / 2.0;
}

/* =========================================================================
 * Clock Offset Estimation from DS-TWR (L5 Algorithm)
 *
 * The ratio of times measured on different devices reveals the
 * clock frequency ratio:
 *
 *   k = f_initiator / f_responder
 *
 * From the DS-TWR timestamps:
 *   k 「 (T_round1 + T_reply2) / (T_reply1 + T_round2)
 *
 * Proof (L4): Let true time be t, measured times:
 *   t'_A = t * (1 + e_A) on initiator
 *   t'_B = t * (1 + e_B) on responder
 *
 * Then:
 *   T_round1 measured on A = T_round1_true * (1 + e_A)
 *   T_reply1 measured on B = T_reply1_true * (1 + e_B)
 *   T_round2 measured on B = T_round2_true * (1 + e_B)
 *   T_reply2 measured on A = T_reply2_true * (1 + e_A)
 *
 * Ratio of sums:
 *   (T_round1_A + T_reply2_A) / (T_reply1_B + T_round2_B)
 *   = (1 + e_A) / (1 + e_B) * (T_rnd1_true + T_rpy2_true) / (T_rpy1_true + T_rnd2_true)
 *
 * For equal reply times (common in symmetric DS-TWR):
 *   T_rnd1_true + T_rpy2_true = T_rpy1_true + T_rnd2_true + 4*T_tof
 *   「 T_rpy1_true + T_rnd2_true   (since T_tof << reply times)
 *
 * Therefore: k 「 (1 + e_A) / (1 + e_B)
 *   clock_offset_ppm = (k - 1) * 1e6 「 (e_A - e_B) * 1e6
 */
double twr_ds_estimate_clock_offset(twr_ds_timestamps_t *ts, double tick_period_ps)
{
    (void)tick_period_ps;
    (void)tick_period_ps;
    double sum_A, sum_B, k;

    if (!ts) return 0.0;

    sum_A = (double)((ts->t_resp_rx - ts->t_poll_tx) +
                     (ts->t_final_tx - ts->t_resp_rx));
    sum_B = (double)((ts->t_resp_tx - ts->t_poll_rx) +
                     (ts->t_final_rx - ts->t_resp_tx));

    if (sum_B <= 0.0) return 0.0;

    k = sum_A / sum_B;
    ts->clock_offset_ppm = (k - 1.0) * 1e6;

    return ts->clock_offset_ppm;
}

double twr_compensate_clock_offset(double distance_m, double clock_offset_ppm)
{
    /* d_corrected = d_measured / (1 + offset_ppm * 1e-6) */
    return distance_m / (1.0 + clock_offset_ppm * 1e-6);
}

double twr_estimate_clock_error_m(double reply_time_us, double clock_offset_ppm)
{
    /* error 「 T_reply * c * |offset| / 2 */
    return reply_time_us * 1e-6 * UWB_C * fabs(clock_offset_ppm) * 1e-6 / 2.0;
}

/* =========================================================================
 * Ranging Burst Statistics (L5: Statistical Filtering)
 *
 * Statistical processing of multiple rapid ranging measurements to
 * improve accuracy and reject outliers.
 *
 * Trimmed Mean (L5 Algorithm):
 *   Sort N distances, discard alpha*N from each end,
 *   compute mean of remaining (1-2*alpha)*N values.
 *   - Robust to outliers (unlike mean)
 *   - More efficient than median for large N
 *   - Typical alpha = 0.25 (trim 25% each end)
 *
 * Skewness Test for NLOS (L5 Algorithm):
 *   skewness = E[(x - mu)^3] / sigma^3
 *   - LOS: symmetric distribution, skewness ~ 0
 *   - NLOS: positive skewness (delays can only increase distance)
 *
 * Kurtosis Test for NLOS:
 *   kurtosis = E[(x - mu)^4] / sigma^4
 *   - LOS: kurtosis ~ 3 (Gaussian-like)
 *   - NLOS: kurtosis > 3 (heavy-tailed due to large positive errors)
 *
 * Reference: Guvenc et al. (2007) "NLOS Identification for UWB Localization"
 * Reference: Huber (1981) "Robust Statistics"
 */
void ranging_burst_compute_stats(const double *distances, int count,
                                 ranging_burst_stats_t *stats)
{
    double *sorted;
    double sum, sum_sq;
    int i, trim_count, start_idx, end_idx;

    if (!distances || count <= 0 || !stats) return;

    memset(stats, 0, sizeof(*stats));
    stats->count = count;

    /* Min, max, mean */
    stats->min_val = distances[0];
    stats->max_val = distances[0];
    sum = 0.0;
    for (i = 0; i < count; i++) {
        if (distances[i] < stats->min_val) stats->min_val = distances[i];
        if (distances[i] > stats->max_val) stats->max_val = distances[i];
        sum += distances[i];
    }
    stats->mean = sum / count;

    /* Standard deviation */
    sum_sq = 0.0;
    for (i = 0; i < count; i++) {
        double diff = distances[i] - stats->mean;
        sum_sq += diff * diff;
    }
    stats->stddev = (count > 1) ? sqrt(sum_sq / (count - 1)) : 0.0;

    /* Median (requires sorted copy) */
    sorted = (double *)malloc(count * sizeof(double));
    if (sorted) {
        memcpy(sorted, distances, count * sizeof(double));
        qsort(sorted, count, sizeof(double), compare_double);
        if (count % 2 == 1) {
            stats->median = sorted[count / 2];
        } else {
            stats->median = (sorted[count / 2 - 1] + sorted[count / 2]) / 2.0;
        }

        /* Trimmed mean (trim UWB_TWR_TRIMMED_MEAN_RATIO from each end) */
        trim_count = (int)(count * UWB_TWR_TRIMMED_MEAN_RATIO);
        if (trim_count < 0) trim_count = 0;
        if (trim_count > count / 3) trim_count = count / 3;

        start_idx = trim_count;
        end_idx = count - trim_count;
        if (end_idx > start_idx) {
            double trim_sum = 0.0;
            int trim_n = 0;
            for (i = start_idx; i < end_idx; i++) {
                trim_sum += sorted[i];
                trim_n++;
            }
            stats->trimmed_mean = (trim_n > 0) ? trim_sum / trim_n : stats->mean;
        } else {
            stats->trimmed_mean = stats->mean;
        }

        free(sorted);
    }
}

int ranging_burst_median_filter(double *distances, int count, int window_size)
{
    double *temp;
    int i, j, k, half_window, valid_count;

    if (!distances || count <= 0 || window_size < 3 || window_size % 2 == 0) {
        return 0;
    }

    temp = (double *)malloc(count * sizeof(double));
    if (!temp) return 0;

    half_window = window_size / 2;
    valid_count = 0;

    for (i = 0; i < count; i++) {
        /* Collect window samples */
        int window_start = i - half_window;
        int window_end = i + half_window;
        double window_vals[9]; /* max window size 9 */
        int wcount = 0;

        if (window_start < 0) window_start = 0;
        if (window_end >= count) window_end = count - 1;

        for (j = window_start; j <= window_end && wcount < 9; j++) {
            window_vals[wcount++] = distances[j];
        }

        /* Sort window and take median */
        for (j = 0; j < wcount - 1; j++) {
            for (k = j + 1; k < wcount; k++) {
                if (window_vals[k] < window_vals[j]) {
                    double tmp = window_vals[j];
                    window_vals[j] = window_vals[k];
                    window_vals[k] = tmp;
                }
            }
        }

        temp[valid_count++] = window_vals[wcount / 2];
    }

    memcpy(distances, temp, valid_count * sizeof(double));
    free(temp);

    return valid_count;
}

double ranging_burst_nlos_score(const double *distances, int count)
{
    double mean, var, std, skewness, kurtosis;
    double m3, m4;
    double score;
    int i;

    if (!distances || count < 4) return 0.0;

    /* Compute mean */
    mean = 0.0;
    for (i = 0; i < count; i++) mean += distances[i];
    mean /= count;

    /* Compute variance, third moment, fourth moment */
    var = 0.0;
    m3 = 0.0;
    m4 = 0.0;
    for (i = 0; i < count; i++) {
        double diff = distances[i] - mean;
        double diff2 = diff * diff;
        var += diff2;
        m3 += diff * diff2;
        m4 += diff2 * diff2;
    }

    if (var < 1e-20) return 0.0;

    var /= count;
    m3 /= count;
    m4 /= count;

    std = sqrt(var);

    /* Skewness = m3 / sigma^3 */
    skewness = m3 / (std * std * std);

    /* Kurtosis = m4 / sigma^4 (excess kurtosis = kurtosis - 3) */
    kurtosis = m4 / (var * var);

    /* NLOS score from Guvenc et al. (2007):
     * score = logistic(skewness, kurtosis)
     * Higher skewness (>0.5) and kurtosis (>3.5) indicate NLOS
     */
    {
        double z_skew = (skewness - 0.1) / 0.8;   /* Normalize skewness */
        double z_kurt = (kurtosis - 3.0) / 2.0;    /* Normalize excess kurtosis */
        double logit;

        if (z_skew < 0) z_skew = 0;
        if (z_kurt < 0) z_kurt = 0;

        /* Weighted combination */
        logit = -2.0 + 1.5 * z_skew + 0.8 * z_kurt;

        /* Logistic function to map to [0, 1] */
        score = 1.0 / (1.0 + exp(-logit));

        /* Clamp */
        if (score > 1.0) score = 1.0;
        if (score < 0.0) score = 0.0;
    }

    return score;
}

double ranging_ewma_filter(const double *history, int count, double alpha)
{
    double filtered;
    int i;

    if (!history || count <= 0) return 0.0;
    if (alpha <= 0.0) alpha = 0.1;
    if (alpha >= 1.0) alpha = 0.9;

    filtered = history[0];

    /* EWMA: s_t = alpha * x_t + (1 - alpha) * s_{t-1} */
    for (i = 1; i < count; i++) {
        filtered = alpha * history[i] + (1.0 - alpha) * filtered;
    }

    return filtered;
}

/* =========================================================================
 * Ranging Error Budget (L4 Fundamental Law)
 *
 * Total RSS error from independent error sources:
 *
 * sigma_total^2 = sigma_noise^2 + sigma_clock^2 + sigma_mp^2 + sigma_ant^2
 *
 * 1. Thermal noise floor (from CRLB):
 *    sigma_noise = c / (2 * pi * sqrt(2) * sqrt(SNR) * B_eff)
 *                 = c / (2 * pi * B_eff * sqrt(2 * SNR))
 *
 * 2. Clock offset error (SS-TWR dominant term):
 *    sigma_clock = T_reply * c * ppm_error * 1e-6 / 2
 *
 * 3. Multipath error (empirical model):
 *    Multipath causes biased range estimates. The error depends on
 *    the relative strength and delay of reflected paths vs. direct path.
 *    For UWB with leading-edge detection: sigma_mp is typically < 0.3m
 *    in residential/office environments.
 *
 * 4. Antenna delay calibration residual:
 *    Typical factory calibration: +/- 2 cm
 *    Can be reduced to < 5 mm with per-device calibration.
 *
 * The RSS combination assumes independent, zero-mean errors.
 * This is valid for noise and clock errors but multipath bias
 * requires separate treatment (see NLOS module).
 */
double twr_error_budget(double snr_linear, double bw_effective,
                        double reply_time_us, double clock_ppm_error,
                        double multipath_error_m, double ant_cal_error_m)
{
    double sigma_noise, sigma_clock;
    double total_var;

    /* 1. Thermal noise */
    sigma_noise = uwb_crlb_distance(snr_linear, bw_effective);

    /* 2. Clock offset */
    sigma_clock = twr_estimate_clock_error_m(reply_time_us, clock_ppm_error);

    /* 3. RSS combination */
    total_var = sigma_noise * sigma_noise
              + sigma_clock * sigma_clock
              + multipath_error_m * multipath_error_m
              + ant_cal_error_m * ant_cal_error_m;

    return sqrt(total_var);
}
