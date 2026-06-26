/**
 * @file tof_tdoa_positioning.c
 * @brief Time-of-Flight and TDOA positioning implementations
 *
 * Implements: TWR/SDS-TWR distance, TOA Gauss-Newton, TDOA Taylor series,
 * NLOS detection/mitigation, AoA estimation/triangulation,
 * UWB link budget, first path detection.
 */

#include <math.h>
#include <string.h>
#include <stdio.h>
#include "../include/tof_tdoa_positioning.h"
#include "../include/indoor_positioning.h"

/* ============================================================================
 * L5 - Two-Way Ranging Distance Estimation
 * ============================================================================ */

double twr_compute_distance(const twr_exchange_t *twr) {
    if (!twr) return -1.0;
    if (twr->tick_period_s <= 0.0) return -1.0;

    double t_round = (twr->t_rr - twr->t_sp) * twr->tick_period_s;
    double t_reply = (twr->t_sr - twr->t_rp) * twr->tick_period_s;

    /* Standard TWR: d = (t_round - t_reply) * c / 2 */
    double tof = (t_round - t_reply) / 2.0;
    if (tof < 0.0) return -1.0;  /* Invalid exchange */

    double distance = tof * SPEED_OF_LIGHT_MPS;
    if (distance > UWB_MAX_RANGE) distance = UWB_MAX_RANGE;
    return distance;
}

double twr_sds_compute_distance(const twr_exchange_t *twr) {
    if (!twr || !twr->is_ss_twr) return -1.0;
    if (twr->tick_period_s <= 0.0) return -1.0;

    /* SDS-TWR uses three messages to cancel clock drift.
     * T_round1 = t_rr - t_sp (round trip at initiator)
     * T_reply1 = t_sr - t_rp (reply time at responder)
     * T_round2 = t_rf - t_sr (round trip at responder)
     * T_reply2 = t_sf - t_rr (reply time at initiator)
     *
     * TOF = (T_round1 * T_round2 - T_reply1 * T_reply2) /
     *       (T_round1 + T_round2 + T_reply1 + T_reply2)
     */
    double T_round1 = (twr->t_rr - twr->t_sp) * twr->tick_period_s;
    double T_reply1 = (twr->t_sr - twr->t_rp) * twr->tick_period_s;
    double T_round2 = (twr->t_rf - twr->t_sr) * twr->tick_period_s;
    double T_reply2 = (twr->t_sf - twr->t_rr) * twr->tick_period_s;

    if (T_round1 <= 0.0 || T_reply1 <= 0.0 || T_round2 <= 0.0 || T_reply2 <= 0.0) {
        return -1.0;
    }

    double num = T_round1 * T_round2 - T_reply1 * T_reply2;
    double den = T_round1 + T_round2 + T_reply1 + T_reply2;

    if (den <= 0.0) return -1.0;

    double tof = num / den;
    if (tof < 0.0) return -1.0;

    double distance = tof * SPEED_OF_LIGHT_MPS;
    if (distance > UWB_MAX_RANGE) distance = UWB_MAX_RANGE;
    return distance;
}

/* ============================================================================
 * L5 - Gauss-Newton TOA Positioning (2D)
 * ============================================================================ */

int toa_positioning_2d(const position2d_t *anchor_positions,
                       const double *toa_distances,
                       int n_anchors,
                       const position2d_t *initial_guess,
                       position2d_t *result,
                       int max_iter, double tol) {
    if (!anchor_positions || !toa_distances || !result || !initial_guess) return -1;
    if (n_anchors < 3) return -1;
    if (max_iter <= 0) max_iter = 30;
    if (tol <= 0.0) tol = 1e-4;

    double x = initial_guess->x;
    double y = initial_guess->y;

    for (int iter = 0; iter < max_iter; iter++) {
        /* Build Jacobian H (Nx2) and residuals */
        double HtH[2][2] = {{0,0},{0,0}};
        double Htr[2] = {0,0};
        double max_residual = 0.0;

        for (int i = 0; i < n_anchors; i++) {
            double dx = x - anchor_positions[i].x;
            double dy = y - anchor_positions[i].y;
            double dist_pred = sqrt(dx*dx + dy*dy);
            if (dist_pred < 1e-10) dist_pred = 1e-10;

            /* Jacobian row: h = [dx/dist, dy/dist] */
            double h0 = dx / dist_pred;
            double h1 = dy / dist_pred;
            double residual = toa_distances[i] - dist_pred;

            HtH[0][0] += h0 * h0;
            HtH[0][1] += h0 * h1;
            HtH[1][0] += h0 * h1;
            HtH[1][1] += h1 * h1;
            Htr[0] += h0 * residual;
            Htr[1] += h1 * residual;

            if (fabs(residual) > max_residual) max_residual = fabs(residual);
        }

        if (max_residual < tol) break;

        /* Solve HtH * delta = Htr */
        double det = HtH[0][0]*HtH[1][1] - HtH[0][1]*HtH[1][0];
        if (fabs(det) < 1e-15) break;

        double dx = (Htr[0]*HtH[1][1] - Htr[1]*HtH[0][1]) / det;
        double dy = (HtH[0][0]*Htr[1] - HtH[0][1]*Htr[0]) / det;

        x += dx;
        y += dy;

        if (sqrt(dx*dx + dy*dy) < tol) break;
    }

    result->x = x;
    result->y = y;
    return 0;
}

/* ============================================================================
 * L5 - Taylor Series TDOA Positioning
 * ============================================================================ */

int tdoa_taylor_series(const position3d_t *anchors,
                       const double *tdoa_s,
                       int n_anchors,
                       double speed,
                       const position3d_t *initial_guess,
                       position3d_t *result,
                       int max_iter, double tol) {
    if (!anchors || !tdoa_s || !result || !initial_guess) return -1;
    if (n_anchors < 4) return -1;
    if (max_iter <= 0) max_iter = 20;
    if (tol <= 0.0) tol = 1e-4;
    if (speed <= 0.0) speed = SPEED_OF_LIGHT_MPS;

    double x = initial_guess->x;
    double y = initial_guess->y;
    double z = initial_guess->z;

    /* Reference anchor is anchors[0] */
    double x0 = anchors[0].x, y0 = anchors[0].y, z0 = anchors[0].z;

    for (int iter = 0; iter < max_iter; iter++) {
        /* Distance to reference anchor */
        double r0 = sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0) + (z-z0)*(z-z0));
        if (r0 < 1e-10) r0 = 1e-10;

        /* Build H (Mx3) and residual b (M) where M = n_anchors - 1 */
        double HtH[3][3] = {{0}};
        double Htb[3] = {0};
        double max_correction = 0.0;
        int M = n_anchors - 1;
        if (M > 0) M = (M > 15) ? 15 : M;

        for (int i = 0; i < M; i++) {
            int idx = i + 1;
            double xi = anchors[idx].x, yi = anchors[idx].y, zi = anchors[idx].z;
            double ri = sqrt((x-xi)*(x-xi) + (y-yi)*(y-yi) + (z-zi)*(z-zi));
            if (ri < 1e-10) ri = 1e-10;

            /* Measured range difference: d_i1 = c * tdoa_{i,1} */
            double d_meas = speed * tdoa_s[i];
            /* Predicted range difference */
            double d_pred = ri - r0;
            /* Residual */
            double b_i = d_meas - d_pred;

            /* Jacobian row for TDOA: H_i = [ (x-xi)/ri - (x-x0)/r0 ,
             *                                (y-yi)/ri - (y-y0)/r0 ,
             *                                (z-zi)/ri - (z-z0)/r0 ] */
            double h0 = (x - xi) / ri - (x - x0) / r0;
            double h1 = (y - yi) / ri - (y - y0) / r0;
            double h2 = (z - zi) / ri - (z - z0) / r0;

            HtH[0][0] += h0*h0; HtH[0][1] += h0*h1; HtH[0][2] += h0*h2;
            HtH[1][0] += h1*h0; HtH[1][1] += h1*h1; HtH[1][2] += h1*h2;
            HtH[2][0] += h2*h0; HtH[2][1] += h2*h1; HtH[2][2] += h2*h2;

            Htb[0] += h0 * b_i;
            Htb[1] += h1 * b_i;
            Htb[2] += h2 * b_i;

            if (fabs(b_i) > max_correction) max_correction = fabs(b_i);
        }

        if (max_correction < tol) break;

        /* Solve 3x3 system */
        double det = HtH[0][0]*(HtH[1][1]*HtH[2][2] - HtH[1][2]*HtH[2][1])
                   - HtH[0][1]*(HtH[1][0]*HtH[2][2] - HtH[1][2]*HtH[2][0])
                   + HtH[0][2]*(HtH[1][0]*HtH[2][1] - HtH[1][1]*HtH[2][0]);
        if (fabs(det) < 1e-15) break;

        double dx = (Htb[0]*(HtH[1][1]*HtH[2][2]-HtH[1][2]*HtH[2][1])
                   - HtH[0][1]*(Htb[1]*HtH[2][2]-HtH[1][2]*Htb[2])
                   + HtH[0][2]*(Htb[1]*HtH[2][1]-HtH[1][1]*Htb[2])) / det;
        double dy = (HtH[0][0]*(Htb[1]*HtH[2][2]-HtH[1][2]*Htb[2])
                   - Htb[0]*(HtH[1][0]*HtH[2][2]-HtH[1][2]*HtH[2][0])
                   + HtH[0][2]*(HtH[1][0]*Htb[2]-Htb[1]*HtH[2][0])) / det;
        double dz = (HtH[0][0]*(HtH[1][1]*Htb[2]-Htb[1]*HtH[2][1])
                   - HtH[0][1]*(HtH[1][0]*Htb[2]-Htb[1]*HtH[2][0])
                   + Htb[0]*(HtH[1][0]*HtH[2][1]-HtH[1][1]*HtH[2][0])) / det;

        x += dx; y += dy; z += dz;
        if (sqrt(dx*dx + dy*dy + dz*dz) < tol) break;
    }

    result->x = x; result->y = y; result->z = z;
    return 0;
}

/* ============================================================================
 * L5 - NLOS Detection
 * ============================================================================ */

int detect_nlos_rssi_distance(double measured_distance, double measured_rssi,
                              const path_loss_model_t *model,
                              double n_sigma_thresh) {
    if (!model) return 0;
    if (measured_distance <= 0.0) return 1;  /* Invalid distance → NLOS */

    double expected_rssi = distance_to_rssi(measured_distance, model);
    double deviation = expected_rssi - measured_rssi;
    /* If measured RSSI is significantly lower than expected → NLOS */
    if (deviation > n_sigma_thresh * model->shadow_std) {
        return 1;
    }
    return 0;
}

void compute_nlos_weights(const uwb_ranging_t *range_measurements,
                          int n_measurements, double *weights) {
    if (!range_measurements || !weights || n_measurements <= 0) return;
    for (int i = 0; i < n_measurements; i++) {
        /* Weight = LOS confidence [0, 1] */
        weights[i] = range_measurements[i].los_confidence;
        if (weights[i] < 0.01) weights[i] = 0.01;  /* Minimum weight */
    }
}

int detect_nlos_residual(const position2d_t *anchor_positions,
                         const double *distances,
                         int n_anchors,
                         const position2d_t *estimated_pos,
                         double *residuals,
                         double threshold,
                         int *nlos_flags) {
    if (!anchor_positions || !distances || !estimated_pos
        || !residuals || !nlos_flags) return -1;

    int n_nlos = 0;
    for (int i = 0; i < n_anchors; i++) {
        double dx = estimated_pos->x - anchor_positions[i].x;
        double dy = estimated_pos->y - anchor_positions[i].y;
        double d_pred = sqrt(dx*dx + dy*dy);
        residuals[i] = fabs(distances[i] - d_pred);
        if (residuals[i] > threshold) {
            nlos_flags[i] = 1;
            n_nlos++;
        } else {
            nlos_flags[i] = 0;
        }
    }
    return n_nlos;
}

/* ============================================================================
 * L5 - Angle of Arrival Estimation
 * ============================================================================ */

double aoa_from_phase(const aoa_measurement_t *meas) {
    if (!meas) return NAN;
    if (meas->wavelength_m <= 0.0 || meas->antenna_spacing_m <= 0.0) return NAN;

    /* theta = arcsin( delta_phi * lambda / (2 * pi * d) ) */
    double arg = meas->phase_difference_rad * meas->wavelength_m
                 / (2.0 * M_PI * meas->antenna_spacing_m);

    /* Clamp to [-1, 1] for arcsin domain */
    if (arg > 1.0) arg = 1.0;
    if (arg < -1.0) arg = -1.0;

    return asin(arg);
}

int aoa_triangulate(const position2d_t *anchor1, const position2d_t *anchor2,
                    double aoa1, double aoa2,
                    position2d_t *result) {
    if (!anchor1 || !anchor2 || !result) return -1;

    /* Represent lines of bearing:
     * L1: through anchor1 with direction (sin(aoa1), cos(aoa1)) in ENU
     *     Actually in standard geometry: aoa from East, CCW
     *     But for indoor positioning, often aoa from North, CW.
     *     We'll use: direction = (sin(aoa), cos(aoa)) for aoa from North CW.
     *
     * Line 1: (x,y) = (x1, y1) + t * (sin(aoa1), cos(aoa1))
     * Line 2: (x,y) = (x2, y2) + s * (sin(aoa2), cos(aoa2))
     *
     * Intersection: solve for t, s
     */
    double dx1 = sin(aoa1), dy1 = cos(aoa1);
    double dx2 = sin(aoa2), dy2 = cos(aoa2);

    /* Determinant of direction vectors: dx1*dy2 - dy1*dx2 */
    double det = dx1 * dy2 - dy1 * dx2;
    if (fabs(det) < 1e-12) {
        /* Lines are parallel or nearly parallel */
        return -1;
    }

    double x12 = anchor2->x - anchor1->x;
    double y12 = anchor2->y - anchor1->y;

    double t = (x12 * dy2 - y12 * dx2) / det;

    result->x = anchor1->x + t * dx1;
    result->y = anchor1->y + t * dy1;
    return 0;
}

int aoa_positioning_stansfield(const position2d_t *anchor_positions,
                               const double *aoa_measurements,
                               int n_anchors,
                               position2d_t *result) {
    if (!anchor_positions || !aoa_measurements || !result) return -1;
    if (n_anchors < 2) return -1;

    /* Stansfield (pseudolinear) estimator:
     * For each anchor i with angle aoa_i:
     *   sin(aoa_i)*x - cos(aoa_i)*y = sin(aoa_i)*x_i - cos(aoa_i)*y_i
     *
     * This gives N linear equations. Solve via least squares.
     */
    double AtA[2][2] = {{0,0},{0,0}};
    double Atb[2] = {0,0};

    for (int i = 0; i < n_anchors; i++) {
        double s = sin(aoa_measurements[i]);
        double c = cos(aoa_measurements[i]);

        double ai0 = s;      /* coefficient for x */
        double ai1 = -c;     /* coefficient for y */
        double bi = s * anchor_positions[i].x - c * anchor_positions[i].y;

        AtA[0][0] += ai0 * ai0;
        AtA[0][1] += ai0 * ai1;
        AtA[1][0] += ai0 * ai1;
        AtA[1][1] += ai1 * ai1;
        Atb[0] += ai0 * bi;
        Atb[1] += ai1 * bi;
    }

    double det = AtA[0][0]*AtA[1][1] - AtA[0][1]*AtA[1][0];
    if (fabs(det) < 1e-15) return -1;

    result->x = (Atb[0]*AtA[1][1] - Atb[1]*AtA[0][1]) / det;
    result->y = (AtA[0][0]*Atb[1] - AtA[0][1]*Atb[0]) / det;
    return 0;
}

/* ============================================================================
 * L6 - UWB Link Budget and Ranging Precision
 * ============================================================================ */

double uwb_link_budget(double tx_power_dbm, double tx_antenna_gain_dbi,
                       double rx_antenna_gain_dbi, double distance,
                       double frequency_mhz, double rx_sensitivity_dbm) {
    if (distance <= 0.0 || frequency_mhz <= 0.0) return -1e9;

    /* Free-space path loss (Friis):
     * FSPL(dB) = 20*log10(d) + 20*log10(f) + 20*log10(4*pi/c)
     *          = 20*log10(d) + 20*log10(f) - 147.55  (f in Hz, d in m)
     * For f in MHz: 20*log10(f_MHz) + 20*log10(1e6) = 20*log10(f_MHz) + 120
     * FSPL = 20*log10(d) + 20*log10(f_MHz) + 120 - 147.55
     *      = 20*log10(d) + 20*log10(f_MHz) - 27.55
     */
    double fspl = 20.0 * log10(distance) + 20.0 * log10(frequency_mhz) - 27.55;
    double rx_power = tx_power_dbm + tx_antenna_gain_dbi
                      + rx_antenna_gain_dbi - fspl;
    return rx_power - rx_sensitivity_dbm;
}

double uwb_ranging_crlb(double bandwidth_hz, double snr_linear) {
    if (bandwidth_hz <= 0.0 || snr_linear <= 0.0) return 1e9;
    /* CRLB for TOA: sigma_d >= c / (2*sqrt(2)*pi*B*sqrt(SNR))
     * This is the fundamental limit of ranging precision */
    return SPEED_OF_LIGHT_MPS / (2.0 * sqrt(2.0) * M_PI * bandwidth_hz * sqrt(snr_linear));
}

/* ============================================================================
 * L6 - First-Path Detection in UWB Channel Impulse Response
 * ============================================================================ */

double uwb_first_path_tof(const double *cir_magnitude,
                          double cir_sample_period_s,
                          int n_samples,
                          double noise_threshold) {
    if (!cir_magnitude || n_samples <= 0 || cir_sample_period_s <= 0.0) return -1.0;

    /* Estimate noise floor from tail of CIR (last 20% of samples) */
    int noise_start = n_samples - n_samples / 5;
    if (noise_start < n_samples / 2) noise_start = n_samples / 2;

    double noise_sum = 0.0;
    int noise_count = 0;
    for (int i = noise_start; i < n_samples; i++) {
        noise_sum += cir_magnitude[i] * cir_magnitude[i];
        noise_count++;
    }
    double noise_floor = (noise_count > 0) ? sqrt(noise_sum / noise_count) : 0.0;

    /* Detection threshold: noise_floor + margin */
    double threshold = noise_floor * noise_threshold;
    if (threshold < 1e-12) threshold = 1e-12;

    /* Find first sample exceeding the detection threshold */
    /* Leading edge detection: look for rising slope and threshold crossing */
    for (int i = 1; i < n_samples; i++) {
        if (cir_magnitude[i] > threshold && cir_magnitude[i] > cir_magnitude[i-1]) {
            /* Found first path. Use linear interpolation for sub-sample accuracy */
            double slope = cir_magnitude[i] - cir_magnitude[i-1];
            double offset = 0.0;
            if (slope > 1e-12) {
                offset = (threshold - cir_magnitude[i-1]) / slope;
            }
            return (i - 1.0 + offset) * cir_sample_period_s;
        }
    }

    return -1.0;  /* No first path detected */
}
