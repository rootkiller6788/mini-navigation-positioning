#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880
#endif
/**
 * mini-uwb-localization: UWB Core Types Implementation
 *
 * Initialization functions and utility computations for UWB
 * localization data structures.
 *
 * Knowledge Coverage: L1 Definitions �� UWB system types
 *                      L4 Fundamental Laws �� CRLB, GDOP
 *                      L3 Math Structures �� Distance geometry
 */

#include "uwb_types.h"
#include <math.h>
#include <string.h>

/* =========================================================================
 * Initialization Functions
 * ========================================================================= */

void uwb_anchor_init(uwb_anchor_t *anchor, uint16_t id, double x, double y, double z)
{
    if (!anchor) return;
    memset(anchor, 0, sizeof(*anchor));
    anchor->id = id;
    anchor->position.x = x;
    anchor->position.y = y;
    anchor->position.z = z;
    anchor->channel = UWB_CHANNEL_5;
    anchor->prf = UWB_PRF_64MHZ;
    anchor->tx_power_dbm = -14.0;
    anchor->antenna_gain_dbi = 2.0;
    anchor->clock_offset_ppm = 0.0;
    anchor->is_active = 1;
    anchor->is_synchronized = 1;
    anchor->last_ranging_ts = 0;
}

void uwb_tag_init(uwb_tag_t *tag, uint16_t id)
{
    if (!tag) return;
    memset(tag, 0, sizeof(*tag));
    tag->id = id;
    tag->position.x = 0.0;
    tag->position.y = 0.0;
    tag->position.z = 0.0;
    tag->velocity.x = 0.0;
    tag->velocity.y = 0.0;
    tag->velocity.z = 0.0;
    tag->pos_cov.var_x = 1.0;
    tag->pos_cov.var_y = 1.0;
    tag->pos_cov.var_z = 1.0;
    tag->pos_cov.cov_xy = 0.0;
    tag->pos_cov.cov_xz = 0.0;
    tag->pos_cov.cov_yz = 0.0;
    tag->clock_offset_ppm = 0.0;
    tag->battery_voltage = 3.3;
    tag->motion_state = 0;
    tag->last_update_ts = 0;
}

void uwb_config_init_default(uwb_config_t *config)
{
    if (!config) return;
    memset(config, 0, sizeof(*config));
    config->channel = UWB_CHANNEL_5;
    config->prf = UWB_PRF_64MHZ;
    config->preamble_len = UWB_PREAMBLE_1024;
    config->datarate = UWB_DATARATE_6M8;
    config->preamble_code = 9;
    config->sfd_id = 1;
    config->tx_power_dbm = -14.0;
    config->antenna_delay_tx_ps = 16450.0;
    config->antenna_delay_rx_ps = 16450.0;
    config->smart_power_enabled = 0;
    config->leading_edge_detection = 1;
    config->ranging_interval_ms = 100;
    config->slot_duration_us = 2000;
    config->response_timeout_us = 20000;
}

void uwb_ranging_meas_init(uwb_ranging_meas_t *meas)
{
    if (!meas) return;
    memset(meas, 0, sizeof(*meas));
    meas->type = UWB_RANGING_SS_TWR;
    meas->quality = UWB_RANGE_QUALITY_POOR;
    meas->distance_variance = 1.0;
}

void uwb_error_metrics_init(uwb_error_metrics_t *metrics)
{
    if (!metrics) return;
    memset(metrics, 0, sizeof(*metrics));
    metrics->min_error = 1e100;
    metrics->crlb = 0.0;
    metrics->gdop_avg = 0.0;
}

/* =========================================================================
 * Distance Geometry
 * ========================================================================= */

double uwb_distance_2d(const uwb_pos2d_t *a, const uwb_pos2d_t *b)
{
    double dx, dy;
    if (!a || !b) return -1.0;
    dx = a->x - b->x;
    dy = a->y - b->y;
    return sqrt(dx * dx + dy * dy);
}

double uwb_distance_3d(const uwb_pos3d_t *a, const uwb_pos3d_t *b)
{
    double dx, dy, dz;
    if (!a || !b) return -1.0;
    dx = a->x - b->x;
    dy = a->y - b->y;
    dz = a->z - b->z;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

/* =========================================================================
 * Channel Information
 * ========================================================================= */

double uwb_channel_frequency_hz(uwb_channel_t channel)
{
    switch (channel) {
    case UWB_CHANNEL_1:  return 3494.4e6;
    case UWB_CHANNEL_2:  return 3993.6e6;
    case UWB_CHANNEL_3:  return 4492.8e6;
    case UWB_CHANNEL_4:  return 3993.6e6;
    case UWB_CHANNEL_5:  return 6489.6e6;
    case UWB_CHANNEL_7:  return 6489.6e6;
    case UWB_CHANNEL_9:  return 7987.2e6;
    case UWB_CHANNEL_11: return 7987.2e6;
    default:             return 0.0;
    }
}

double uwb_channel_bandwidth_hz(uwb_channel_t channel)
{
    switch (channel) {
    case UWB_CHANNEL_1:  return 499.2e6;
    case UWB_CHANNEL_2:  return 499.2e6;
    case UWB_CHANNEL_3:  return 499.2e6;
    case UWB_CHANNEL_4:  return 1331.2e6;
    case UWB_CHANNEL_5:  return 499.2e6;
    case UWB_CHANNEL_7:  return 1081.6e6;
    case UWB_CHANNEL_9:  return 499.2e6;
    case UWB_CHANNEL_11: return 1331.2e6;
    default:             return 0.0;
    }
}

/* =========================================================================
 * Cramer-Rao Lower Bound for TOA-based Ranging
 *
 * CRLB = c / (2 * sqrt(2) * pi * SNR * B_eff)
 *
 * Derivation (L4 Fundamental Law):
 * The time-delay estimation for a signal s(t) with AWGN has:
 *   var(tau_hat) >= 1 / (SNR * B_rms^2)
 * where B_rms = 2*pi * sqrt(integral(f^2|S(f)|^2 df) / integral(|S(f)|^2 df))
 *
 * For a signal with flat spectrum over B_eff:
 *   B_rms^2 = (pi^2 / 3) * B_eff^2
 *   var(tau) >= 3 / (pi^2 * SNR * B_eff^2)
 *
 * The standard CRLB form for UWB:
 *   sigma_d >= c / (2 * sqrt(2) * pi * SNR * B_eff)
 *
 * Reference: Gezici et al. (2005) "Localization via Ultra-Wideband Radios"
 *             IEEE Signal Processing Magazine
 *
 * Complexity: O(1)
 * @param snr_linear    SNR (linear scale, not dB)
 * @param bw_effective  Effective signal bandwidth [Hz]
 * @return              CRLB for distance estimation [m]
 */
double uwb_crlb_distance(double snr_linear, double bw_effective)
{
    double denominator;

    if (snr_linear <= 0.0 || bw_effective <= 0.0) {
        return 1e100; /* degenerate case: infinite bound */
    }

    /* CRLB = c / (2*sqrt(2)*pi * SNR * B_eff) */
    denominator = 2.0 * M_SQRT2 * M_PI * snr_linear * bw_effective;
    if (denominator < 1e-20) return 1e100;

    return UWB_C / denominator;
}

/* =========================================================================
 * Geometric Dilution of Precision (GDOP)
 *
 * GDOP = sqrt(trace((H^T H)^{-1}))
 *
 * where H is the N x D design matrix of unit vectors from tag to each
 * anchor. Each row of H is:
 *   h_i = [(x - a_ix)/r_i, (y - a_iy)/r_i, (z - a_iz)/r_i]
 *
 * GDOP quantifies how anchor geometry amplifies ranging errors:
 *   sigma_position = GDOP * sigma_range
 *
 * Good geometry: GDOP < 2   (anchors well-distributed)
 * Fair geometry: 2 <= GDOP < 5
 * Poor geometry: GDOP >= 5 (anchors nearly coplanar/collinear)
 *
 * Reference: Kaplan & Hegarty (2017) "Understanding GPS/GNSS Principles"
 *             Chapter 7: GPS Performance and Error Effects
 *
 * Complexity: O(N * D^2 + D^3)
 */
double uwb_compute_gdop(const uwb_pos3d_t *anchors, int num_anchors,
                        const uwb_pos3d_t *tag_pos)
{
    double H[12]; (void)H;     /* max 4 anchors x 3 dims = 12 */
    double HTH[9];    /* 3 x 3 */
    double HTH_inv[9];
    double range;
    int i, j, k;
    (void)j; (void)k;
    double trace_val;
    double det;
    int dim = 3; (void)dim;

    if (!anchors || !tag_pos || num_anchors < 4) {
        return 1e100;
    }

    /* Build H^T * H as a 3x3 matrix */
    memset(HTH, 0, 9 * sizeof(double));

    for (i = 0; i < num_anchors; i++) {
        double dx = tag_pos->x - anchors[i].x;
        double dy = tag_pos->y - anchors[i].y;
        double dz = tag_pos->z - anchors[i].z;
        range = sqrt(dx * dx + dy * dy + dz * dz);

        if (range < 1e-6) {
            /* Tag coincident with anchor �� avoid division by zero */
            continue;
        }

        /* Unit vector components */
        double ux = dx / range;
        double uy = dy / range;
        double uz = dz / range;

        /* Rank-1 update: HTH += u * u^T */
        HTH[0] += ux * ux;  HTH[1] += ux * uy;  HTH[2] += ux * uz;
        HTH[3] += uy * ux;  HTH[4] += uy * uy;  HTH[5] += uy * uz;
        HTH[6] += uz * ux;  HTH[7] += uz * uy;  HTH[8] += uz * uz;
    }

    /* Invert 3x3 symmetric matrix HTH
     * Using adjugate method: A^{-1} = adj(A) / det(A) */
    det = HTH[0] * (HTH[4] * HTH[8] - HTH[5] * HTH[7])
        - HTH[1] * (HTH[3] * HTH[8] - HTH[5] * HTH[6])
        + HTH[2] * (HTH[3] * HTH[7] - HTH[4] * HTH[6]);

    if (fabs(det) < 1e-15) {
        return 1e100; /* Singular �� degenerate anchor geometry */
    }

    HTH_inv[0] =  (HTH[4] * HTH[8] - HTH[5] * HTH[7]) / det;
    HTH_inv[1] = -(HTH[1] * HTH[8] - HTH[2] * HTH[7]) / det;
    HTH_inv[2] =  (HTH[1] * HTH[5] - HTH[2] * HTH[4]) / det;
    HTH_inv[3] = -(HTH[3] * HTH[8] - HTH[5] * HTH[6]) / det;
    HTH_inv[4] =  (HTH[0] * HTH[8] - HTH[2] * HTH[6]) / det;
    HTH_inv[5] = -(HTH[0] * HTH[5] - HTH[2] * HTH[3]) / det;
    HTH_inv[6] =  (HTH[3] * HTH[7] - HTH[4] * HTH[6]) / det;
    HTH_inv[7] = -(HTH[0] * HTH[7] - HTH[1] * HTH[6]) / det;
    HTH_inv[8] =  (HTH[0] * HTH[4] - HTH[1] * HTH[3]) / det;

    /* GDOP = sqrt(trace of inverse) */
    trace_val = HTH_inv[0] + HTH_inv[4] + HTH_inv[8];
    if (trace_val < 0.0) trace_val = 0.0;

    return sqrt(trace_val);
}
