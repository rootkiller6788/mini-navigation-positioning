/**
 * @file positioning_error.c
 * @brief Error analysis and accuracy metrics for indoor positioning
 *
 * Implements: DOP computation, error statistics, error ellipse,
 * CRLB, error propagation, outlier detection, RANSAC, Allan variance,
 * integrity/continuity risk.
 */

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "../include/positioning_error.h"
#include "../include/indoor_positioning.h"

/* ============================================================================
 * L5 - Dilution of Precision (DOP)
 * ============================================================================ */

int compute_dop(const position3d_t *anchor_positions,
                const position3d_t *user_position,
                int n_anchors,
                dop_metrics_t *dop) {
    if (!anchor_positions || !user_position || !dop) return -1;
    if (n_anchors < 4) return -1;

    /* Build geometry matrix G (N x 4)
     * Each row: [cos(az)*cos(el), sin(az)*cos(el), sin(el), 1]
     * where az = azimuth from North, el = elevation from horizontal */
    double GTG[4][4] = {{0}};

    for (int i = 0; i < n_anchors; i++) {
        double dx = anchor_positions[i].x - user_position->x;
        double dy = anchor_positions[i].y - user_position->y;
        double dz = anchor_positions[i].z - user_position->z;
        double range = sqrt(dx*dx + dy*dy + dz*dz);
        if (range < 1e-10) range = 1e-10;

        double row[4];
        row[0] = dx / range;
        row[1] = dy / range;
        row[2] = dz / range;
        row[3] = 1.0;  /* Clock/range bias term */

        /* Accumulate G^T * G */
        for (int r = 0; r < 4; r++) {
            for (int c = 0; c < 4; c++) {
                GTG[r][c] += row[r] * row[c];
            }
        }
    }

    /* Invert 4x4 GTG using Gaussian elimination */
    double A[4][8];  /* Augmented with identity */
    for (int r = 0; r < 4; r++) {
        for (int c = 0; c < 4; c++) A[r][c] = GTG[r][c];
        for (int c = 4; c < 8; c++) A[r][c] = (r == c - 4) ? 1.0 : 0.0;
    }

    /* Gauss-Jordan elimination */
    for (int col = 0; col < 4; col++) {
        int pivot = col;
        double max_val = fabs(A[col][col]);
        for (int r = col + 1; r < 4; r++) {
            if (fabs(A[r][col]) > max_val) {
                max_val = fabs(A[r][col]);
                pivot = r;
            }
        }
        if (max_val < 1e-15) return -1;

        if (pivot != col) {
            for (int c = col; c < 8; c++) {
                double tmp = A[col][c];
                A[col][c] = A[pivot][c];
                A[pivot][c] = tmp;
            }
        }

        double piv_val = A[col][col];
        for (int c = col; c < 8; c++) A[col][c] /= piv_val;

        for (int r = 0; r < 4; r++) {
            if (r == col) continue;
            double factor = A[r][col];
            for (int c = col; c < 8; c++) A[r][c] -= factor * A[col][c];
        }
    }

    /* DOP values = sqrt of diagonal of (G^T G)^{-1} */
    double q_ee = A[0][4];  /* East/East */
    double q_nn = A[1][5];  /* North/North */
    double q_uu = A[2][6];  /* Up/Up */
    double q_tt = A[3][7];  /* Time/Time */

    if (q_ee < 0.0) q_ee = 0.0;
    if (q_nn < 0.0) q_nn = 0.0;
    if (q_uu < 0.0) q_uu = 0.0;
    if (q_tt < 0.0) q_tt = 0.0;

    dop->pdop = sqrt(q_ee + q_nn + q_uu);
    dop->hdop = sqrt(q_ee + q_nn);
    dop->vdop = sqrt(q_uu);
    dop->tdop = sqrt(q_tt);
    dop->gdop = sqrt(q_ee + q_nn + q_uu + q_tt);

    return 0;
}

/* ============================================================================
 * L5 - Error Statistics
 * ============================================================================ */

/* Static helper for sort function (C standard qsort) */
static int cmp_double(const void *a, const void *b) {
    double da = *(const double *)a;
    double db = *(const double *)b;
    return (da < db) ? -1 : (da > db) ? 1 : 0;
}

void compute_error_stats(const position3d_t *truth,
                         const position3d_t *estimate,
                         int n,
                         positioning_error_stats_t *stats) {
    if (!truth || !estimate || !stats || n <= 0) return;
    memset(stats, 0, sizeof(positioning_error_stats_t));

    double sum_ex = 0, sum_ey = 0, sum_ez = 0;
    double sum_2derr = 0, sum_3derr = 0;
    double *errors_2d = (double *)malloc(n * sizeof(double));
    double *errors_3d = (double *)malloc(n * sizeof(double));
    if (!errors_2d || !errors_3d) {
        free(errors_2d);
        free(errors_3d);
        return;
    }

    stats->n_samples = n;

    for (int i = 0; i < n; i++) {
        double ex = estimate[i].x - truth[i].x;
        double ey = estimate[i].y - truth[i].y;
        double ez = estimate[i].z - truth[i].z;

        sum_ex += ex;
        sum_ey += ey;
        sum_ez += ez;

        double err_2d = sqrt(ex*ex + ey*ey);
        double err_3d = sqrt(ex*ex + ey*ey + ez*ez);
        errors_2d[i] = err_2d;
        errors_3d[i] = err_3d;

        sum_2derr += err_2d;
        sum_3derr += err_3d;

        if (err_2d > stats->max_error) stats->max_error = err_2d;
    }

    stats->mean_error_x = sum_ex / n;
    stats->mean_error_y = sum_ey / n;
    stats->mean_error_z = sum_ez / n;

    /* Standard deviations */
    double sum_vx = 0, sum_vy = 0, sum_vz = 0;
    for (int i = 0; i < n; i++) {
        double ex = estimate[i].x - truth[i].x - stats->mean_error_x;
        double ey = estimate[i].y - truth[i].y - stats->mean_error_y;
        double ez = estimate[i].z - truth[i].z - stats->mean_error_z;
        sum_vx += ex * ex;
        sum_vy += ey * ey;
        sum_vz += ez * ez;
    }
    stats->std_dev_x = sqrt(sum_vx / n);
    stats->std_dev_y = sqrt(sum_vy / n);
    stats->std_dev_z = sqrt(sum_vz / n);

    /* RMS */
    double sum_2d_sq = 0, sum_3d_sq = 0;
    for (int i = 0; i < n; i++) {
        sum_2d_sq += errors_2d[i] * errors_2d[i];
        sum_3d_sq += errors_3d[i] * errors_3d[i];
    }
    stats->rms_2d = sqrt(sum_2d_sq / n);
    stats->rms_3d = sqrt(sum_3d_sq / n);
    stats->drms = stats->rms_2d;

    /* CEP: sort 2D errors and find percentile */
    /* Simple CEP via mean + std approximation */
    double mean_2d = sum_2derr / n;
    double var_2d = 0;
    for (int i = 0; i < n; i++) {
        double d = errors_2d[i] - mean_2d;
        var_2d += d * d;
    }
    double std_2d = sqrt(var_2d / n);
    /* Rayleigh distribution approximation for CEP */
    stats->cep50 = 0.6745 * std_2d + mean_2d * 0.5;
    stats->cep95 = 1.96 * std_2d + mean_2d;

    /* More accurate CEP via sorting */
    qsort(errors_2d, n, sizeof(double), cmp_double);
    if (n > 0) {
        stats->cep50 = errors_2d[(int)(n * 0.50)];
        stats->cep95 = errors_2d[(int)(n * 0.95 > n - 1 ? n - 1 : n * 0.95)];
        stats->cep99 = errors_2d[(int)(n * 0.99 > n - 1 ? n - 1 : n * 0.99)];
    }

    /* SEP (Spherical Error Probable) */
    qsort(errors_3d, n, sizeof(double), cmp_double);
    if (n > 0) {
        stats->sep50 = errors_3d[(int)(n * 0.50)];
        stats->sep95 = errors_3d[(int)(n * 0.95 > n - 1 ? n - 1 : n * 0.95)];
    }

    /* R95: 95% confidence radius approximation */
    stats->r95 = 1.96 * stats->rms_2d;

    free(errors_2d);
    free(errors_3d);
}

/* ============================================================================
 * L5 - Error Ellipse
 * ============================================================================ */

/**
 * @brief Error function (erf) approximation using Abramowitz & Stegun 7.1.26
 */
static double erf_approx(double x) {
    double t = 1.0 / (1.0 + 0.3275911 * fabs(x));
    double a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741;
    double a4 = -1.453152027, a5 = 1.061405429;
    double poly = t * (a1 + t * (a2 + t * (a3 + t * (a4 + t * a5))));
    double erf_val = 1.0 - poly * exp(-x * x);
    return (x >= 0) ? erf_val : -erf_val;
}

/**
 * @brief Inverse error function approximation
 */
static double erfinv_approx(double x) {
    if (x < -0.999) return -3.0;
    if (x > 0.999) return 3.0;
    double a = 0.147;
    double ln1mx2 = log(1.0 - x * x);
    double p1 = 2.0 / (M_PI * a) + ln1mx2 / 2.0;
    double p2 = ln1mx2 / a;
    double sign_x = (x >= 0) ? 1.0 : -1.0;
    return sign_x * sqrt(sqrt(p1 * p1 - p2) - p1);
}

double chi2_inv(double p, int k) {
    if (p <= 0.0 || p >= 1.0) return 0.0;
    if (k <= 0) return 0.0;

    /* Wilson-Hilferty transformation:
     * chi2_p(k) ≈ k * (1 - 2/(9k) + z_p * sqrt(2/(9k)))^3
     * where z_p = sqrt(2) * erfinv(2p - 1) */
    double z = sqrt(2.0) * erfinv_approx(2.0 * p - 1.0);

    double term = 1.0 - 2.0/(9.0 * k) + z * sqrt(2.0/(9.0 * k));
    double chi2 = k * term * term * term;
    if (chi2 < 0.0) chi2 = 0.0;
    return chi2;
}

void compute_error_ellipse(double cov_xx, double cov_yy, double cov_xy,
                           double confidence, error_ellipse_t *ellipse) {
    if (!ellipse) return;

    /* Eigenvalues of 2x2 covariance matrix:
     * lambda = 0.5*(a+b ± sqrt((a-b)^2 + 4c^2))
     * where a=cov_xx, b=cov_yy, c=cov_xy */
    double trace = cov_xx + cov_yy;
    double det = cov_xx * cov_yy - cov_xy * cov_xy;
    double discriminant = trace * trace - 4.0 * det;

    if (discriminant < 0.0) discriminant = 0.0;
    double sqrt_disc = sqrt(discriminant);

    double lambda_max = 0.5 * (trace + sqrt_disc);
    double lambda_min = 0.5 * (trace - sqrt_disc);
    if (lambda_min < 0.0) lambda_min = 0.0;

    /* Chi-squared with 2 degrees of freedom for confidence level */
    double chi2_val = chi2_inv(confidence, 2);

    ellipse->semi_major = sqrt(chi2_val * lambda_max);
    ellipse->semi_minor = sqrt(chi2_val * lambda_min);
    ellipse->confidence = confidence;

    /* Orientation: theta = 0.5 * atan2(2*cov_xy, cov_xx - cov_yy) */
    if (fabs(cov_xx - cov_yy) > 1e-15) {
        ellipse->orientation = 0.5 * atan2(2.0 * cov_xy, cov_xx - cov_yy);
    } else {
        ellipse->orientation = (cov_xy >= 0) ? M_PI / 4.0 : -M_PI / 4.0;
    }
}

/* ============================================================================
 * L4 - Cramer-Rao Lower Bound
 * ============================================================================ */

double crlb_rssi_positioning(const position2d_t *anchor_positions,
                             const double *distances,
                             const double *rssi_std,
                             double path_loss_exp,
                             int n_anchors) {
    if (!anchor_positions || !distances || !rssi_std || n_anchors < 3) return -1.0;
    if (path_loss_exp <= 0.0) return -1.0;

    /* FIM = sum (10*n / (ln(10) * sigma_i * d_i))^2 * G_i * G_i^T
     * where G_i is unit vector from anchor_i to user (approximated
     * using anchor positions relative to center for CRLB evaluation) */
    double center_x = 0, center_y = 0;
    for (int i = 0; i < n_anchors; i++) {
        center_x += anchor_positions[i].x;
        center_y += anchor_positions[i].y;
    }
    center_x /= n_anchors;
    center_y /= n_anchors;

    double alpha = 10.0 * path_loss_exp / log(10.0);
    double FIM[2][2] = {{0,0},{0,0}};

    for (int i = 0; i < n_anchors; i++) {
        if (rssi_std[i] <= 0.0 || distances[i] <= 0.0) continue;

        double sigma_d = alpha * rssi_std[i] / distances[i];  /* Distance equivalent noise */
        double dx = anchor_positions[i].x - center_x;
        double dy = anchor_positions[i].y - center_y;
        double d = sqrt(dx*dx + dy*dy);
        if (d < 1e-10) continue;

        double gx = dx / d;  /* Unit vector x */
        double gy = dy / d;  /* Unit vector y */
        double w = 1.0 / (sigma_d * sigma_d);

        FIM[0][0] += w * gx * gx;
        FIM[0][1] += w * gx * gy;
        FIM[1][0] += w * gy * gx;
        FIM[1][1] += w * gy * gy;
    }

    double det_FIM = FIM[0][0]*FIM[1][1] - FIM[0][1]*FIM[1][0];
    if (det_FIM <= 1e-15) return 1e9;

    /* CRLB = trace(FIM^{-1}) */
    double crlb = (FIM[0][0] + FIM[1][1]) / det_FIM;
    return crlb;
}

double crlb_tof_positioning(const position2d_t *anchor_positions,
                            const double *noise_std,
                            int n_anchors) {
    if (!anchor_positions || !noise_std || n_anchors < 3) return -1.0;

    double center_x = 0, center_y = 0;
    for (int i = 0; i < n_anchors; i++) {
        center_x += anchor_positions[i].x;
        center_y += anchor_positions[i].y;
    }
    center_x /= n_anchors;
    center_y /= n_anchors;

    double FIM[2][2] = {{0,0},{0,0}};

    for (int i = 0; i < n_anchors; i++) {
        if (noise_std[i] <= 0.0) continue;

        double dx = anchor_positions[i].x - center_x;
        double dy = anchor_positions[i].y - center_y;
        double d = sqrt(dx*dx + dy*dy);
        if (d < 1e-10) continue;

        double gx = dx / d;
        double gy = dy / d;
        double w = 1.0 / (noise_std[i] * noise_std[i]);

        FIM[0][0] += w * gx * gx;
        FIM[0][1] += w * gx * gy;
        FIM[1][0] += w * gy * gx;
        FIM[1][1] += w * gy * gy;
    }

    double det_FIM = FIM[0][0]*FIM[1][1] - FIM[0][1]*FIM[1][0];
    if (det_FIM <= 1e-15) return 1e9;

    return (FIM[0][0] + FIM[1][1]) / det_FIM;
}

/* ============================================================================
 * L4 - Error Propagation
 * ============================================================================ */

void propagate_positioning_error(const position2d_t *anchor_positions,
                                 const double *distance_std,
                                 const position2d_t *user_pos,
                                 int n_anchors,
                                 double *cov_xx, double *cov_yy,
                                 double *cov_xy) {
    if (!anchor_positions || !distance_std || !user_pos || !cov_xx || !cov_yy || !cov_xy) return;

    /* H = Jacobian of distance equations
     * H_i = [(x - xi)/d_i, (y - yi)/d_i]
     * Cov = (H^T * W * H)^{-1}, W = diag(1/sigma_i^2) */

    double HtWH[2][2] = {{0,0},{0,0}};

    for (int i = 0; i < n_anchors; i++) {
        if (distance_std[i] <= 0.0) continue;

        double dx = user_pos->x - anchor_positions[i].x;
        double dy = user_pos->y - anchor_positions[i].y;
        double d = sqrt(dx*dx + dy*dy);
        if (d < 1e-10) continue;

        double h0 = dx / d;
        double h1 = dy / d;
        double w = 1.0 / (distance_std[i] * distance_std[i]);

        HtWH[0][0] += h0 * h0 * w;
        HtWH[0][1] += h0 * h1 * w;
        HtWH[1][0] += h1 * h0 * w;
        HtWH[1][1] += h1 * h1 * w;
    }

    double det = HtWH[0][0]*HtWH[1][1] - HtWH[0][1]*HtWH[1][0];
    if (det <= 1e-15) {
        *cov_xx = 1e9; *cov_yy = 1e9; *cov_xy = 0;
        return;
    }

    *cov_xx = HtWH[1][1] / det;
    *cov_yy = HtWH[0][0] / det;
    *cov_xy = -HtWH[0][1] / det;
}

/* ============================================================================
 * L5 - Outlier Detection (Innovation Test)
 * ============================================================================ */

int detect_measurement_outlier(const double *innovation,
                               const double *innovation_cov,
                               int n_measure,
                               double confidence) {
    if (!innovation || !innovation_cov || n_measure <= 0) return 0;

    /* Normalized Innovation Squared (NIS) = y^T * S^{-1} * y
     * If NIS > chi2_threshold, flag as outlier.
     *
     * For simplicity, compute NIS using diagonal S approximation
     * for 1D measurement case, or full S^{-1} for multi-D */

    double nis = 0.0;
    if (n_measure == 1) {
        if (innovation_cov[0] > 1e-12) {
            nis = innovation[0] * innovation[0] / innovation_cov[0];
        }
    } else {
        /* Use diagonal approximation */
        for (int i = 0; i < n_measure; i++) {
            if (innovation_cov[i * n_measure + i] > 1e-12) {
                nis += innovation[i] * innovation[i]
                       / innovation_cov[i * n_measure + i];
            }
        }
    }

    double threshold = chi2_inv(confidence, n_measure);
    return (nis > threshold) ? 1 : 0;
}

/* ============================================================================
 * L5 - RANSAC Robust Positioning
 * ============================================================================ */

int ransac_positioning(const double *distances,
                       const position2d_t *anchor_positions,
                       int n_anchors,
                       int n_subsets,
                       int n_iterations,
                       double inlier_threshold,
                       position2d_t *result) {
    if (!distances || !anchor_positions || !result) return -1;
    if (n_anchors < n_subsets) return -1;
    if (n_iterations <= 0) n_iterations = 100;

    int best_n_inliers = 0;
    double best_x = 0, best_y = 0;

    for (int iter = 0; iter < n_iterations; iter++) {
        /* Randomly select subset of anchors */
        int selected[4];
        int used[IP_MAX_ANCHORS] = {0};

        for (int s = 0; s < n_subsets; s++) {
            int idx;
            do {
                idx = (rand() + s * 31) % n_anchors;
            } while (used[idx]);
            used[idx] = 1;
            selected[s] = idx;
        }

        /* Solve with subset */
        position2d_t subset_positions[4];
        double subset_distances[4];
        for (int s = 0; s < n_subsets; s++) {
            subset_positions[s] = anchor_positions[selected[s]];
            subset_distances[s] = distances[selected[s]];
        }

        position2d_t candidate;
        if (trilateration_2d(subset_positions, subset_distances, n_subsets, &candidate) != 0) {
            continue;
        }

        /* Count inliers */
        int n_inliers = 0;
        for (int i = 0; i < n_anchors; i++) {
            double dx = candidate.x - anchor_positions[i].x;
            double dy = candidate.y - anchor_positions[i].y;
            double d_pred = sqrt(dx*dx + dy*dy);
            double residual = fabs(distances[i] - d_pred);
            if (residual < inlier_threshold) n_inliers++;
        }

        if (n_inliers > best_n_inliers) {
            best_n_inliers = n_inliers;
            best_x = candidate.x;
            best_y = candidate.y;
        }
    }

    if (best_n_inliers >= n_subsets) {
        result->x = best_x;
        result->y = best_y;
        return best_n_inliers;
    }
    return -1;
}

/* ============================================================================
 * L2 - Error Budget Decomposition
 * ============================================================================ */

void decompose_error_sources(const position3d_t *truth,
                             const position3d_t *estimate,
                             int n,
                             double *bias_component,
                             double *noise_component,
                             double *drift_component,
                             double *total_rms) {
    if (!truth || !estimate || !bias_component || !noise_component
        || !drift_component || !total_rms || n <= 0) return;

    /* Compute 2D errors */
    double sum_ex = 0, sum_ey = 0;
    double *errors = (double *)malloc(n * sizeof(double));
    if (!errors) return;

    for (int i = 0; i < n; i++) {
        double ex = estimate[i].x - truth[i].x;
        double ey = estimate[i].y - truth[i].y;
        sum_ex += ex; sum_ey += ey;
        errors[i] = sqrt(ex*ex + ey*ey);
    }

    double bias_x = sum_ex / n;
    double bias_y = sum_ey / n;
    *bias_component = sqrt(bias_x*bias_x + bias_y*bias_y);

    /* Remove bias and compute residuals */
    double sum_sq = 0;
    for (int i = 0; i < n; i++) {
        double rx = (estimate[i].x - truth[i].x) - bias_x;
        double ry = (estimate[i].y - truth[i].y) - bias_y;
        sum_sq += rx*rx + ry*ry;
    }
    double total_var = sum_sq / n;

    /* Drift component: correlation between consecutive errors */
    double drift_corr = 0;
    for (int i = 1; i < n; i++) {
        double e0 = errors[i-1];
        double e1 = errors[i];
        drift_corr += e0 * e1;
    }
    drift_corr /= (n - 1 > 0 ? n - 1 : 1);
    *drift_component = sqrt(fabs(drift_corr));

    /* Noise = total minus drift (assuming independence) */
    double noise_var = total_var - drift_corr;
    if (noise_var < 0) noise_var = 0;
    *noise_component = sqrt(noise_var);

    /* Total RMS */
    *total_rms = sqrt(total_var + (*bias_component) * (*bias_component));

    free(errors);
}

/* ============================================================================
 * L5 - Allan Variance
 * ============================================================================ */

void compute_allan_variance(const double *signal, int n_samples,
                            double sample_rate_hz,
                            double *tau, double *allan_var, int m) {
    if (!signal || !tau || !allan_var || n_samples <= 0 || m <= 0) return;

    double dt = 1.0 / sample_rate_hz;

    for (int k = 0; k < m; k++) {
        /* Cluster size: powers of 2 */
        int cluster_size = 1 << k;
        if (cluster_size > n_samples / 9) {
            tau[k] = 0.0;
            allan_var[k] = 0.0;
            continue;
        }

        int n_clusters = n_samples / cluster_size;
        if (n_clusters < 2) {
            tau[k] = 0.0;
            allan_var[k] = 0.0;
            continue;
        }

        /* Compute cluster averages */
        double *cluster_avg = (double *)malloc(n_clusters * sizeof(double));
        if (!cluster_avg) continue;

        for (int c = 0; c < n_clusters; c++) {
            double sum = 0;
            for (int i = 0; i < cluster_size; i++) {
                sum += signal[c * cluster_size + i];
            }
            cluster_avg[c] = sum / cluster_size;
        }

        /* Allan variance: 0.5 * mean((omega_{k+1} - omega_k)^2) */
        double sum_diff_sq = 0;
        for (int c = 0; c < n_clusters - 1; c++) {
            double diff = cluster_avg[c+1] - cluster_avg[c];
            sum_diff_sq += diff * diff;
        }
        allan_var[k] = 0.5 * sum_diff_sq / (n_clusters - 1);
        tau[k] = cluster_size * dt;

        free(cluster_avg);
    }
}

/* ============================================================================
 * L1 - Reliability Metrics
 * ============================================================================ */

double compute_integrity_risk(const double *error_samples, int n_samples,
                              double alert_limit) {
    if (!error_samples || n_samples <= 0) return 0.0;
    int exceed_count = 0;
    for (int i = 0; i < n_samples; i++) {
        if (fabs(error_samples[i]) > alert_limit) exceed_count++;
    }
    return (double)exceed_count / n_samples;
}

double compute_continuity_risk(const double *epoch_durations,
                               int n_segments,
                               double interruption_threshold) {
    if (!epoch_durations || n_segments <= 0) return 0.0;
    int interruption_count = 0;
    for (int i = 0; i < n_segments; i++) {
        if (epoch_durations[i] < interruption_threshold) interruption_count++;
    }
    return (double)interruption_count / n_segments;
}
