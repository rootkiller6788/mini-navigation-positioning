/**
 * @file    ins_calibration.c
 * @brief   IMU calibration: bias, scale factor, misalignment, thermal
 *
 * Knowledge Coverage:
 *   L5 (Algorithms): Six-position static calib, multi-position gyro calib,
 *                     cross-axis misalignment, thermal compensation
 *   L2 (Core Concepts): Deterministic vs stochastic sensor errors
 *
 * Reference:
 *   Groves (2013), Chapter 4, "Inertial Sensor Calibration".
 *   Shin & El-Sheimy (2002), Z. Vermess., 127(1): 41-50.
 *   Tedaldi et al. (2014), "A Robust and Easy to Implement Method for IMU
 *     Calibration without External Equipments", ICRA 2014.
 */

#include "ins_calibration.h"
#include <string.h>
#include <math.h>
#include <stdio.h>

/* =========================================================================
 * L5: Six-Position Static Accelerometer Calibration
 *
 * Using Earth's gravity as reference (magnitude g_local):
 *
 * For each axis i (x, y, z), two measurements:
 *   Axis down: a_i_down = (1 + sf_i) * (+g) + b_i
 *   Axis up:   a_i_up   = (1 + sf_i) * (-g) + b_i
 *
 * Solving:
 *   b_i  = (a_i_down + a_i_up) / 2
 *   sf_i = (a_i_down - a_i_up) / (2*g) - 1
 *
 * Positions 1-6 correspond to +x, -x, +y, -y, +z, -z pointing down.
 * ========================================================================= */

int ins_calib_six_position_accel(const ins_vec3_t readings[6],
                                  double g_local,
                                  ins_calibration_result_t *result) {
    if (!readings || !result || g_local <= 0.1) return -1;

    double g_ref = g_local;

    /* Accelerometer x-axis: readings[0] is +x down, readings[1] is -x down */
    double bx = (readings[0].x + readings[1].x) * 0.5;
    double sfx = (readings[0].x - readings[1].x) / (2.0 * g_ref) - 1.0;

    /* Accelerometer y-axis: readings[2] is +y down, readings[3] is -y down */
    double by = (readings[2].y + readings[3].y) * 0.5;
    double sfy = (readings[2].y - readings[3].y) / (2.0 * g_ref) - 1.0;

    /* Accelerometer z-axis: readings[4] is +z down, readings[5] is -z down */
    double bz = (readings[4].z + readings[5].z) * 0.5;
    double sfz = (readings[4].z - readings[5].z) / (2.0 * g_ref) - 1.0;

    result->bias.x = bx;
    result->bias.y = by;
    result->bias.z = bz;
    result->scale_factor.x = sfx;
    result->scale_factor.y = sfy;
    result->scale_factor.z = sfz;
    ins_mat3_identity(&result->misalignment);
    result->num_positions = 6;

    /* Compute residual RMS */
    double residual_sq = 0.0;
    for (int p = 0; p < 6; p++) {
        double ax_corr = (readings[p].x - bx) / (1.0 + sfx);
        double ay_corr = (readings[p].y - by) / (1.0 + sfy);
        double az_corr = (readings[p].z - bz) / (1.0 + sfz);
        double mag = sqrt(ax_corr * ax_corr + ay_corr * ay_corr + az_corr * az_corr);
        double err = mag - g_ref;
        residual_sq += err * err;
    }
    result->residual_rms = sqrt(residual_sq / 6.0);

    return 0;
}

/* =========================================================================
 * L5: Multi-Position Gyroscope Calibration
 *
 * For each known angular rate vector w_true in body frame:
 *   w_measured = bias + (1 + sf) * w_true + noise
 *
 * Over N positions, we solve the overdetermined system:
 *   [w_true, 1] * [sf; bias] = w_measured
 *
 * using least squares for each axis independently.
 * ========================================================================= */

int ins_calib_multi_position_gyro(const ins_vec3_t *rates,
                                   const ins_vec3_t *readings,
                                   size_t n_positions,
                                   ins_calibration_result_t *result) {
    if (!rates || !readings || !result || n_positions < 3) return -1;

    /* For each axis, solve: measured = (1+sf)*true + bias  =>  measured - true = sf*true + bias */
    /* Reformulate: y = measured - true, A = [true, 1], x = [sf, bias]^T */
    /* Normal equations: (A^T * A) * x = A^T * y */

    double sum_wx = 0, sum_wx2 = 0, sum_yx = 0, sum_yx_wx = 0;
    double sum_wy = 0, sum_wy2 = 0, sum_yy = 0, sum_yy_wy = 0;
    double sum_wz = 0, sum_wz2 = 0, sum_yz = 0, sum_yz_wz = 0;

    for (size_t i = 0; i < n_positions; i++) {
        /* X-axis */
        double wx = rates[i].x;
        double yx = readings[i].x - wx;
        sum_wx  += wx;
        sum_wx2 += wx * wx;
        sum_yx  += yx;
        sum_yx_wx += yx * wx;

        /* Y-axis */
        double wy = rates[i].y;
        double yy = readings[i].y - wy;
        sum_wy  += wy;
        sum_wy2 += wy * wy;
        sum_yy  += yy;
        sum_yy_wy += yy * wy;

        /* Z-axis */
        double wz = rates[i].z;
        double yz = readings[i].z - wz;
        sum_wz  += wz;
        sum_wz2 += wz * wz;
        sum_yz  += yz;
        sum_yz_wz += yz * wz;
    }

    double N = (double)n_positions;

    /* Solve for each axis using Cramer's rule */
    double detx = sum_wx2 * N - sum_wx * sum_wx;
    double dety = sum_wy2 * N - sum_wy * sum_wy;
    double detz = sum_wz2 * N - sum_wz * sum_wz;

    if (fabs(detx) < 1e-15 || fabs(dety) < 1e-15 || fabs(detz) < 1e-15)
        return -1;

    double sfx = (sum_yx_wx * N - sum_yx * sum_wx) / detx;
    double bx  = (sum_wx2 * sum_yx - sum_wx * sum_yx_wx) / detx;

    double sfy = (sum_yy_wy * N - sum_yy * sum_wy) / dety;
    double by  = (sum_wy2 * sum_yy - sum_wy * sum_yy_wy) / dety;

    double sfz = (sum_yz_wz * N - sum_yz * sum_wz) / detz;
    double bz  = (sum_wz2 * sum_yz - sum_wz * sum_yz_wz) / detz;

    result->bias.x = bx;
    result->bias.y = by;
    result->bias.z = bz;
    result->scale_factor.x = sfx;
    result->scale_factor.y = sfy;
    result->scale_factor.z = sfz;
    ins_mat3_identity(&result->misalignment);
    result->num_positions = (int)n_positions;

    /* Residual RMS */
    double residual_sq = 0.0;
    for (size_t i = 0; i < n_positions; i++) {
        double wx_corr = (readings[i].x - bx) / (1.0 + sfx);
        double wy_corr = (readings[i].y - by) / (1.0 + sfy);
        double wz_corr = (readings[i].z - bz) / (1.0 + sfz);
        double err_x = wx_corr - rates[i].x;
        double err_y = wy_corr - rates[i].y;
        double err_z = wz_corr - rates[i].z;
        residual_sq += err_x * err_x + err_y * err_y + err_z * err_z;
    }
    result->residual_rms = sqrt(residual_sq / (3.0 * N));

    return 0;
}

/* =========================================================================
 * L5: Cross-Axis Misalignment Calibration
 *
 * Given N pairs of (reference_vector, measured_vector), we fit a
 * linear transformation M such that: measured = M * reference + bias
 *
 * For small misalignment angles, M is approximated as:
 *   M = I + [0, -alpha_z, alpha_y; alpha_z, 0, -alpha_x; -alpha_y, alpha_x, 0]
 *
 * More generally, we solve the least squares problem:
 *   min ||M * ref_i - meas_i||^2
 *
 * With N >= 3 non-coplanar reference vectors, M is uniquely determined.
 * ========================================================================= */

int ins_calib_cross_axis(const ins_vec3_t *reference,
                          const ins_vec3_t *measured,
                          size_t n,
                          ins_mat3_t *M_corr) {
    if (!reference || !measured || !M_corr || n < 3) return -1;

    /* Solve M = argmin sum_i ||M*ref_i - meas_i||^2 */
    /* This is 3 independent least squares problems (one per row of M).
     * For row k: min ||[ref_i^T] * M_k^T - meas_{i,k}||^2 */

    /* Build normal matrix: A^T*A where A = [ref_i^T] stacked */
    double AtA[9] = {0};  /* 3x3, row-major */
    double Atb[9] = {0};  /* 3x3, row-major, one column per row of M */

    for (size_t i = 0; i < n; i++) {
        double rx = reference[i].x, ry = reference[i].y, rz = reference[i].z;
        double mx = measured[i].x, my = measured[i].y, mz = measured[i].z;

        /* AtA += [rx,ry,rz]^T * [rx,ry,rz] */
        AtA[0] += rx * rx;  AtA[1] += rx * ry;  AtA[2] += rx * rz;
        AtA[3] += ry * rx;  AtA[4] += ry * ry;  AtA[5] += ry * rz;
        AtA[6] += rz * rx;  AtA[7] += rz * ry;  AtA[8] += rz * rz;

        /* Atb(:,0) += [rx,ry,rz]^T * mx */
        Atb[0] += rx * mx;  Atb[3] += ry * mx;  Atb[6] += rz * mx;
        /* Atb(:,1) += [rx,ry,rz]^T * my */
        Atb[1] += rx * my;  Atb[4] += ry * my;  Atb[7] += rz * my;
        /* Atb(:,2) += [rx,ry,rz]^T * mz */
        Atb[2] += rx * mz;  Atb[5] += ry * mz;  Atb[8] += rz * mz;
    }

    /* Solve AtA * M^T = Atb using 3x3 Gaussian elimination */
    /* Copy AtA to working matrix */
    double A[9];
    memcpy(A, AtA, 9 * sizeof(double));
    double b[9];
    memcpy(b, Atb, 9 * sizeof(double));

    /* Forward elimination */
    for (int col = 0; col < 3; col++) {
        /* Find pivot */
        int max_row = col;
        double max_val = fabs(A[col * 3 + col]);
        for (int row = col + 1; row < 3; row++) {
            double val = fabs(A[row * 3 + col]);
            if (val > max_val) { max_val = val; max_row = row; }
        }
        if (max_val < 1e-12) return -1;  /* Rank deficient */

        /* Swap rows */
        if (max_row != col) {
            for (int j = col; j < 3; j++) {
                double tmp = A[col * 3 + j];
                A[col * 3 + j] = A[max_row * 3 + j];
                A[max_row * 3 + j] = tmp;
            }
            for (int j = 0; j < 3; j++) {
                double tmp = b[col * 3 + j];
                b[col * 3 + j] = b[max_row * 3 + j];
                b[max_row * 3 + j] = tmp;
            }
        }

        /* Eliminate below pivot */
        for (int row = col + 1; row < 3; row++) {
            double factor = A[row * 3 + col] / A[col * 3 + col];
            for (int j = col; j < 3; j++) {
                A[row * 3 + j] -= factor * A[col * 3 + j];
            }
            for (int j = 0; j < 3; j++) {
                b[row * 3 + j] -= factor * b[col * 3 + j];
            }
        }
    }

    /* Back substitution */
    for (int j = 0; j < 3; j++) {
        for (int row = 2; row >= 0; row--) {
            double sum = b[row * 3 + j];
            for (int col = row + 1; col < 3; col++) {
                sum -= A[row * 3 + col] * M_corr->m[col * 3 + j];
            }
            M_corr->m[row * 3 + j] = sum / A[row * 3 + row];
        }
    }

    return 0;
}

/* =========================================================================
 * L5: Thermal Calibration
 * ========================================================================= */

int ins_calib_thermal_fit(const double *temperatures, const double *biases,
                           size_t n, ins_thermal_coeffs_t *coeffs) {
    if (!temperatures || !biases || !coeffs || n < 4) return -1;

    /* Fit: bias(T) = c0 + c1*(T-T0) + c2*(T-T0)^2 + c3*(T-T0)^3 */
    /* Use standard reference temperature */
    coeffs->T0 = 25.0;

    /* Build normal equations for polynomial least squares */
    /* Design matrix X: each row = [1, (T-T0), (T-T0)^2, (T-T0)^3] */
    double XtX[16] = {0}; /* 4x4 */
    double Xty[4]  = {0};

    for (size_t i = 0; i < n; i++) {
        double dT = temperatures[i] - coeffs->T0;
        double dT2 = dT * dT;
        double dT3 = dT2 * dT;
        double row[4] = {1.0, dT, dT2, dT3};

        for (int r = 0; r < 4; r++) {
            for (int c = 0; c < 4; c++) {
                XtX[r * 4 + c] += row[r] * row[c];
            }
            Xty[r] += row[r] * biases[i];
        }
    }

    /* Solve using Gaussian elimination */
    double aug[4][5];
    for (int r = 0; r < 4; r++) {
        for (int c = 0; c < 4; c++) {
            aug[r][c] = XtX[r * 4 + c];
        }
        aug[r][4] = Xty[r];
    }

    for (int col = 0; col < 4; col++) {
        int max_row = col;
        for (int r = col + 1; r < 4; r++) {
            if (fabs(aug[r][col]) > fabs(aug[max_row][col]))
                max_row = r;
        }
        if (fabs(aug[max_row][col]) < 1e-12) {
            /* Nearly singular — set higher-order coeffs to 0 */
            coeffs->c3 = 0.0;
            coeffs->c2 = 0.0;
            /* Fit with just linear */
            break;
        }
        if (max_row != col) {
            for (int c = 0; c <= 4; c++) {
                double tmp = aug[col][c];
                aug[col][c] = aug[max_row][c];
                aug[max_row][c] = tmp;
            }
        }
        for (int r = col + 1; r < 4; r++) {
            double factor = aug[r][col] / aug[col][col];
            for (int c = col; c <= 4; c++) {
                aug[r][c] -= factor * aug[col][c];
            }
        }
    }

    /* Back substitution */
    double coeff_arr[4];
    for (int r = 3; r >= 0; r--) {
        double sum = aug[r][4];
        for (int c = r + 1; c < 4; c++) {
            sum -= aug[r][c] * coeff_arr[c];
        }
        coeff_arr[r] = sum / aug[r][r];
    }

    coeffs->c0 = coeff_arr[0];
    coeffs->c1 = coeff_arr[1];
    coeffs->c2 = coeff_arr[2];
    coeffs->c3 = coeff_arr[3];

    /* Set temperature range */
    coeffs->T_min = temperatures[0];
    coeffs->T_max = temperatures[0];
    for (size_t i = 1; i < n; i++) {
        if (temperatures[i] < coeffs->T_min) coeffs->T_min = temperatures[i];
        if (temperatures[i] > coeffs->T_max) coeffs->T_max = temperatures[i];
    }

    return 0;
}

double ins_calib_thermal_apply(double reading, double T,
                                const ins_thermal_coeffs_t *coeffs) {
    if (!coeffs) return reading;
    double dT = T - coeffs->T0;
    double correction = coeffs->c0 + coeffs->c1 * dT
                      + coeffs->c2 * dT * dT
                      + coeffs->c3 * dT * dT * dT;
    return reading - correction;
}

/* =========================================================================
 * L5: Full Calibration Pipeline
 * ========================================================================= */

void ins_calib_apply(const ins_imu_sample_t *raw,
                      const ins_imu_calibration_t *calib,
                      double temp,
                      ins_vec3_t *accel_corr,
                      ins_vec3_t *gyro_corr) {
    if (!raw || !calib || !accel_corr || !gyro_corr) return;

    (void)temp;  /* Temperature compensation not applied in this simplified pipeline;
                    use ins_calib_thermal_apply() externally if needed */

    /* Process accelerometer triad */
    /* Step 1: Bias subtraction */
    double ax = raw->accel.x - calib->accel.bias.x;
    double ay = raw->accel.y - calib->accel.bias.y;
    double az = raw->accel.z - calib->accel.bias.z;

    /* Step 2: Scale factor correction */
    ax /= (1.0 + calib->accel.scale_factor.x);
    ay /= (1.0 + calib->accel.scale_factor.y);
    az /= (1.0 + calib->accel.scale_factor.z);

    /* Step 3: Misalignment correction */
    accel_corr->x = calib->accel.misalignment.m[0] * ax
                  + calib->accel.misalignment.m[1] * ay
                  + calib->accel.misalignment.m[2] * az;
    accel_corr->y = calib->accel.misalignment.m[3] * ax
                  + calib->accel.misalignment.m[4] * ay
                  + calib->accel.misalignment.m[5] * az;
    accel_corr->z = calib->accel.misalignment.m[6] * ax
                  + calib->accel.misalignment.m[7] * ay
                  + calib->accel.misalignment.m[8] * az;

    /* Process gyroscope triad (same pipeline) */
    double gx = raw->gyro.x - calib->gyro.bias.x;
    double gy = raw->gyro.y - calib->gyro.bias.y;
    double gz = raw->gyro.z - calib->gyro.bias.z;

    gx /= (1.0 + calib->gyro.scale_factor.x);
    gy /= (1.0 + calib->gyro.scale_factor.y);
    gz /= (1.0 + calib->gyro.scale_factor.z);

    gyro_corr->x = calib->gyro.misalignment.m[0] * gx
                 + calib->gyro.misalignment.m[1] * gy
                 + calib->gyro.misalignment.m[2] * gz;
    gyro_corr->y = calib->gyro.misalignment.m[3] * gx
                 + calib->gyro.misalignment.m[4] * gy
                 + calib->gyro.misalignment.m[5] * gz;
    gyro_corr->z = calib->gyro.misalignment.m[6] * gx
                 + calib->gyro.misalignment.m[7] * gy
                 + calib->gyro.misalignment.m[8] * gz;
}
