/**
 * @file    ins_integration.c
 * @brief   GPS/INS loosely-coupled integration with error-state Kalman filter
 *
 * Knowledge Coverage:
 *   L5 (Algorithms): Error-state Kalman filter (15-state), GPS update
 *   L6 (Canonical Problems): Loosely coupled INS/GPS integration
 *   L7 (Applications): ZUPT for pedestrian navigation, GNSS outage prediction
 *
 * Reference:
 *   Groves (2013), Chapters 14-16, "INS/GNSS Integration".
 *   Farrell (2008), "Aided Navigation: GPS with High Rate Sensors",
 *     McGraw-Hill, Chapters 8-11.
 *   Skog et al. (2010), IEEE Trans. Biomed. Eng., 57(11): 2657-2666.
 *
 * Course Mapping:
 *   Stanford AA272 - GPS (Kalman filter for navigation)
 *   MIT 2.171 - Precision Machine Design (sensor fusion)
 *   Michigan EECS 455 - Communications (estimation)
 */

#include "ins_integration.h"
#include "ins_attitude.h"
#include <string.h>
#include <math.h>
#include <stdio.h>

/* =========================================================================
 * L5: Kalman Filter Initialization
 *
 * The 15-state error-state vector:
 *   x[0..2]   = position error (N, E, D) [m]
 *   x[3..5]   = velocity error (N, E, D) [m/s]
 *   x[6..8]   = attitude error (psi_x, psi_y, psi_z) [rad]
 *   x[9..11]  = accelerometer bias error [m/s^2]
 *   x[12..14] = gyroscope bias error [rad/s]
 *
 * Initial covariance P is diagonal with variances corresponding to
 * expected initial navigation uncertainties.
 * ========================================================================= */

void ins_kf_init(ins_kf_state_t *kf,
                 double sigma_pos, double sigma_vel, double sigma_att) {
    if (!kf) return;

    memset(kf, 0, sizeof(*kf));

    /* Initialize P as diagonal */
    double sigma2_pos = sigma_pos * sigma_pos;
    double sigma2_vel = sigma_vel * sigma_vel;
    double sigma2_att = sigma_att * sigma_att;

    /* Typical values for accel/gyro bias uncertainties */
    double sigma2_acc_bias = 0.01 * 0.01;  /* 1 cm/s^2 */
    double sigma2_gyro_bias = 1e-5;        /* ~5.7 deg/hr */

    for (int i = 0; i < 3; i++) {
        kf->P[i * 15 + i] = sigma2_pos;
        kf->P[(3 + i) * 15 + (3 + i)] = sigma2_vel;
        kf->P[(6 + i) * 15 + (6 + i)] = sigma2_att;
        kf->P[(9 + i) * 15 + (9 + i)] = sigma2_acc_bias;
        kf->P[(12 + i) * 15 + (12 + i)] = sigma2_gyro_bias;
    }

    kf->time = 0.0;
    kf->initialized = 1;
}

/* =========================================================================
 * L5: Kalman Filter Predict Step
 *
 * Propagates the error state covariance P using the state transition
 * matrix Phi from the Psi-angle error model.
 *
 * P_pred = Phi * P * Phi^T + Q_d
 *
 * The process noise Q_d accounts for sensor random walk:
 *   Q_d(accel bias) = Q_accel_psd * dt
 *   Q_d(gyro bias)  = Q_gyro_psd * dt
 *
 * The error state itself is NOT propagated (x stays at zero after
 * closed-loop correction, or accumulates small residual errors).
 * ========================================================================= */

void ins_kf_predict(ins_kf_state_t *kf,
                    const ins_mat3_t *C_body_ned,
                    const ins_vec3_t *f_body,
                    double lat, double alt, double dt,
                    double Q_accel_psd, double Q_gyro_psd) {
    if (!kf || !kf->initialized || !C_body_ned || !f_body || dt <= 0) return;

    /* Dummy state vector for computing Phi (we only need Phi, not state propagation) */
    ins_vec3_t dummy_pos = {0}, dummy_vel = {0}, dummy_psi = {0};

    /* Get Phi matrix */
    double Phi[81];
    ins_psi_error_propagate(&dummy_pos, &dummy_vel, &dummy_psi,
                             C_body_ned, f_body, lat, alt, dt, Phi);

    /* P_pred = Phi * P * Phi^T */
    double P_temp[225];

    /* P_temp = Phi * P */
    for (int i = 0; i < 15; i++) {
        for (int j = 0; j < 15; j++) {
            double sum = 0.0;
            for (int k = 0; k < 15; k++) {
                sum += Phi[i * 15 + k] * kf->P[k * 15 + j];
            }
            P_temp[i * 15 + j] = sum;
        }
    }

    /* P = P_temp * Phi^T */
    for (int i = 0; i < 15; i++) {
        for (int j = 0; j < 15; j++) {
            double sum = 0.0;
            for (int k = 0; k < 15; k++) {
                sum += P_temp[i * 15 + k] * Phi[j * 15 + k];
            }
            kf->P[i * 15 + j] = sum;
        }
    }

    /* Add process noise: Q_d diagonal */
    double q_acc = Q_accel_psd * dt;
    double q_gyr = Q_gyro_psd * dt;
    for (int i = 9; i < 12; i++) {
        kf->P[i * 15 + i] += q_acc;
    }
    for (int i = 12; i < 15; i++) {
        kf->P[i * 15 + i] += q_gyr;
    }

    kf->time += dt;
}

/* =========================================================================
 * L5: Kalman Filter GPS Update (Loosely Coupled)
 *
 * Measurement vector: z = [pos_gps - pos_ins, vel_gps - vel_ins]^T (6x1)
 *
 * Measurement model: z = H * x + v
 *   H = [I_3, 0_3, 0_3, 0_3, 0_3;
 *        0_3, I_3, 0_3, 0_3, 0_3]
 *
 * Measurement noise R (6x6 diagonal from GPS std):
 *   R = diag(sigma_pos^2, sigma_vel^2)
 * ========================================================================= */

int ins_kf_update_gps_pos_vel(ins_kf_state_t *kf,
                               const ins_gps_measurement_t *gps,
                               const ins_nav_solution_t *ins) {
    if (!kf || !kf->initialized || !gps || !ins) return -1;

    /* Innovation: z = measurement - prediction */
    double z[6];

    /* Position innovation (in meters, approximately) */
    double M = ins_meridian_radius(ins->pos.lat);
    double N = ins_prime_vertical_radius(ins->pos.lat);
    double cos_lat = cos(ins->pos.lat);

    z[0] = (gps->pos.lat - ins->pos.lat) * (M + ins->pos.alt);
    z[1] = (gps->pos.lon - ins->pos.lon) * (N + ins->pos.alt) * cos_lat;
    z[2] = gps->pos.alt - ins->pos.alt;

    /* Velocity innovation */
    z[3] = gps->vel_ned.x - ins->vel_ned.x;
    z[4] = gps->vel_ned.y - ins->vel_ned.y;
    z[5] = gps->vel_ned.z - ins->vel_ned.z;

    /* H matrix: H = [I_3, 0, 0; 0, I_3, 0] for 6x15 */
    /* S = H * P * H^T */
    double S[36] = {0};

    /* Build S from relevant P sub-blocks */
    for (int i = 0; i < 6; i++) {
        int row_p = (i < 3) ? i : (i - 3 + 3);
        for (int j = 0; j < 6; j++) {
            int col_p = (j < 3) ? j : (j - 3 + 3);
            S[i * 6 + j] = kf->P[row_p * 15 + col_p];
        }
    }

    /* Add measurement noise R */
    for (int i = 0; i < 3; i++) {
        S[i * 6 + i] += gps->pos_std[i] * gps->pos_std[i];
        S[(3 + i) * 6 + (3 + i)] += gps->vel_std[i] * gps->vel_std[i];
    }

    /* Compute Kalman gain K = P * H^T * S^(-1) */
    /* First compute P * H^T (15x6) = columns 0..2 and 3..5 of P */
    double P_HT[90] = {0};
    for (int i = 0; i < 15; i++) {
        for (int j = 0; j < 3; j++) {
            P_HT[i * 6 + j] = kf->P[i * 15 + j];
            P_HT[i * 6 + (3 + j)] = kf->P[i * 15 + (3 + j)];
        }
    }

    /* Solve S * K^T = (P*H^T)^T for each row of K^T */
    /* We solve S * x = b for 15 right-hand sides */
    double K[90] = {0};  /* 15x6, K(transpose of what we need) */
    for (int col = 0; col < 15; col++) {
        /* For each column of K^T (row of K), solve S * k_col = P_HT_col */
        /* Use Gaussian elimination on S (6x6) */
        double S_work[36];
        memcpy(S_work, S, 36 * sizeof(double));
        double b_work[6];
        for (int i = 0; i < 6; i++) {
            b_work[i] = P_HT[col * 6 + i];
        }

        /* Forward elimination */
        for (int r = 0; r < 6; r++) {
            int max_r = r;
            for (int r2 = r + 1; r2 < 6; r2++) {
                if (fabs(S_work[r2 * 6 + r]) > fabs(S_work[max_r * 6 + r]))
                    max_r = r2;
            }
            if (fabs(S_work[max_r * 6 + r]) < 1e-14) continue;

            if (max_r != r) {
                for (int c = r; c < 6; c++) {
                    double tmp = S_work[r * 6 + c];
                    S_work[r * 6 + c] = S_work[max_r * 6 + c];
                    S_work[max_r * 6 + c] = tmp;
                }
                double tmp = b_work[r];
                b_work[r] = b_work[max_r];
                b_work[max_r] = tmp;
            }

            double pivot = S_work[r * 6 + r];
            for (int r2 = r + 1; r2 < 6; r2++) {
                double factor = S_work[r2 * 6 + r] / pivot;
                for (int c = r; c < 6; c++) {
                    S_work[r2 * 6 + c] -= factor * S_work[r * 6 + c];
                }
                b_work[r2] -= factor * b_work[r];
            }
        }

        /* Back substitution */
        double k_row[6];
        for (int r = 5; r >= 0; r--) {
            double sum = b_work[r];
            for (int c = r + 1; c < 6; c++) {
                sum -= S_work[r * 6 + c] * k_row[c];
            }
            if (fabs(S_work[r * 6 + r]) > 1e-14) {
                k_row[r] = sum / S_work[r * 6 + r];
            } else {
                k_row[r] = 0.0;
            }
        }

        /* Store: K[col, :] = k_row */
        for (int i = 0; i < 6; i++) {
            K[col * 6 + i] = k_row[i];
        }
    }

    /* State update: x = x + K * z */
    for (int i = 0; i < 15; i++) {
        double Kz = 0.0;
        for (int j = 0; j < 6; j++) {
            Kz += K[i * 6 + j] * z[j];
        }
        kf->x[i] += Kz;
    }

    /* Covariance update: P = (I - K*H) * P */
    /* First compute I - K*H */
    double I_KH[225];
    memcpy(I_KH, kf->P, 225 * sizeof(double));

    /* P = P - K * H * P = P - K * [P(0..2,:); P(3..5,:)] */
    double P_new[225];
    memcpy(P_new, kf->P, 225 * sizeof(double));
    for (int i = 0; i < 15; i++) {
        for (int j = 0; j < 15; j++) {
            double sum = 0.0;
            for (int k = 0; k < 3; k++) {
                sum += K[i * 6 + k] * kf->P[k * 15 + j];
                sum += K[i * 6 + (3 + k)] * kf->P[(3 + k) * 15 + j];
            }
            P_new[i * 15 + j] -= sum;
        }
    }
    memcpy(kf->P, P_new, 225 * sizeof(double));

    return 0;
}

/* =========================================================================
 * L5: Apply Kalman Correction to Navigation Solution
 * ========================================================================= */

void ins_kf_apply_correction(ins_kf_state_t *kf, ins_nav_solution_t *ins) {
    if (!kf || !ins) return;

    /* Position correction */
    double M = ins_meridian_radius(ins->pos.lat);
    double N = ins_prime_vertical_radius(ins->pos.lat);
    double cos_lat = cos(ins->pos.lat);

    ins->pos.lat -= kf->x[0] / (M + ins->pos.alt);
    ins->pos.lon -= kf->x[1] / ((N + ins->pos.alt) * cos_lat);
    ins->pos.alt -= kf->x[2];

    /* Velocity correction */
    ins->vel_ned.x -= kf->x[3];
    ins->vel_ned.y -= kf->x[4];
    ins->vel_ned.z -= kf->x[5];

    /* Attitude correction: rotate by -psi */
    /* q_corrected = q_ins * q_delta where q_delta ~= [1, -psi_x/2, -psi_y/2, -psi_z/2] */
    ins_quat_t q_ins, q_corr, q_result;
    q_ins.w = ins->q[0]; q_ins.x = ins->q[1]; q_ins.y = ins->q[2]; q_ins.z = ins->q[3];

    double hx = -kf->x[6] * 0.5;
    double hy = -kf->x[7] * 0.5;
    double hz = -kf->x[8] * 0.5;
    double hn = sqrt(hx * hx + hy * hy + hz * hz);
    if (hn < 1e-15) {
        ins_quat_identity(&q_corr);
    } else {
        double sin_hn = sin(hn);
        q_corr.w = cos(hn);
        q_corr.x = (hx / hn) * sin_hn;
        q_corr.y = (hy / hn) * sin_hn;
        q_corr.z = (hz / hn) * sin_hn;
    }
    ins_quat_mul(&q_ins, &q_corr, &q_result);
    ins_quat_normalize(&q_result);
    ins->q[0] = q_result.w; ins->q[1] = q_result.x;
    ins->q[2] = q_result.y; ins->q[3] = q_result.z;

    /* Update Euler angles */
    ins_euler_t att;
    ins_quat_to_euler(&q_result, &att);
    ins->roll = att.roll;
    ins->pitch = att.pitch;
    ins->yaw = att.yaw;

    /* Reset position, velocity, attitude errors to zero (closed-loop) */
    for (int i = 0; i < 9; i++) {
        kf->x[i] = 0.0;
    }
}

/* =========================================================================
 * L6: Loosely Coupled Integration Loop
 * ========================================================================= */

int ins_integrate_loose(ins_mech_state_t *mech_state,
                         ins_kf_state_t *kf,
                         const ins_imu_sample_t *imu_data,
                         size_t num_imu,
                         const ins_gps_measurement_t *gps_data,
                         size_t num_gps,
                         size_t imu_per_gps,
                         ins_nav_solution_t *output) {
    if (!mech_state || !kf || !imu_data || !output || num_imu == 0)
        return -1;

    /* Make mutable copies */
    ins_mech_state_t ms;
    memcpy(&ms, mech_state, sizeof(ms));
    ins_kf_state_t kf_local;
    memcpy(&kf_local, kf, sizeof(kf_local));

    /* Process noise PSDs (typical tactical-grade IMU) */
    double Qa = 1e-6;   /* Accel VRW: 0.001 m/s/sqrt(hr) ~ 1.7e-5 m/s^2/sqrt(Hz) */
    double Qg = 1e-12;  /* Gyro ARW: 0.001 deg/sqrt(hr) ~ 1.7e-7 rad/sqrt(Hz) */

    size_t gps_idx = 0;

    for (size_t i = 0; i < num_imu; i++) {
        /* INS mechanization propagation */
        int ret = ins_mech_step(&ms, &imu_data[i]);
        if (ret != 0) return -1;

        /* Get current INS solution */
        ins_nav_solution_t nav;
        ins_mech_get_solution(&ms, &nav);

        /* Kalman filter predict */
        ins_mat3_t C_bn;
        ins_quat_to_dcm(&ms.quat, &C_bn);
        ins_kf_predict(&kf_local, &C_bn, &imu_data[i].accel,
                       ms.pos.lat, ms.pos.alt, imu_data[i].dt, Qa, Qg);

        /* GPS update if it's time */
        if (gps_data && gps_idx < num_gps &&
            (i + 1) % imu_per_gps == 0) {
            ins_kf_update_gps_pos_vel(&kf_local, &gps_data[gps_idx], &nav);
            ins_kf_apply_correction(&kf_local, &nav);

            /* Apply corrections back to mechanization state */
            ms.pos = nav.pos;
            ins_vec3_copy(&nav.vel_ned, &ms.vel_ned);
            ms.quat.w = nav.q[0]; ms.quat.x = nav.q[1];
            ms.quat.y = nav.q[2]; ms.quat.z = nav.q[3];

            gps_idx++;
        }

        output[i] = nav;
    }

    return 0;
}

/* =========================================================================
 * L6: Zero-Velocity Detection (SHOE)
 * ========================================================================= */

int ins_zupt_detect(const ins_vec3_t *accel, const ins_vec3_t *gyro,
                     size_t window_len,
                     double sigma_a, double sigma_w,
                     double threshold) {
    if (!accel || !gyro || window_len < 1) return 0;

    /* Compute average gravity direction */
    double g_avg_x = 0, g_avg_y = 0, g_avg_z = 0;
    for (size_t i = 0; i < window_len; i++) {
        g_avg_x += accel[i].x;
        g_avg_y += accel[i].y;
        g_avg_z += accel[i].z;
    }
    double g_norm = sqrt(g_avg_x * g_avg_x + g_avg_y * g_avg_y + g_avg_z * g_avg_z);
    if (g_norm < 1e-12) return 0;
    g_avg_x /= g_norm;
    g_avg_y /= g_norm;
    g_avg_z /= g_norm;

    /* Compute test statistic */
    double T = 0.0;
    for (size_t i = 0; i < window_len; i++) {
        /* Accelerometer deviation from gravity */
        double ax = accel[i].x, ay = accel[i].y, az = accel[i].z;
        double a_mag = sqrt(ax * ax + ay * ay + az * az);
        double a_dev = a_mag - INS_GRAVITY_EQUATOR;
        double a_term = (a_dev * a_dev) / (sigma_a * sigma_a);

        /* Gyroscope energy */
        double wx = gyro[i].x, wy = gyro[i].y, wz = gyro[i].z;
        double w_term = (wx * wx + wy * wy + wz * wz) / (sigma_w * sigma_w);

        T += a_term + w_term;
    }
    T /= (double)window_len;

    return (T < threshold) ? 1 : 0;
}

/* =========================================================================
 * L6: ZUPT Kalman Update
 * ========================================================================= */

int ins_kf_update_zupt(ins_kf_state_t *kf, const ins_nav_solution_t *ins) {
    if (!kf || !ins) return -1;

    /* Innovation: measured velocity should be zero */
    double z[3];
    z[0] = -ins->vel_ned.x;
    z[1] = -ins->vel_ned.y;
    z[2] = -ins->vel_ned.z;

    /* Measurement noise for ZUPT: ~0.01 m/s std */
    double zupt_var = 0.0001;  /* (0.01 m/s)^2 */

    /* H = [0_3, I_3, 0_3, 0_3, 0_3] */
    /* S = HPH^T + R = P(3..5, 3..5) + R */
    double S[9];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            S[i * 3 + j] = kf->P[(3 + i) * 15 + (3 + j)];
            if (i == j) S[i * 3 + j] += zupt_var;
        }
    }

    /* K = P*H^T * S^(-1) = P(:, 3..5) * S^(-1) */
    /* Compute S inverse (3x3) */
    double det = S[0] * (S[4] * S[8] - S[5] * S[7])
               - S[1] * (S[3] * S[8] - S[5] * S[6])
               + S[2] * (S[3] * S[7] - S[4] * S[6]);

    if (fabs(det) < 1e-20) return -1;

    double Sinv[9];
    double inv_det = 1.0 / det;
    Sinv[0] = (S[4] * S[8] - S[5] * S[7]) * inv_det;
    Sinv[1] = (S[2] * S[7] - S[1] * S[8]) * inv_det;
    Sinv[2] = (S[1] * S[5] - S[2] * S[4]) * inv_det;
    Sinv[3] = (S[5] * S[6] - S[3] * S[8]) * inv_det;
    Sinv[4] = (S[0] * S[8] - S[2] * S[6]) * inv_det;
    Sinv[5] = (S[2] * S[3] - S[0] * S[5]) * inv_det;
    Sinv[6] = (S[3] * S[7] - S[4] * S[6]) * inv_det;
    Sinv[7] = (S[1] * S[6] - S[0] * S[7]) * inv_det;
    Sinv[8] = (S[0] * S[4] - S[1] * S[3]) * inv_det;

    /* K = P(:,3..5) * Sinv */
    double K_zupt[45] = {0};  /* 15x3 */
    for (int i = 0; i < 15; i++) {
        for (int j = 0; j < 3; j++) {
            double sum = 0.0;
            for (int k = 0; k < 3; k++) {
                sum += kf->P[i * 15 + (3 + k)] * Sinv[k * 3 + j];
            }
            K_zupt[i * 3 + j] = sum;
        }
    }

    /* State update */
    for (int i = 0; i < 15; i++) {
        double Kz = 0.0;
        for (int j = 0; j < 3; j++) {
            Kz += K_zupt[i * 3 + j] * z[j];
        }
        kf->x[i] += Kz;
    }

    /* Covariance update: P = (I - K*H) * P */
    /* K*H only affects columns 3..5 of P */
    double P_new[225];
    memcpy(P_new, kf->P, 225 * sizeof(double));
    for (int i = 0; i < 15; i++) {
        for (int j = 0; j < 15; j++) {
            double sum = 0.0;
            for (int k = 0; k < 3; k++) {
                sum += K_zupt[i * 3 + k] * kf->P[(3 + k) * 15 + j];
            }
            P_new[i * 15 + j] -= sum;
        }
    }
    memcpy(kf->P, P_new, 225 * sizeof(double));

    return 0;
}

/* =========================================================================
 * L7: GNSS Outage Budget Estimation
 *
 * Computes maximum time that INS can operate without GNSS while
 * maintaining specified position accuracy.
 *
 * Uses the cubic-growth gyro dominance model:
 *   pos_error(t) = (1/6) * gyro_bias * g * t^3
 *
 * Solving for t: t = (6 * max_error / (gyro_bias * g))^(1/3)
 *
 * Realistic examples:
 *   Navigation-grade (0.005 deg/hr gyro): 1 nm drift ~ 6 hours
 *   Tactical-grade (1 deg/hr gyro): 1 nm drift ~ 6-8 minutes
 *   Consumer-grade (3600 deg/hr gyro): 1 nm drift ~ seconds (effectively 0)
 * ========================================================================= */

double ins_gnss_outage_budget(const ins_imu_error_model_t *model,
                               double max_pos_error) {
    if (!model || max_pos_error <= 0) return 0.0;

    /* Use worst-axis gyro bias */
    double gyro_bias_max = 0.0;
    const ins_axis_error_t *gyro[3] = {&model->gyro_x, &model->gyro_y, &model->gyro_z};
    for (int i = 0; i < 3; i++) {
        if (fabs(gyro[i]->bias_offset) > gyro_bias_max) {
            gyro_bias_max = fabs(gyro[i]->bias_offset);
        }
    }
    if (gyro_bias_max < 1e-15) return 1e10;  /* Perfect gyros — infinite time */

    double g = INS_GRAVITY_EQUATOR;
    return pow(6.0 * max_pos_error / (gyro_bias_max * g), 1.0 / 3.0);
}
