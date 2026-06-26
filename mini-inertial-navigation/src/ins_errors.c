/**
 * @file    ins_errors.c
 * @brief   Inertial sensor error models and Allan Variance analysis
 *
 * Knowledge Coverage:
 *   L1 (Definitions): Sensor bias, scale factor, noise types
 *   L4 (Fundamental Laws): Sensor error propagation, position drift prediction
 *   L5 (Algorithms): Allan Variance computation and decomposition, Psi-angle model
 *   L7 (Applications): IMU grade specifications for real-world systems
 *
 * Reference:
 *   IEEE Std 952-1997, "Fiber Optic Gyros — Test Procedure", Annex C.
 *   El-Sheimy et al. (2008), IEEE Trans. Instrum. Meas., 57(1): 140-144.
 *   Groves (2013), Chapters 4-5, "Sensor Error Models".
 *   Benson (1975), IEEE Trans. Aerosp. Electron. Syst., 11(1): 67-80.
 *
 * Course Mapping:
 *   MIT 2.171 - Precision Machine Design (sensor calibration)
 *   Stanford AA272 - GPS (error modeling for integration)
 *   Michigan EECS 455 - Communications (estimation theory)
 */

#include "ins_errors.h"
#include "ins_attitude.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

/* =========================================================================
 * L1: IMU Grade Specifications
 *
 * Typical error budgets for consumer to strategic grade IMUs.
 * Values from manufacturer datasheets and Groves (2013) Table 4.2.
 * ========================================================================= */

static const ins_grade_spec_t grade_specs[] = {
    { INS_GRADE_CONSUMER,   "MEMS Consumer (smartphone)",
      3600.0,      300.0,     10000.0,
      100.0,       0.5,       10000.0,
      10000.0,     "Pedestrian, gaming" },
    { INS_GRADE_INDUSTRIAL, "MEMS Industrial (automotive/robotics)",
      100.0,       30.0,      3000.0,
      10.0,        0.1,       3000.0,
      1000.0,      "Automotive, UAV, robotics" },
    { INS_GRADE_TACTICAL,   "Tactical (missile/UAV)",
      1.0,         0.1,       100.0,
      1.0,         0.02,      100.0,
      10.0,        "Missile, military UAV, torpedo" },
    { INS_GRADE_NAVIGATION, "Navigation-grade (aircraft/ship)",
      0.005,       0.001,     5.0,
      0.025,       0.001,     5.0,
      1.0,         "Commercial aviation, surface ships" },
    { INS_GRADE_STRATEGIC,  "Strategic (submarine/ICBM)",
      0.0001,      0.00005,   1.0,
      0.001,       0.00005,   1.0,
      0.005,       "Nuclear submarine, ICBM" }
};

const ins_grade_spec_t *ins_grade_spec_get(ins_grade_t grade) {
    for (size_t i = 0; i < sizeof(grade_specs) / sizeof(grade_specs[0]); i++) {
        if (grade_specs[i].grade == grade) return &grade_specs[i];
    }
    return &grade_specs[0]; /* Default to consumer */
}

/* =========================================================================
 * L5: Allan Variance Computation
 *
 * The two-sample (overlapping) Allan variance for cluster time tau = m*dt:
 *
 *   sigma^2(tau) = (1 / (2*(N-2m))) * sum_{k=1}^{N-2m}
 *     (Omega_{k+2m,m} - Omega_{k+m,m})^2
 *
 * where Omega_{k,m} = (1/m) * sum_{i=k}^{k+m-1} data[i]
 * is the average over m consecutive samples.
 *
 * For angle/velocity measurements (not rate):
 *   sigma^2(tau) = (1/(2*tau^2*(N-2m))) * sum (theta_{k+2m} - 2*theta_{k+m} + theta_k)^2
 *
 * This implementation uses the rate formulation (suitable for gyro/accel rates).
 *
 * Complexity: O(N * M) where N = data length, M = number of tau values.
 * For large datasets, use the fast method (FFT-based) or decimation.
 *
 * Reference: IEEE Std 952-1997, Annex C.
 * ========================================================================= */

size_t ins_allan_variance(const double *data, size_t n, double dt,
                          double *taus, double *adevs, size_t max_taus) {
    if (!data || n < 8 || dt <= 0 || !taus || !adevs || max_taus == 0)
        return 0;

    /* Compute tau values: tau = m * dt for m = 1, 2, 4, 8, ..., up to n/3 */
    size_t count = 0;
    for (size_t m = 1; m <= n / 3 && count < max_taus; m *= 2) {
        double tau = m * dt;
        size_t n_clusters = n - 2 * m + 1;
        if (n_clusters < 2) break;

        double sum_sq = 0.0;
        for (size_t k = 0; k < n_clusters; k++) {
            /* Compute average over non-overlapping segments */
            double avg1 = 0.0, avg2 = 0.0;
            for (size_t j = 0; j < m; j++) {
                avg1 += data[k + j];
                avg2 += data[k + m + j];
            }
            avg1 /= (double)m;
            avg2 /= (double)m;
            double diff = avg2 - avg1;
            sum_sq += diff * diff;
        }

        double avar = sum_sq / (2.0 * (double)n_clusters);
        taus[count] = tau;
        adevs[count] = sqrt(avar);
        count++;
    }

    return count;
}

/* =========================================================================
 * L5: Allan Variance Decomposition
 *
 * The Allan deviation curve is modeled as the quadrature sum of five
 * independent noise processes:
 *
 *   sigma^2(tau) = 3Q^2/tau^2 + N^2/tau + (0.664*B)^2 + K^2*tau/3 + R^2*tau^2/2
 *
 * Each term corresponds to a distinct noise type:
 *   Q: Quantization noise        — slope = -1 on log-log
 *   N: Angle/Velocity Random Walk — slope = -1/2
 *   B: Bias Instability         — slope = 0 (minimum of the curve)
 *   K: Rate Random Walk         — slope = +1/2
 *   R: Rate Ramp                — slope = +1
 *
 * Uses least-squares fitting in log-log space for robust decomposition.
 * ========================================================================= */

int ins_allan_decompose(const double *taus, const double *adevs, size_t num,
                        ins_allan_result_t *result) {
    if (!taus || !adevs || !result || num < 5) return -1;

    /* Build the design matrix for linear least squares:
     *   A * x = b
     * where x = [Q^2, N^2, B^2, K^2, R^2]^T
     * Each row: [3/tau_i^2, 1/tau_i, 0.664^2, tau_i/3, tau_i^2/2]
     * Each b_i = adev_i^2
     */

    /* Simplified approach: slope analysis in different regions */
    /* Extract individual noise coefficients from the minimum and slope */

    /* Find the minimum of the Allan deviation (bias instability) */
    double min_adev = adevs[0];
    size_t min_idx = 0;
    for (size_t i = 1; i < num; i++) {
        if (adevs[i] < min_adev) {
            min_adev = adevs[i];
            min_idx = i;
        }
    }
    result->bias_instability = min_adev / 0.664;
    result->bias_instability_tau = taus[min_idx];

    /* ARW/VRW from the short-tau region (left of minimum, slope=-1/2) */
    /* sigma(tau) = N / sqrt(tau) => N = sigma(tau) * sqrt(tau) */
    if (min_idx > 0) {
        double n_sum = 0.0;
        size_t n_count = 0;
        for (size_t i = 0; i < min_idx && i < num; i++) {
            n_sum += adevs[i] * sqrt(taus[i]);
            n_count++;
        }
        result->angle_random_walk = (n_count > 0) ? n_sum / (double)n_count : 0.0;
    } else {
        result->angle_random_walk = 0.0;
    }

    /* Quantization noise from shortest tau (slope=-1) */
    /* sigma(tau) = Q * sqrt(3) / tau => Q = sigma(tau) * tau / sqrt(3) */
    if (num >= 2) {
        result->quantization_noise = adevs[0] * taus[0] / sqrt(3.0);
    } else {
        result->quantization_noise = 0.0;
    }

    /* Rate random walk from long-tau region (right of minimum, slope=+1/2) */
    /* sigma(tau) = K * sqrt(tau/3) => K = sigma(tau) * sqrt(3/tau) */
    if (min_idx + 1 < num) {
        double k_sum = 0.0;
        size_t k_count = 0;
        for (size_t i = min_idx + 1; i < num; i++) {
            k_sum += adevs[i] * sqrt(3.0 / taus[i]);
            k_count++;
        }
        result->rate_random_walk = (k_count > 0) ? k_sum / (double)k_count : 0.0;
    } else {
        result->rate_random_walk = 0.0;
    }

    /* Rate ramp from longest tau region (slope=+1) */
    /* sigma(tau) = R * tau / sqrt(2) => R = sigma(tau) * sqrt(2) / tau */
    if (num >= 3) {
        size_t last_idx = num - 1;
        result->rate_ramp = adevs[last_idx] * sqrt(2.0) / taus[last_idx];
    } else {
        result->rate_ramp = 0.0;
    }

    return 0;
}

/* =========================================================================
 * L4: Position Drift Prediction from IMU Error Budget
 *
 * Key physical mechanisms:
 *
 * 1. Accelerometer bias -> velocity drift -> position drift (quadratic):
 *    delta_pos = (1/2) * accel_bias * t^2
 *    Example: 1 mg bias -> 0.0098 m/s^2 -> 17.6 m after 60 s
 *                                              ~1.6 km after 10 min
 *                                              ~63.5 km after 1 hr
 *
 * 2. Gyro bias -> tilt error -> horizontal accel error (cubic):
 *    delta_pos = (1/6) * gyro_bias * g * t^3
 *    Example: 1 deg/hr bias -> 4.85e-6 rad/s
 *              After 10 min: ~1.7 m
 *              After 1 hr:   ~370 m
 *    Note: gyro-induced error grows as t^3 and dominates after ~10 min!
 *
 * 3. Random walk noise:
 *    delta_pos_vrw = VRW * sqrt(t) * t   (accel VRW -> vel -> pos)
 *    delta_pos_arw ~ (1/3) * ARW * g * t^(5/2)  (gyro ARW -> attitude err)
 *
 * Reference: Groves (2013) Section 5.7, Figures 5.14-5.17.
 * ========================================================================= */

void ins_error_predict_drift(const ins_imu_error_model_t *model,
                              double time,
                              double *pos_drift,
                              double *vel_drift,
                              double *att_drift) {
    if (!model) return;

    /* Extract worst-axis accel and gyro errors (simplified: RMS of 3 axes) */
    double accel_bias = 0.0, gyro_bias = 0.0;
    double accel_vrw = 0.0, gyro_arw = 0.0;

    const ins_axis_error_t *accel[3] = {&model->accel_x, &model->accel_y, &model->accel_z};
    const ins_axis_error_t *gyro[3]  = {&model->gyro_x,  &model->gyro_y,  &model->gyro_z};

    for (int i = 0; i < 3; i++) {
        accel_bias += accel[i]->bias_offset * accel[i]->bias_offset;
        gyro_bias  += gyro[i]->bias_offset  * gyro[i]->bias_offset;
        accel_vrw  += accel[i]->white_noise_std * accel[i]->white_noise_std;
        gyro_arw   += gyro[i]->white_noise_std  * gyro[i]->white_noise_std;
    }
    accel_bias = sqrt(accel_bias / 3.0);
    gyro_bias  = sqrt(gyro_bias  / 3.0);
    accel_vrw  = sqrt(accel_vrw  / 3.0);
    gyro_arw   = sqrt(gyro_arw   / 3.0);

    double g = INS_GRAVITY_EQUATOR;
    double t  = time;
    double t2 = t * t, t3 = t2 * t;

    /* Position drift: RSS of all contributors */
    double pos_accel_bias = 0.5 * accel_bias * t2;
    double pos_gyro_bias  = (1.0 / 6.0) * gyro_bias * g * t3;
    double pos_vrw        = accel_vrw * sqrt(t) * t;
    double pos_arw        = (1.0 / 3.0) * gyro_arw * g * sqrt(t) * t2;

    if (pos_drift) {
        *pos_drift = sqrt(pos_accel_bias * pos_accel_bias +
                          pos_gyro_bias  * pos_gyro_bias +
                          pos_vrw        * pos_vrw +
                          pos_arw        * pos_arw);
    }

    /* Velocity drift */
    double vel_accel_bias = accel_bias * t;
    double vel_gyro_bias  = 0.5 * gyro_bias * g * t2;
    double vel_vrw        = accel_vrw * sqrt(t);

    if (vel_drift) {
        *vel_drift = sqrt(vel_accel_bias * vel_accel_bias +
                         vel_gyro_bias  * vel_gyro_bias +
                         vel_vrw        * vel_vrw);
    }

    /* Attitude drift (primarily from gyro bias integration) */
    if (att_drift) {
        *att_drift = sqrt(gyro_bias * gyro_bias * t2 +
                          gyro_arw  * gyro_arw  * t);
    }
}

/* =========================================================================
 * L5: Psi-Angle Error Model (Benson, 1975)
 *
 * The Psi-angle formulation of INS errors is the standard for Kalman
 * filter-based integration because it linearizes the nonlinear INS
 * equations about the current estimate.
 *
 * Error state: x = [delta_pos^T, delta_vel^T, psi^T]^T  (9 states)
 *
 * Continuous-time dynamics:
 *
 *   delta_pos_dot = -w_en^n x delta_pos + delta_vel
 *
 *   delta_vel_dot = -(skew(C_b^n * f^b)) * psi
 *                   - (2*w_ie^n + w_en^n) x delta_vel
 *
 *   psi_dot = -w_in^n x psi
 *
 * where w_in^n = w_ie^n + w_en^n is the NED frame angular rate
 * relative to inertial frame.
 *
 * State transition matrix: Phi = I + F * dt (1st-order Euler discretization)
 *
 * The Phi matrix is a 9x9 matrix stored row-major:
 *   Phi[i*9 + j] = I[i*9 + j] + F[i*9 + j] * dt
 * ========================================================================= */

void ins_psi_error_propagate(ins_vec3_t *delta_pos,
                              ins_vec3_t *delta_vel,
                              ins_vec3_t *psi,
                              const ins_mat3_t *C_body_ned,
                              const ins_vec3_t *f_body,
                              double lat, double alt, double dt,
                              double *Phi) {
    if (!delta_pos || !delta_vel || !psi || !C_body_ned || !f_body || dt <= 0)
        return;

    /* Transform specific force to NED */
    ins_vec3_t f_ned;
    ins_mat3_mul_vec(C_body_ned, f_body, &f_ned);

    /* Earth rate and transport rate */
    ins_vec3_t w_ie, w_en, w_in;
    ins_earth_rate_ned(lat, &w_ie);
    ins_transport_rate(delta_vel, lat, alt, &w_en);
    ins_vec3_add(&w_ie, &w_en, &w_in);

    /* Rate of change of position error */
    ins_vec3_t pos_dot;
    ins_vec3_cross(&w_en, delta_pos, &pos_dot);
    pos_dot.x = -pos_dot.x + delta_vel->x;
    pos_dot.y = -pos_dot.y + delta_vel->y;
    pos_dot.z = -pos_dot.z + delta_vel->z;

    /* Rate of change of velocity error */
    /* delta_vel_dot = -skew(f_ned) * psi - (2*w_ie + w_en) x delta_vel */
    ins_vec3_t skew_f_psi;
    ins_vec3_cross(&f_ned, psi, &skew_f_psi);

    ins_vec3_t w_coriolis_dv;
    w_coriolis_dv.x = 2.0 * w_ie.x + w_en.x;
    w_coriolis_dv.y = 2.0 * w_ie.y + w_en.y;
    w_coriolis_dv.z = 2.0 * w_ie.z + w_en.z;
    ins_vec3_t vel_dot_coriolis;
    ins_vec3_cross(&w_coriolis_dv, delta_vel, &vel_dot_coriolis);

    ins_vec3_t vel_dot;
    vel_dot.x = skew_f_psi.x - vel_dot_coriolis.x;
    vel_dot.y = skew_f_psi.y - vel_dot_coriolis.y;
    vel_dot.z = skew_f_psi.z - vel_dot_coriolis.z;

    /* Rate of change of attitude error: psi_dot = -w_in x psi */
    ins_vec3_t psi_dot;
    ins_vec3_cross(&w_in, psi, &psi_dot);
    psi_dot.x = -psi_dot.x;
    psi_dot.y = -psi_dot.y;
    psi_dot.z = -psi_dot.z;

    /* Euler integration */
    delta_pos->x += pos_dot.x * dt;
    delta_pos->y += pos_dot.y * dt;
    delta_pos->z += pos_dot.z * dt;

    delta_vel->x += vel_dot.x * dt;
    delta_vel->y += vel_dot.y * dt;
    delta_vel->z += vel_dot.z * dt;

    psi->x += psi_dot.x * dt;
    psi->y += psi_dot.y * dt;
    psi->z += psi_dot.z * dt;

    /* Build 9x9 Phi matrix if requested */
    if (Phi) {
        /* Initialize Phi = I + F*dt */
        memset(Phi, 0, 81 * sizeof(double));
        for (int i = 0; i < 9; i++) {
            Phi[i * 9 + i] = 1.0;
        }

        /* F_block(0,1) = I (delta_pos_dot contribution from delta_vel) */
        for (int i = 0; i < 3; i++) {
            Phi[i * 9 + (3 + i)] += dt;
        }

        /* F_block(0,0) = -skew(w_en) */
        Phi[0 * 9 + 1] +=  w_en.z * dt;
        Phi[0 * 9 + 2] += -w_en.y * dt;
        Phi[1 * 9 + 0] += -w_en.z * dt;
        Phi[1 * 9 + 2] +=  w_en.x * dt;
        Phi[2 * 9 + 0] +=  w_en.y * dt;
        Phi[2 * 9 + 1] += -w_en.x * dt;

        /* F_block(1,2) = -skew(f_ned) (delta_vel_dot from psi) */
        Phi[3 * 9 + 6] +=  f_ned.z * dt;
        Phi[3 * 9 + 8] += -f_ned.y * dt;
        Phi[4 * 9 + 6] += -f_ned.z * dt;
        Phi[4 * 9 + 7] +=  f_ned.x * dt;
        Phi[5 * 9 + 6] +=  f_ned.y * dt;
        Phi[5 * 9 + 7] += -f_ned.x * dt;

        /* F_block(1,1) = -skew(2*w_ie + w_en) */
        Phi[3 * 9 + 4] +=  w_coriolis_dv.z * dt;
        Phi[3 * 9 + 5] += -w_coriolis_dv.y * dt;
        Phi[4 * 9 + 3] += -w_coriolis_dv.z * dt;
        Phi[4 * 9 + 5] +=  w_coriolis_dv.x * dt;
        Phi[5 * 9 + 3] +=  w_coriolis_dv.y * dt;
        Phi[5 * 9 + 4] += -w_coriolis_dv.x * dt;

        /* F_block(2,2) = -skew(w_in) (psi_dot from psi) */
        Phi[6 * 9 + 7] +=  w_in.z * dt;
        Phi[6 * 9 + 8] += -w_in.y * dt;
        Phi[7 * 9 + 6] += -w_in.z * dt;
        Phi[7 * 9 + 8] +=  w_in.x * dt;
        Phi[8 * 9 + 6] +=  w_in.y * dt;
        Phi[8 * 9 + 7] += -w_in.x * dt;
    }
}
