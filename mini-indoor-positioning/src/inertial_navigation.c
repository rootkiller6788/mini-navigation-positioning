/**
 * @file inertial_navigation.c
 * @brief Inertial navigation and dead reckoning algorithms
 *
 * Implements: IMU calibration, strapdown INS mechanization,
 * ZUPT correction, step detection, PDR, Madgwick/Mahony AHRS,
 * quaternion integration (1st-order and RK4).
 */

#include <math.h>
#include <string.h>
#include <stdio.h>
#include "../include/inertial_navigation.h"

/* ============================================================================
 * IMU Calibration
 * ============================================================================ */

void imu_apply_calibration(const imu_sample_t *raw,
                           const imu_calibration_t *cal,
                           imu_sample_t *corrected) {
    if (!raw || !cal || !corrected) return;

    /* Corrected_accel = scale * (misalign_inv * (raw_accel - bias))
     * For simplicity, we assume misalignment is small and apply:
     * corrected = scale * (raw - bias) + misalign_correction */

    for (int i = 0; i < 3; i++) {
        double raw_arr[3] = {raw->accel_x, raw->accel_y, raw->accel_z};
        double gyro_arr[3] = {raw->gyro_x, raw->gyro_y, raw->gyro_z};
        double mag_arr[3] = {raw->mag_x, raw->mag_y, raw->mag_z};

        /* Remove bias, apply scale */
        double a_corr = cal->accel_scale[i] * (raw_arr[i] - cal->accel_bias[i]);
        double g_corr = cal->gyro_scale[i] * (gyro_arr[i] - cal->gyro_bias[i]);
        double m_corr = cal->mag_scale[i] * (mag_arr[i] - cal->mag_bias[i]);

        /* Apply misalignment correction (cross-axis) */
        a_corr += cal->accel_misalign[i][0] * raw_arr[0]
                + cal->accel_misalign[i][1] * raw_arr[1]
                + cal->accel_misalign[i][2] * raw_arr[2];
        g_corr += cal->gyro_misalign[i][0] * gyro_arr[0]
                + cal->gyro_misalign[i][1] * gyro_arr[1]
                + cal->gyro_misalign[i][2] * gyro_arr[2];
        m_corr += cal->mag_misalign[i][0] * mag_arr[0]
                + cal->mag_misalign[i][1] * mag_arr[1]
                + cal->mag_misalign[i][2] * mag_arr[2];

        if (i == 0) {
            corrected->accel_x = a_corr;
            corrected->gyro_x = g_corr;
            corrected->mag_x = m_corr;
        } else if (i == 1) {
            corrected->accel_y = a_corr;
            corrected->gyro_y = g_corr;
            corrected->mag_y = m_corr;
        } else {
            corrected->accel_z = a_corr;
            corrected->gyro_z = g_corr;
            corrected->mag_z = m_corr;
        }
    }
    corrected->timestamp_us = raw->timestamp_us;
    corrected->mag_valid = raw->mag_valid;
}

/* ============================================================================
 * Quaternion Operations
 * ============================================================================ */

void quaternion_multiply(const quaternion_t *q1, const quaternion_t *q2,
                         quaternion_t *result) {
    if (!q1 || !q2 || !result) return;
    /* Hamilton product: q1 ⊗ q2 */
    result->w = q1->w * q2->w - q1->x * q2->x - q1->y * q2->y - q1->z * q2->z;
    result->x = q1->w * q2->x + q1->x * q2->w + q1->y * q2->z - q1->z * q2->y;
    result->y = q1->w * q2->y - q1->x * q2->z + q1->y * q2->w + q1->z * q2->x;
    result->z = q1->w * q2->z + q1->x * q2->y - q1->y * q2->x + q1->z * q2->w;
}

void quaternion_conjugate(const quaternion_t *q, quaternion_t *result) {
    if (!q || !result) return;
    result->w =  q->w;
    result->x = -q->x;
    result->y = -q->y;
    result->z = -q->z;
}

/* ============================================================================
 * L5 - Quaternion Integration (1st-order)
 * ============================================================================ */

void quaternion_integrate_1st_order(quaternion_t *q,
                                   double gx, double gy, double gz, double dt) {
    if (!q) return;
    /* q_{k+1} = q_k ⊗ [1, 0.5*wx*dt, 0.5*wy*dt, 0.5*wz*dt] */
    quaternion_t dq;
    double angle = 0.5 * dt;
    dq.w = 1.0;
    dq.x = gx * angle;
    dq.y = gy * angle;
    dq.z = gz * angle;

    /* Normalize increment to avoid scale errors */
    double norm = sqrt(dq.w*dq.w + dq.x*dq.x + dq.y*dq.y + dq.z*dq.z);
    if (norm > 1e-12) {
        dq.w /= norm;
        dq.x /= norm;
        dq.y /= norm;
        dq.z /= norm;
    }

    quaternion_t temp;
    quaternion_multiply(q, &dq, &temp);
    *q = temp;
    quaternion_normalize(q);
}

/* ============================================================================
 * L5 - Quaternion Integration (Runge-Kutta 4th-order)
 * ============================================================================ */

void quaternion_integrate_rk4(quaternion_t *q,
                              double gx, double gy, double gz, double dt) {
    if (!q) return;

    /* RK4 for quaternion: dq/dt = 0.5 * q ⊗ omega
     * where omega = (0, gx, gy, gz) is the angular velocity quaternion */

    /* omega as quaternion (pure) */
    double wx = gx, wy = gy, wz = gz;

    /* Compute q_dot = 0.5 * q ⊗ omega */
    #define QDOT(qw, qx, qy, qz) \
        do { \
            double _w = -0.5 * (qx*wx + qy*wy + qz*wz); \
            double _x =  0.5 * (qw*wx + qy*wz - qz*wy); \
            double _y =  0.5 * (qw*wy + qz*wx - qx*wz); \
            double _z =  0.5 * (qw*wz + qx*wy - qy*wx); \
            (qw) = _w; (qx) = _x; (qy) = _y; (qz) = _z; \
        } while(0)

    double qw = q->w, qx = q->x, qy = q->y, qz = q->z;

    /* k1 */
    double k1w = qw, k1x = qx, k1y = qy, k1z = qz;
    QDOT(k1w, k1x, k1y, k1z);

    /* k2 = q_dot(q + dt/2 * k1) */
    double k2w = qw + 0.5*dt*k1w;
    double k2x = qx + 0.5*dt*k1x;
    double k2y = qy + 0.5*dt*k1y;
    double k2z = qz + 0.5*dt*k1z;
    QDOT(k2w, k2x, k2y, k2z);

    /* k3 */
    double k3w = qw + 0.5*dt*k2w;
    double k3x = qx + 0.5*dt*k2x;
    double k3y = qy + 0.5*dt*k2y;
    double k3z = qz + 0.5*dt*k2z;
    QDOT(k3w, k3x, k3y, k3z);

    /* k4 */
    double k4w = qw + dt*k3w;
    double k4x = qx + dt*k3x;
    double k4y = qy + dt*k3y;
    double k4z = qz + dt*k3z;
    QDOT(k4w, k4x, k4y, k4z);

    /* q_new = q + dt/6*(k1 + 2*k2 + 2*k3 + k4) */
    q->w = qw + dt/6.0 * (k1w + 2.0*k2w + 2.0*k3w + k4w);
    q->x = qx + dt/6.0 * (k1x + 2.0*k2x + 2.0*k3x + k4x);
    q->y = qy + dt/6.0 * (k1y + 2.0*k2y + 2.0*k3y + k4y);
    q->z = qz + dt/6.0 * (k1z + 2.0*k2z + 2.0*k3z + k4z);

    quaternion_normalize(q);
    #undef QDOT
}

/* ============================================================================
 * L2/L5 - Strapdown INS
 * ============================================================================ */

void ins_init(ins_state_t *ins, const position3d_t *initial_pos,
              double initial_heading, double sample_rate_hz) {
    if (!ins) return;
    memset(ins, 0, sizeof(ins_state_t));
    ins->nav.pos = *initial_pos;
    ins->nav.yaw = initial_heading;
    ins->nav.roll = 0.0;
    ins->nav.pitch = 0.0;
    /* Initialize attitude quaternion from heading */
    euler_to_quaternion(0.0, 0.0, initial_heading, &ins->quat);
    ins->dt = 1.0 / sample_rate_hz;
    ins->pos_cov[0] = 1.0;  /* Initial position uncertainty: 1m */
    ins->pos_cov[1] = 1.0;
    ins->pos_cov[2] = 1.0;
}

void ins_mechanize(ins_state_t *ins, const imu_sample_t *imu) {
    if (!ins || !imu) return;

    /* Time step */
    double dt = ins->dt;
    if (ins->last_time > 0) {
        dt = (imu->timestamp_us - ins->last_time) * 1e-6;
        if (dt <= 0.0 || dt > 1.0) dt = ins->dt;  /* Sanity */
    }
    ins->last_time = imu->timestamp_us;

    /* 1. Attitude update: integrate gyro to update quaternion */
    quaternion_integrate_1st_order(&ins->quat, imu->gyro_x, imu->gyro_y, imu->gyro_z, dt);

    /* 2. Transform specific force from body frame to navigation frame */
    double accel_body[3] = {imu->accel_x, imu->accel_y, imu->accel_z};
    double accel_nav[3];
    quaternion_rotate_vector(accel_body, &ins->quat, accel_nav);

    /* 3. Velocity update: subtract gravity from vertical channel (ENU: z is up)
     * a_nav = f_nav - [0, 0, g] where g = 9.80665 m/s^2 downward
     * In ENU: gravity is -g in z direction */
    accel_nav[2] -= IMU_GRAVITY_MSS;

    /* 4. Position update: integrate velocity */
    ins->nav.vel.vx += accel_nav[0] * dt;
    ins->nav.vel.vy += accel_nav[1] * dt;
    ins->nav.vel.vz += accel_nav[2] * dt;

    ins->nav.pos.x += ins->nav.vel.vx * dt;
    ins->nav.pos.y += ins->nav.vel.vy * dt;
    ins->nav.pos.z += ins->nav.vel.vz * dt;

    /* Update Euler angles for convenience */
    quaternion_to_euler(&ins->quat, &ins->nav.roll, &ins->nav.pitch, &ins->nav.yaw);

    /* Propagate position uncertainty (simplified: grows quadratically with dt) */
    double drift_rate = 1e-4;  /* 0.1 mm/s drift */
    ins->pos_cov[0] += drift_rate * dt;
    ins->pos_cov[1] += drift_rate * dt;
    ins->pos_cov[2] += drift_rate * dt;
}

/* ============================================================================
 * L5 - Zero Velocity Update (ZUPT)
 * ============================================================================ */

int detect_zero_velocity(const imu_sample_t *imu,
                         double accel_thresh, double gyro_thresh) {
    if (!imu) return 0;
    double accel_mag = sqrt(imu->accel_x * imu->accel_x +
                            imu->accel_y * imu->accel_y +
                            imu->accel_z * imu->accel_z);
    double gyro_mag = sqrt(imu->gyro_x * imu->gyro_x +
                           imu->gyro_y * imu->gyro_y +
                           imu->gyro_z * imu->gyro_z);

    /* Acceleration magnitude should be close to gravity (9.81 ± threshold) */
    /* Gyroscope magnitude should be near zero */
    if (fabs(accel_mag - IMU_GRAVITY_MSS) < accel_thresh && gyro_mag < gyro_thresh) {
        return 1;
    }
    return 0;
}

void ins_apply_zupt(ins_state_t *ins, int zupt_detected) {
    if (!ins) return;
    if (!zupt_detected) return;

    /* When zero velocity is detected, set velocity to zero */
    ins->nav.vel.vx = 0.0;
    ins->nav.vel.vy = 0.0;
    ins->nav.vel.vz = 0.0;

    /* Use the known-zero-velocity to estimate and correct gyro bias
     * (simplified correction: damp position uncertainty) */
    ins->pos_cov[0] *= 0.5;
    ins->pos_cov[1] *= 0.5;
    ins->pos_cov[2] *= 0.5;
}

/* ============================================================================
 * L5 - Madgwick AHRS Algorithm
 * ============================================================================ */

void madgwick_init(madgwick_ahrs_t *ahrs, double beta, double sample_period) {
    if (!ahrs) return;
    /* Initialize quaternion to identity (no rotation) */
    ahrs->q.w = 1.0; ahrs->q.x = 0.0; ahrs->q.y = 0.0; ahrs->q.z = 0.0;
    ahrs->beta = beta;
    ahrs->sample_period = sample_period;
}

void madgwick_update_imu(madgwick_ahrs_t *ahrs,
                         double gx, double gy, double gz,
                         double ax, double ay, double az) {
    if (!ahrs) return;

    double q1 = ahrs->q.w, q2 = ahrs->q.x, q3 = ahrs->q.y, q4 = ahrs->q.z;
    double beta = ahrs->beta;

    /* Normalize accelerometer measurement */
    double recip_norm;
    double norm_a = sqrt(ax*ax + ay*ay + az*az);
    if (norm_a < 1e-12) return;  /* Avoid division by zero */
    recip_norm = 1.0 / norm_a;
    ax *= recip_norm;
    ay *= recip_norm;
    az *= recip_norm;

    /* Gradient descent algorithm corrective step
     * Objective function: f = [2*(q2*q4 - q1*q3) - ax, ...]
     * Jacobian: J = ∂f/∂q */

    /* Compute the gradient descent step (simplified from Madgwick paper) */
    double _2q1 = 2.0 * q1;
    double _2q2 = 2.0 * q2;
    double _2q3 = 2.0 * q3;
    double _2q4 = 2.0 * q4;

    /* f = estimated gravity direction in body frame */
    double f1 = _2q2 * q4 - _2q1 * q3 - ax;
    double f2 = _2q1 * q2 + _2q3 * q4 - ay;
    double f3 = 1.0 - _2q2 * q2 - _2q3 * q3 - az;

    /* Jacobian J^T * f (simplified: gradient = J^T * f) */
    double j11 = -_2q3;
    double j12 =  _2q4;
    double j13 = -_2q1;
    double j14 =  _2q2;
    double j21 =  _2q2;
    double j22 =  _2q1;
    double j23 =  _2q4;
    double j24 =  _2q3;
    double j31 =  _2q1;
    double j32 = -_2q2;
    double j33 = -_2q3;
    double j34 =  _2q4;

    /* Gradient: ∇f = J^T * f */
    double s1 = j11*f1 + j21*f2 + j31*f3;
    double s2 = j12*f1 + j22*f2 + j32*f3;
    double s3 = j13*f1 + j23*f2 + j33*f3;
    double s4 = j14*f1 + j24*f2 + j34*f3;

    /* Normalize step direction */
    double norm_s = sqrt(s1*s1 + s2*s2 + s3*s3 + s4*s4);
    if (norm_s < 1e-12) {
        /* Skip correction if gradient is zero */
        s1 = s2 = s3 = s4 = 0.0;
    } else {
        s1 /= norm_s; s2 /= norm_s; s3 /= norm_s; s4 /= norm_s;
    }

    /* Compute rate of change of quaternion: q_dot = 0.5*q⊗omega - beta*gradient */
    double q_dot1 = 0.5 * (-q2*gx - q3*gy - q4*gz) - beta * s1;
    double q_dot2 = 0.5 * ( q1*gx + q3*gz - q4*gy) - beta * s2;
    double q_dot3 = 0.5 * ( q1*gy - q2*gz + q4*gx) - beta * s3;
    double q_dot4 = 0.5 * ( q1*gz + q2*gy - q3*gx) - beta * s4;

    /* Integrate */
    q1 += q_dot1 * ahrs->sample_period;
    q2 += q_dot2 * ahrs->sample_period;
    q3 += q_dot3 * ahrs->sample_period;
    q4 += q_dot4 * ahrs->sample_period;

    /* Normalize */
    double norm_q = sqrt(q1*q1 + q2*q2 + q3*q3 + q4*q4);
    if (norm_q < 1e-12) return;
    ahrs->q.w = q1 / norm_q;
    ahrs->q.x = q2 / norm_q;
    ahrs->q.y = q3 / norm_q;
    ahrs->q.z = q4 / norm_q;
}

/* ============================================================================
 * L5 - Madgwick AHRS with Magnetometer (MARG mode)
 * ============================================================================ */

void madgwick_update_marg(madgwick_ahrs_t *ahrs,
                          double gx, double gy, double gz,
                          double ax, double ay, double az,
                          double mx, double my, double mz) {
    if (!ahrs) return;

    double q1 = ahrs->q.w, q2 = ahrs->q.x, q3 = ahrs->q.y, q4 = ahrs->q.z;
    double beta = ahrs->beta;

    /* Normalize accelerometer */
    double norm_a = sqrt(ax*ax + ay*ay + az*az);
    if (norm_a < 1e-12) norm_a = 1e-12;
    ax /= norm_a; ay /= norm_a; az /= norm_a;

    /* Normalize magnetometer */
    double norm_m = sqrt(mx*mx + my*my + mz*mz);
    if (norm_m < 1e-12) norm_m = 1e-12;
    mx /= norm_m; my /= norm_m; mz /= norm_m;

    /* Reference direction of Earth's magnetic field (in navigation frame)
     * h = q ⊗ (0, mx, my, mz) ⊗ q* */
    double _2q1 = 2.0 * q1, _2q2 = 2.0 * q2, _2q3 = 2.0 * q3, _2q4 = 2.0 * q4;
    double _2q1mx = 2.0 * q1 * mx;
    double _2q1my = 2.0 * q1 * my;
    double _2q1mz = 2.0 * q1 * mz;
    double _2q2mx = 2.0 * q2 * mx;
    double hx = mx * q1*q1 - _2q1my * q4 + _2q1mz * q3 + mx * q2*q2
              + _2q2 * my * q3 + _2q2 * mz * q4 - mx * q3*q3 - mx * q4*q4;
    double hy = _2q1mx * q4 + my * q1*q1 - _2q1mz * q2 + _2q2mx * q3
              - my * q2*q2 + my * q3*q3 + _2q3 * mz * q4 - my * q4*q4;

    /* Reference magnetic field direction (North component only) */
    double _2bx = sqrt(hx*hx + hy*hy);
    double _2bz = -_2q1mx * q3 + _2q1my * q2 + mz * q1*q1 + _2q2mx * q4
                - mz * q2*q2 + _2q3 * my * q4 - mz * q3*q3 + mz * q4*q4;

    /* Gradient descent for accelerometer part (same as IMU) */
    double f1 = _2q2*q4 - _2q1*q3 - ax;
    double f2 = _2q1*q2 + _2q3*q4 - ay;
    double f3 = 1.0 - _2q2*q2 - _2q3*q3 - az;

    /* Magnetometer part */
    double f4 = _2bx*(0.5 - q3*q3 - q4*q4) + _2bz*(q2*q4 - q1*q3) - mx;
    double f5 = _2bx*(q2*q3 - q1*q4) + _2bz*(q1*q2 + q3*q4) - my;
    double f6 = _2bx*(q1*q3 + q2*q4) + _2bz*(0.5 - q2*q2 - q3*q3) - mz;

    /* Compute gradient step (combined) */
    double j11 = -_2q3, j12 = _2q4, j13 = -_2q1, j14 = _2q2;
    double j21 = _2q2,  j22 = _2q1, j23 = _2q4, j24 = _2q3;

    /* Mag Jacobian elements */
    double j41 = -_2bz*q3,                   j42 = _2bz*q4;
    double j43 = -4.0*_2bx*q3 - _2bz*q1,     j44 = -4.0*_2bx*q4 + _2bz*q2;
    double j51 = -_2bx*q4 + _2bz*q2,          j52 = _2bx*q3 + _2bz*q1;
    double j53 = _2bx*q2 + _2bz*q4,           j54 = -_2bx*q1 + _2bz*q3;
    double j61 = _2bx*q3,                     j62 = _2bx*q4 - 4.0*_2bz*q2;
    double j63 = _2bx*q1 - 4.0*_2bz*q3,       j64 = _2bx*q2;

    double s1 = j11*f1 + j21*f2 + j41*f4 + j51*f5 + j61*f6;
    double s2 = j12*f1 + j22*f2 + j42*f4 + j52*f5 + j62*f6;
    double s3 = j13*f1 + j23*f2 + j43*f4 + j53*f5 + j63*f6;
    double s4 = j14*f1 + j24*f2 + j44*f4 + j54*f5 + j64*f6;

    /* Normalize step */
    double norm_s = sqrt(s1*s1 + s2*s2 + s3*s3 + s4*s4);
    if (norm_s < 1e-12) norm_s = 1e-12;
    s1 /= norm_s; s2 /= norm_s; s3 /= norm_s; s4 /= norm_s;

    /* Gyro contribution */
    double q_dot1 = 0.5*(-q2*gx - q3*gy - q4*gz) - beta*s1;
    double q_dot2 = 0.5*( q1*gx + q3*gz - q4*gy) - beta*s2;
    double q_dot3 = 0.5*( q1*gy - q2*gz + q4*gx) - beta*s3;
    double q_dot4 = 0.5*( q1*gz + q2*gy - q3*gx) - beta*s4;

    q1 += q_dot1 * ahrs->sample_period;
    q2 += q_dot2 * ahrs->sample_period;
    q3 += q_dot3 * ahrs->sample_period;
    q4 += q_dot4 * ahrs->sample_period;

    double norm_q = sqrt(q1*q1 + q2*q2 + q3*q3 + q4*q4);
    if (norm_q < 1e-12) norm_q = 1e-12;
    ahrs->q.w = q1/norm_q; ahrs->q.x = q2/norm_q;
    ahrs->q.y = q3/norm_q; ahrs->q.z = q4/norm_q;
}

/* ============================================================================
 * L5 - Mahony AHRS Algorithm
 * ============================================================================ */

void mahony_init(mahony_ahrs_t *ahrs, double kp, double ki, double sample_period) {
    if (!ahrs) return;
    ahrs->q.w = 1.0; ahrs->q.x = 0.0; ahrs->q.y = 0.0; ahrs->q.z = 0.0;
    ahrs->kp = kp;
    ahrs->ki = ki;
    ahrs->integral_fb[0] = 0.0;
    ahrs->integral_fb[1] = 0.0;
    ahrs->integral_fb[2] = 0.0;
    ahrs->sample_period = sample_period;
}

void mahony_update_imu(mahony_ahrs_t *ahrs,
                       double gx, double gy, double gz,
                       double ax, double ay, double az) {
    if (!ahrs) return;

    /* Normalize accelerometer */
    double norm_a = sqrt(ax*ax + ay*ay + az*az);
    if (norm_a < 1e-12) norm_a = 1e-12;
    ax /= norm_a; ay /= norm_a; az /= norm_a;

    /* Estimated direction of gravity (from current orientation) */
    double q1 = ahrs->q.w, q2 = ahrs->q.x, q3 = ahrs->q.y, q4 = ahrs->q.z;
    double vx = 2.0 * (q2*q4 - q1*q3);
    double vy = 2.0 * (q1*q2 + q3*q4);
    double vz = q1*q1 - q2*q2 - q3*q3 + q4*q4;

    /* Error is cross product between estimated and measured gravity direction */
    double ex = (ay*vz - az*vy);
    double ey = (az*vx - ax*vz);
    double ez = (ax*vy - ay*vx);

    /* Integral feedback */
    if (ahrs->ki > 0.0) {
        ahrs->integral_fb[0] += ahrs->ki * ex * ahrs->sample_period;
        ahrs->integral_fb[1] += ahrs->ki * ey * ahrs->sample_period;
        ahrs->integral_fb[2] += ahrs->ki * ez * ahrs->sample_period;
    } else {
        ahrs->integral_fb[0] = 0.0;
        ahrs->integral_fb[1] = 0.0;
        ahrs->integral_fb[2] = 0.0;
    }

    /* Corrected angular rate = gyro + kp*error + integral_feedback */
    gx += ahrs->kp * ex + ahrs->integral_fb[0];
    gy += ahrs->kp * ey + ahrs->integral_fb[1];
    gz += ahrs->kp * ez + ahrs->integral_fb[2];

    /* Integrate quaternion using corrected gyro rates */
    quaternion_integrate_1st_order(&ahrs->q, gx, gy, gz, ahrs->sample_period);
}

/* ============================================================================
 * L6 - Pedestrian Dead Reckoning (PDR)
 * ============================================================================ */

void pdr_init(pdr_state_t *pdr, position2d_t initial_pos,
              double initial_heading, double stride_length) {
    if (!pdr) return;
    memset(pdr, 0, sizeof(pdr_state_t));
    pdr->position = initial_pos;
    pdr->heading = initial_heading;
    pdr->step_det.step_threshold = 0.5;  /* 0.5 m/s^2 above gravity */
    pdr->step_det.step_timeout_ms = 250.0;  /* Max 4 steps per second */
    pdr->step_det.stride_length = stride_length;
}

int pdr_process_accel(pdr_state_t *pdr, double accel_magnitude,
                      double heading, double timestamp_ms) {
    if (!pdr) return 0;

    /* Step detection: peak detection on filtered acceleration magnitude.
     * Detects when magnitude crosses threshold from above to below
     * (negative slope zero-crossing of bandpass-filtered signal,
     *  simplified here as threshold crossing with timeout). */

    step_detector_t *sd = &pdr->step_det;

    if (timestamp_ms - sd->last_step_time < sd->step_timeout_ms) {
        /* Too soon after last step */
        sd->accel_mag_prev = accel_magnitude;
        return 0;
    }

    /* Detect if acceleration magnitude exceeds threshold */
    if (accel_magnitude > (IMU_GRAVITY_MSS + sd->step_threshold)
        && sd->accel_mag_prev <= (IMU_GRAVITY_MSS + sd->step_threshold)) {
        /* Rising edge detected → step event */
        sd->last_step_time = timestamp_ms;
        sd->step_count++;

        /* Update position using current heading and stride length */
        pdr->position.x += sd->stride_length * sin(heading);
        pdr->position.y += sd->stride_length * cos(heading);
        pdr->heading = heading;
        pdr->total_distance += sd->stride_length;
        sd->accel_mag_prev = accel_magnitude;
        return 1;
    }

    sd->accel_mag_prev = accel_magnitude;
    return 0;
}

/* ============================================================================
 * Stride Length Estimation Models
 * ============================================================================ */

double pdr_stride_length_weinberg(double accel_max, double accel_min, double K) {
    /* Weinberg model: L = K * (a_max - a_min)^(1/4)
     * a_max, a_min should be in m/s^2 after removing gravity */
    double diff = accel_max - accel_min;
    if (diff < 0.01) diff = 0.01;
    return K * pow(diff, 0.25);
}

double pdr_stride_length_kim(double avg_abs_accel, double K) {
    /* Kim model: L = K * (avg_abs_accel)^(1/3) */
    if (avg_abs_accel < 0.01) avg_abs_accel = 0.01;
    return K * cbrt(avg_abs_accel);
}
