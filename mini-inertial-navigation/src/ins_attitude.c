/**
 * @file    ins_attitude.c
 * @brief   Attitude representation: quaternions, Euler angles, DCM
 *
 * Knowledge Coverage:
 *   L1 (Definitions): Quaternion, Euler angle, DCM types
 *   L2 (Core Concepts): Quaternion algebra, rotation composition
 *   L3 (Math Structures): Quaternion kinematics, SO(3) Lie algebra
 *   L4 (Fundamental Laws): Quaternion propagation preserves unit norm
 *   L5 (Algorithms): Quaternion integration (Euler, exact, Bortz coning)
 *
 * Reference:
 *   Kuipers (1999), "Quaternions and Rotation Sequences", Princeton.
 *   Markley & Crassidis (2014), "Fundamentals of Spacecraft Attitude
 *     Determination and Control", Springer.
 *   Bortz (1971), IEEE Trans. Aerosp. Electron. Syst., 7(1): 61-66.
 *   Savage (1998), J. Guid. Control Dyn., 21(1): 19-28.
 *
 * Course Mapping:
 *   MIT 6.832 - Underactuated Robotics (3D rotations)
 *   Stanford AA279 - Spacecraft Attitude Dynamics
 *   Berkeley EE C128 - Mechatronics (IMU orientation)
 *   ETH 151-0854 - Trajectory Generation (quaternion splines)
 */

#include "ins_attitude.h"
#include "ins_core.h"
#include <string.h>
#include <math.h>
#include <stdio.h>

/* =========================================================================
 * L2: Quaternion Core Operations
 * ========================================================================= */

void ins_quat_identity(ins_quat_t *q) {
    if (!q) return;
    q->w = 1.0;
    q->x = q->y = q->z = 0.0;
}

void ins_quat_copy(const ins_quat_t *src, ins_quat_t *dst) {
    if (!src || !dst) return;
    dst->w = src->w;
    dst->x = src->x;
    dst->y = src->y;
    dst->z = src->z;
}

double ins_quat_norm_sq(const ins_quat_t *q) {
    if (!q) return 0.0;
    return q->w * q->w + q->x * q->x + q->y * q->y + q->z * q->z;
}

double ins_quat_norm(const ins_quat_t *q) {
    if (!q) return 0.0;
    return sqrt(ins_quat_norm_sq(q));
}

int ins_quat_normalize(ins_quat_t *q) {
    if (!q) return -1;
    double n_sq = ins_quat_norm_sq(q);
    if (n_sq < 1e-15) return -1;
    double inv_n = 1.0 / sqrt(n_sq);
    q->w *= inv_n;
    q->x *= inv_n;
    q->y *= inv_n;
    q->z *= inv_n;
    return 0;
}

void ins_quat_mul(const ins_quat_t *p, const ins_quat_t *q, ins_quat_t *r) {
    if (!p || !q || !r) return;
    r->w = p->w * q->w - p->x * q->x - p->y * q->y - p->z * q->z;
    r->x = p->w * q->x + p->x * q->w + p->y * q->z - p->z * q->y;
    r->y = p->w * q->y - p->x * q->z + p->y * q->w + p->z * q->x;
    r->z = p->w * q->z + p->x * q->y - p->y * q->x + p->z * q->w;
}

void ins_quat_conjugate(const ins_quat_t *q, ins_quat_t *qc) {
    if (!q || !qc) return;
    qc->w = q->w;
    qc->x = -q->x;
    qc->y = -q->y;
    qc->z = -q->z;
}

void ins_quat_inverse(const ins_quat_t *q, ins_quat_t *qi) {
    if (!q || !qi) return;
    ins_quat_conjugate(q, qi);
    double n_sq = ins_quat_norm_sq(q);
    if (n_sq > 1e-15) {
        double inv_n_sq = 1.0 / n_sq;
        qi->w *= inv_n_sq;
        qi->x *= inv_n_sq;
        qi->y *= inv_n_sq;
        qi->z *= inv_n_sq;
    }
}

/* =========================================================================
 * L2: Quaternion <-> DCM Conversion
 *
 * For unit quaternion q = (w, x, y, z), the DCM is:
 *
 *     [ 1-2(y^2+z^2)   2(xy-wz)       2(xz+wy)     ]
 * C = [ 2(xy+wz)       1-2(x^2+z^2)   2(yz-wx)     ]
 *     [ 2(xz-wy)       2(yz+wx)       1-2(x^2+y^2) ]
 *
 * ========================================================================= */

void ins_quat_to_dcm(const ins_quat_t *q, ins_mat3_t *C) {
    if (!q || !C) return;
    double w = q->w, x = q->x, y = q->y, z = q->z;
    double xx = x * x, yy = y * y, zz = z * z;
    double xy = x * y, xz = x * z, yz = y * z;
    double wx = w * x, wy = w * y, wz = w * z;

    C->m[0] = 1.0 - 2.0 * (yy + zz);
    C->m[1] = 2.0 * (xy - wz);
    C->m[2] = 2.0 * (xz + wy);

    C->m[3] = 2.0 * (xy + wz);
    C->m[4] = 1.0 - 2.0 * (xx + zz);
    C->m[5] = 2.0 * (yz - wx);

    C->m[6] = 2.0 * (xz - wy);
    C->m[7] = 2.0 * (yz + wx);
    C->m[8] = 1.0 - 2.0 * (xx + yy);
}

void ins_dcm_to_quat(const ins_mat3_t *C, ins_quat_t *q) {
    if (!C || !q) return;

    /* Shepperd (1978) method: select formula based on largest trace element */
    double trace = C->m[0] + C->m[4] + C->m[8];
    double w, x, y, z;

    if (trace > 0.0) {
        double s = sqrt(trace + 1.0);
        w = s * 0.5;
        s = 0.5 / s;
        x = (C->m[7] - C->m[5]) * s;
        y = (C->m[2] - C->m[6]) * s;
        z = (C->m[3] - C->m[1]) * s;
    } else {
        int i = 0;
        if (C->m[4] > C->m[0]) i = 1;
        if (C->m[8] > C->m[4] && C->m[8] > C->m[0]) i = 2;

        switch (i) {
        case 0: {
            double s = sqrt(C->m[0] - C->m[4] - C->m[8] + 1.0);
            x = s * 0.5;
            s = 0.5 / s;
            w = (C->m[7] - C->m[5]) * s;
            y = (C->m[1] + C->m[3]) * s;
            z = (C->m[2] + C->m[6]) * s;
            break;
        }
        case 1: {
            double s = sqrt(C->m[4] - C->m[0] - C->m[8] + 1.0);
            y = s * 0.5;
            s = 0.5 / s;
            w = (C->m[2] - C->m[6]) * s;
            x = (C->m[1] + C->m[3]) * s;
            z = (C->m[5] + C->m[7]) * s;
            break;
        }
        default: {
            double s = sqrt(C->m[8] - C->m[0] - C->m[4] + 1.0);
            z = s * 0.5;
            s = 0.5 / s;
            w = (C->m[3] - C->m[1]) * s;
            x = (C->m[2] + C->m[6]) * s;
            y = (C->m[5] + C->m[7]) * s;
            break;
        }
        }
    }
    q->w = w; q->x = x; q->y = y; q->z = z;
}

/* =========================================================================
 * L2: Quaternion <-> Euler Angles (3-2-1 Yaw-Pitch-Roll)
 *
 * Roll (phi), Pitch (theta), Yaw (psi) using aerospace NED:
 *
 *   DCM = R_x(roll) * R_y(pitch) * R_z(yaw)
 *
 * Quaternion composition:
 *   q = q_z(yaw) * q_y(pitch) * q_x(roll)
 *
 *   q_z(psi) = [cos(psi/2), 0, 0, sin(psi/2)]
 *   q_y(theta) = [cos(theta/2), 0, sin(theta/2), 0]
 *   q_x(phi) = [cos(phi/2), sin(phi/2), 0, 0]
 * ========================================================================= */

void ins_quat_to_euler(const ins_quat_t *q, ins_euler_t *e) {
    if (!q || !e) return;
    double w = q->w, x = q->x, y = q->y, z = q->z;

    /* Roll (phi) */
    double sinr_cosp = 2.0 * (w * x + y * z);
    double cosr_cosp = 1.0 - 2.0 * (x * x + y * y);
    e->roll = atan2(sinr_cosp, cosr_cosp);

    /* Pitch (theta) */
    double sinp = 2.0 * (w * y - z * x);
    if (fabs(sinp) >= 1.0) {
        e->pitch = copysign(M_PI / 2.0, sinp);
    } else {
        e->pitch = asin(sinp);
    }

    /* Yaw (psi) */
    double siny_cosp = 2.0 * (w * z + x * y);
    double cosy_cosp = 1.0 - 2.0 * (y * y + z * z);
    e->yaw = atan2(siny_cosp, cosy_cosp);
}

void ins_euler_to_quat(const ins_euler_t *e, ins_quat_t *q) {
    if (!e || !q) return;
    double half_r = e->roll * 0.5;
    double half_p = e->pitch * 0.5;
    double half_y = e->yaw * 0.5;

    double cr = cos(half_r), sr = sin(half_r);
    double cp = cos(half_p), sp = sin(half_p);
    double cy = cos(half_y), sy = sin(half_y);

    /* q = q_z(yaw) * q_y(pitch) * q_x(roll) */
    /* q_z: [cy, 0, 0, sy] */
    /* q_y: [cp, 0, sp, 0] */
    /* q_x: [cr, sr, 0, 0] */

    /* q_y * q_x first */
    double qyx_w = cp * cr;
    double qyx_x = cp * sr;
    double qyx_y = sp * cr;
    double qyx_z = -sp * sr;

    /* q_z * (q_y * q_x) */
    q->w = cy * qyx_w - sy * qyx_z;
    q->x = cy * qyx_x - sy * qyx_y;
    q->y = cy * qyx_y + sy * qyx_x;
    q->z = cy * qyx_z + sy * qyx_w;
}

/* =========================================================================
 * L2: Euler <-> DCM
 * ========================================================================= */

void ins_euler_to_dcm(const ins_euler_t *e, ins_mat3_t *C) {
    if (!e || !C) return;
    double cr = cos(e->roll), sr = sin(e->roll);
    double cp = cos(e->pitch), sp = sin(e->pitch);
    double cy = cos(e->yaw), sy = sin(e->yaw);

    /* C = R_x(roll) * R_y(pitch) * R_z(yaw) — body-to-NED */
    /* Row 0 */
    C->m[0] = cp * cy;
    C->m[1] = cp * sy;
    C->m[2] = -sp;
    /* Row 1 */
    C->m[3] = sr * sp * cy - cr * sy;
    C->m[4] = sr * sp * sy + cr * cy;
    C->m[5] = sr * cp;
    /* Row 2 */
    C->m[6] = cr * sp * cy + sr * sy;
    C->m[7] = cr * sp * sy - sr * cy;
    C->m[8] = cr * cp;
}

void ins_dcm_to_euler(const ins_mat3_t *C, ins_euler_t *e) {
    if (!C || !e) return;
    /* Pitch: from C[0][2] = -sin(pitch) */
    if (fabs(C->m[2]) >= 1.0) {
        e->pitch = copysign(M_PI / 2.0, -C->m[2]);
        e->roll = 0.0;  /* Gimbal lock: roll and yaw degenerate */
        e->yaw = atan2(-C->m[3], C->m[0]);
    } else {
        e->pitch = asin(-C->m[2]);
        e->roll = atan2(C->m[5], C->m[8]);
        e->yaw = atan2(C->m[1], C->m[0]);
    }
}

/* =========================================================================
 * L2: Vector Rotation by Quaternion
 *
 * v_rot = q * v * q^(-1)
 *
 * Efficient implementation using Rodrigues formula / triple product:
 *   v_rot = v + 2*r x (r x v + w*v)
 * where q = [w, r] with r = [x, y, z]
 *
 * This avoids constructing the full DCM and is more efficient
 * (30 multiplications vs 45 for DCM approach).
 * ========================================================================= */

void ins_quat_rotate_vector(const ins_quat_t *q, const ins_vec3_t *v,
                             ins_vec3_t *v_rot) {
    if (!q || !v || !v_rot) return;
    double w = q->w;
    ins_vec3_t r; r.x = q->x; r.y = q->y; r.z = q->z;

    /* t = 2 * r x v */
    ins_vec3_t t;
    ins_vec3_cross(&r, v, &t);
    ins_vec3_scale(&t, 2.0, &t);

    /* v_rot = v + w*t + r x t */
    ins_vec3_t r_cross_t;
    ins_vec3_cross(&r, &t, &r_cross_t);

    v_rot->x = v->x + w * t.x + r_cross_t.x;
    v_rot->y = v->y + w * t.y + r_cross_t.y;
    v_rot->z = v->z + w * t.z + r_cross_t.z;
}

/* =========================================================================
 * L3: Quaternion Kinematics
 *
 * The differential equation for quaternion propagation from
 * angular velocity omega = [wx, wy, wz]^T is:
 *
 *   dq/dt = 0.5 * q * Omega
 *
 * where Omega is the pure quaternion [0, wx, wy, wz].
 *
 * Theorem: Unit norm preservation
 *   d/dt(|q|^2) = 2q·(dq/dt) = q·(q*Omega) = |q|^2 * Omega.w = 0
 *   since Omega.w = 0 for pure vector quaternion.
 *   Hence |q(t)| = const for all t.
 * ========================================================================= */

void ins_quat_kinematics(const ins_quat_t *q, const ins_vec3_t *omega,
                          ins_quat_t *dq_dt) {
    if (!q || !omega || !dq_dt) return;
    double wx = omega->x, wy = omega->y, wz = omega->z;

    dq_dt->w = -0.5 * (wx * q->x + wy * q->y + wz * q->z);
    dq_dt->x =  0.5 * (wx * q->w + wz * q->y - wy * q->z);
    dq_dt->y =  0.5 * (wy * q->w - wz * q->x + wx * q->z);
    dq_dt->z =  0.5 * (wz * q->w + wy * q->x - wx * q->y);
}

/* =========================================================================
 * L5: Quaternion Integration Methods
 * ========================================================================= */

void ins_quat_update_euler(ins_quat_t *q, const ins_vec3_t *omega, double dt) {
    if (!q || !omega || dt <= 0) return;

    ins_quat_t dq, q_new;
    ins_quat_kinematics(q, omega, &dq);

    q_new.w = q->w + dq.w * dt;
    q_new.x = q->x + dq.x * dt;
    q_new.y = q->y + dq.y * dt;
    q_new.z = q->z + dq.z * dt;

    /* Renormalize: |q|^2 after Euler step = 1 + |omega|^2*dt^2/4 */
    ins_quat_normalize(&q_new);
    ins_quat_copy(&q_new, q);
}

void ins_quat_update_exact(ins_quat_t *q, const ins_vec3_t *omega, double dt) {
    if (!q || !omega || dt <= 0) return;

    double wx = omega->x, wy = omega->y, wz = omega->z;
    double w_norm = sqrt(wx * wx + wy * wy + wz * wz);
    double half_angle = w_norm * dt * 0.5;

    ins_quat_t q_delta;
    /* Taylor expansion for small angles to avoid sin(x)/x singularity */
    if (half_angle < 1e-7) {
        double s = dt * 0.5 * (1.0 - half_angle * half_angle / 6.0);
        q_delta.w = 1.0 - half_angle * half_angle * 0.5;
        q_delta.x = wx * s;
        q_delta.y = wy * s;
        q_delta.z = wz * s;
    } else {
        double sin_h = sin(half_angle);
        double cos_h = cos(half_angle);
        double scale = sin_h / w_norm;
        q_delta.w = cos_h;
        q_delta.x = wx * scale;
        q_delta.y = wy * scale;
        q_delta.z = wz * scale;
    }

    /* q_new = q * q_delta (rotation composition: first delta, then q) */
    ins_quat_t q_new;
    ins_quat_mul(q, &q_delta, &q_new);
    ins_quat_normalize(&q_new);
    ins_quat_copy(&q_new, q);
}

void ins_quat_update_coning(ins_quat_t *q, const ins_vec3_t *omega_prev,
                             const ins_vec3_t *omega_curr, double dt) {
    if (!q || !omega_prev || !omega_curr || dt <= 0) return;

    /* Bortz (1971) rotation vector with coning compensation */
    /* phi = omega_curr * dt + (1/12) * (omega_prev x omega_curr) * dt^2 */

    ins_vec3_t cross_term;
    ins_vec3_cross(omega_prev, omega_curr, &cross_term);

    double dt2 = dt * dt;
    ins_vec3_t phi;
    phi.x = omega_curr->x * dt + cross_term.x * dt2 / 12.0;
    phi.y = omega_curr->y * dt + cross_term.y * dt2 / 12.0;
    phi.z = omega_curr->z * dt + cross_term.z * dt2 / 12.0;

    /* q_delta = rotation vector to quaternion */
    double phi_norm = ins_vec3_norm(&phi);
    double half_phi = phi_norm * 0.5;

    ins_quat_t q_delta;
    if (phi_norm < 1e-12) {
        ins_quat_identity(&q_delta);
    } else if (half_phi < 1e-7) {
        double s = 0.5 * (1.0 - half_phi * half_phi / 6.0);
        q_delta.w = 1.0 - half_phi * half_phi * 0.5;
        q_delta.x = phi.x * s;
        q_delta.y = phi.y * s;
        q_delta.z = phi.z * s;
    } else {
        double sin_h = sin(half_phi);
        double cos_h = cos(half_phi);
        double scale = sin_h / phi_norm;
        q_delta.w = cos_h;
        q_delta.x = phi.x * scale;
        q_delta.y = phi.y * scale;
        q_delta.z = phi.z * scale;
    }

    ins_quat_t q_new;
    ins_quat_mul(q, &q_delta, &q_new);
    ins_quat_normalize(&q_new);
    ins_quat_copy(&q_new, q);
}

/* =========================================================================
 * L5: Sculling Compensation
 *
 * When a vehicle undergoes simultaneous linear and angular vibration,
 * a net velocity error accumulates ("sculling" effect).
 *
 * Two-sample algorithm (Savage, 1998):
 *   delta_v_scull = (2/3) * (dalpha1 x dv2 + dv1 x dalpha2)
 * ========================================================================= */

void ins_sculling_compensation(const ins_vec3_t *delta_alpha1,
                                const ins_vec3_t *delta_alpha2,
                                const ins_vec3_t *delta_v1,
                                const ins_vec3_t *delta_v2,
                                ins_vec3_t *delta_v_scull) {
    if (!delta_alpha1 || !delta_alpha2 || !delta_v1 || !delta_v2 || !delta_v_scull)
        return;

    ins_vec3_t term1, term2, temp;
    ins_vec3_cross(delta_alpha1, delta_v2, &term1);
    ins_vec3_cross(delta_v1, delta_alpha2, &term2);
    ins_vec3_add(&term1, &term2, &temp);
    ins_vec3_scale(&temp, 2.0 / 3.0, delta_v_scull);
}

/* =========================================================================
 * L3: Angle Wrapping
 * ========================================================================= */

double ins_angle_wrap(double angle) {
    while (angle > M_PI)  angle -= 2.0 * M_PI;
    while (angle < -M_PI) angle += 2.0 * M_PI;
    return angle;
}

double ins_angle_diff(double a, double b) {
    return ins_angle_wrap(a - b);
}
