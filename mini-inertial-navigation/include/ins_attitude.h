#ifndef INS_ATTITUDE_H
#define INS_ATTITUDE_H
#include "ins_core.h"
#ifdef __cplusplus
extern "C" {
#endif

/* =========================================================================
 * L1: Quaternion Definition
 *
 * Quaternion representation of attitude (body-to-NED):
 *   q = [q0, q1, q2, q3] = cos(theta/2) + u*sin(theta/2)
 *   where u is the unit rotation axis and theta is the rotation magnitude.
 *
 * Reference: Kuipers (1999), "Quaternions and Rotation Sequences",
 *   Princeton University Press.
 * Reference: Markley & Crassidis (2014), "Fundamentals of Spacecraft
 *   Attitude Determination and Control", Springer.
 * ========================================================================= */

/* -------------------------------------------------------------------------
 * L1: Attitude Representation Types
 * ------------------------------------------------------------------------- */

/** Quaternion [w, x, y, z] where w is scalar part */
typedef struct {
    double w, x, y, z;
} ins_quat_t;

/** Euler angles [rad] following aerospace NED convention */
typedef struct {
    double roll;
    double pitch;
    double yaw;
} ins_euler_t;

/* -------------------------------------------------------------------------
 * L2: Quaternion Core Operations
 * ------------------------------------------------------------------------- */

/** Set quaternion to identity (zero rotation): q = [1, 0, 0, 0] */
void ins_quat_identity(ins_quat_t *q);

/** Copy quaternion */
void ins_quat_copy(const ins_quat_t *src, ins_quat_t *dst);

/** Compute quaternion norm squared: |q|^2 = w^2 + x^2 + y^2 + z^2 */
double ins_quat_norm_sq(const ins_quat_t *q);

/** Compute quaternion norm: |q| = sqrt(w^2 + x^2 + y^2 + z^2) */
double ins_quat_norm(const ins_quat_t *q);

/**
 * Normalize quaternion to unit length.
 * Returns 0 on success, -1 if norm is too small.
 * Essential for maintaining valid rotation representation.
 */
int ins_quat_normalize(ins_quat_t *q);

/**
 * Quaternion multiplication: r = p * q
 * Represents rotation composition: first q, then p.
 *
 * Formula:
 *   r.w = p.w*q.w - p.x*q.x - p.y*q.y - p.z*q.z
 *   r.x = p.w*q.x + p.x*q.w + p.y*q.z - p.z*q.y
 *   r.y = p.w*q.y - p.x*q.z + p.y*q.w + p.z*q.x
 *   r.z = p.w*q.z + p.x*q.y - p.y*q.x + p.z*q.w
 *
 * Complexity: 16 multiplications, 12 additions
 */
void ins_quat_mul(const ins_quat_t *p, const ins_quat_t *q, ins_quat_t *r);

/**
 * Quaternion conjugate: q* = [w, -x, -y, -z]
 * For unit quaternion: q* = q^(-1)
 */
void ins_quat_conjugate(const ins_quat_t *q, ins_quat_t *qc);

/**
 * Quaternion inverse: q^(-1) = q* / |q|^2
 * For unit quaternions this equals the conjugate.
 */
void ins_quat_inverse(const ins_quat_t *q, ins_quat_t *qi);

/* -------------------------------------------------------------------------
 * L2: Quaternion to/from Rotation Representations
 * ------------------------------------------------------------------------- */

/**
 * Convert quaternion to Direction Cosine Matrix (DCM).
 * For unit quaternion q = [w, x, y, z]:
 *
 *     [ 1-2(yy+zz)   2(xy-wz)     2(xz+wy)   ]
 * C = [ 2(xy+wz)     1-2(xx+zz)   2(yz-wx)   ]
 *     [ 2(xz-wy)     2(yz+wx)     1-2(xx+yy) ]
 *
 * This is the body-to-NED transformation: v_ned = C * v_body
 */
void ins_quat_to_dcm(const ins_quat_t *q, ins_mat3_t *C);

/**
 * Convert DCM to quaternion using Shepperd's method (1978)
 * for numerical stability.
 */
void ins_dcm_to_quat(const ins_mat3_t *C, ins_quat_t *q);

/**
 * Convert quaternion to Euler angles (roll, pitch, yaw).
 * NED convention with 3-2-1 (yaw-pitch-roll) rotation sequence.
 *
 *   roll  = atan2(2(q.w*q.x + q.y*q.z), 1 - 2(q.x^2 + q.y^2))
 *   pitch = asin(2(q.w*q.y - q.z*q.x))
 *   yaw   = atan2(2(q.w*q.z + q.x*q.y), 1 - 2(q.y^2 + q.z^2))
 *
 * Pitch bounded to [-pi/2, pi/2]; gimbal lock at pitch = +/- pi/2.
 */
void ins_quat_to_euler(const ins_quat_t *q, ins_euler_t *e);

/**
 * Convert Euler angles to quaternion.
 * Uses 3-2-1 (yaw-pitch-roll) rotation sequence.
 */
void ins_euler_to_quat(const ins_euler_t *e, ins_quat_t *q);

/* -------------------------------------------------------------------------
 * L2: Vector Rotation by Quaternion
 * ------------------------------------------------------------------------- */

/**
 * Rotate a 3D vector by quaternion: v_rotated = q * v * q^(-1)
 * Efficient implementation avoiding full quaternion multiply.
 *
 * @param q      Unit quaternion representing rotation
 * @param v      Input vector
 * @param v_rot  Output rotated vector
 */
void ins_quat_rotate_vector(const ins_quat_t *q, const ins_vec3_t *v,
                             ins_vec3_t *v_rot);

/* -------------------------------------------------------------------------
 * L3: Quaternion Kinematics
 * ------------------------------------------------------------------------- */

/**
 * Quaternion time derivative from angular velocity.
 *
 * dq/dt = 0.5 * q * omega_q
 * where omega_q = [0, wx, wy, wz]
 *
 * Explicitly:
 *   dq0/dt = -0.5 * (wx*q1 + wy*q2 + wz*q3)
 *   dq1/dt =  0.5 * (wx*q0 + wz*q2 - wy*q3)
 *   dq2/dt =  0.5 * (wy*q0 - wz*q1 + wx*q3)
 *   dq3/dt =  0.5 * (wz*q0 + wy*q1 - wx*q2)
 *
 * This is the fundamental equation for strapdown attitude propagation.
 *
 * @param q       Current quaternion
 * @param omega   Angular velocity [rad/s] in body frame
 * @param dq_dt   Output quaternion derivative
 */
void ins_quat_kinematics(const ins_quat_t *q, const ins_vec3_t *omega,
                          ins_quat_t *dq_dt);

/* -------------------------------------------------------------------------
 * L5: Quaternion Integration Algorithms
 * ------------------------------------------------------------------------- */

/**
 * First-order quaternion update (Euler method).
 * q(t+dt) = q(t) + dq/dt * dt, followed by normalization.
 */
void ins_quat_update_euler(ins_quat_t *q, const ins_vec3_t *omega, double dt);

/**
 * Exact quaternion update for constant angular velocity.
 * q(t+dt) = q(t) * [cos(|omega|*dt/2), sin(|omega|*dt/2)*omega/|omega|]
 */
void ins_quat_update_exact(ins_quat_t *q, const ins_vec3_t *omega, double dt);

/**
 * Third-order quaternion update with coning compensation (Bortz, 1971).
 * phi = omega*dt + 1/12 * cross(omega_prev, omega) * dt^2
 * q(t+dt) = q(t) * q(phi)
 *
 * Reference: Bortz (1971), IEEE Trans. Aerosp. Electron. Syst., 7(1): 61-66.
 */
void ins_quat_update_coning(ins_quat_t *q, const ins_vec3_t *omega_prev,
                             const ins_vec3_t *omega_curr, double dt);

/* -------------------------------------------------------------------------
 * L5: Sculling Compensation
 * ------------------------------------------------------------------------- */

/**
 * Sculling velocity compensation.
 *
 * Two-sample algorithm (Savage, 1998):
 *   delta_v_scull = 2/3 * (delta_alpha1 x delta_v2 + delta_v1 x delta_alpha2)
 *
 * Reference: Savage (1998), J. Guid. Control Dyn., 21(1): 19-28.
 */
void ins_sculling_compensation(const ins_vec3_t *delta_alpha1,
                                const ins_vec3_t *delta_alpha2,
                                const ins_vec3_t *delta_v1,
                                const ins_vec3_t *delta_v2,
                                ins_vec3_t *delta_v_scull);

/* -------------------------------------------------------------------------
 * L3: Euler Angle Utilities
 * ------------------------------------------------------------------------- */

/** Convert Euler angles to DCM using 3-2-1 (yaw-pitch-roll) sequence. */
void ins_euler_to_dcm(const ins_euler_t *e, ins_mat3_t *C);

/** Convert DCM to Euler angles. Handles gimbal-lock singularity. */
void ins_dcm_to_euler(const ins_mat3_t *C, ins_euler_t *e);

/** Wrap angle to [-pi, pi] range. */
double ins_angle_wrap(double angle);

/** Angular difference clipped to [-pi, pi]. */
double ins_angle_diff(double a, double b);

#ifdef __cplusplus
}
#endif
#endif /* INS_ATTITUDE_H */
