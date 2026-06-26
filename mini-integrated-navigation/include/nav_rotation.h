/**
 * @file nav_rotation.h
 * @brief Rotation Representations for Integrated Navigation
 *
 * L2 Core Concepts: Quaternion, DCM, Euler angle representations.
 * L3 Math Structures: SO(3) rotation group operations.
 *
 * Reference: Kuipers, "Quaternions and Rotation Sequences"
 *            Titterton & Weston, "Strapdown Inertial Navigation Technology"
 */

#ifndef NAV_ROTATION_H
#define NAV_ROTATION_H

#include "nav_common.h"

typedef struct {
    NAV_PRECISION w;
    NAV_PRECISION x;
    NAV_PRECISION y;
    NAV_PRECISION z;
} nav_quat_t;

typedef struct {
    NAV_PRECISION m[9];
} nav_dcm_t;

typedef struct {
    NAV_PRECISION roll;
    NAV_PRECISION pitch;
    NAV_PRECISION yaw;
} nav_euler_t;

typedef struct {
    NAV_PRECISION axis_x;
    NAV_PRECISION axis_y;
    NAV_PRECISION axis_z;
    NAV_PRECISION angle;
} nav_axis_angle_t;

void nav_quat_identity(nav_quat_t *q);
void nav_quat_set(nav_quat_t *q, NAV_PRECISION w, NAV_PRECISION x,
                  NAV_PRECISION y, NAV_PRECISION z);
void nav_quat_copy(nav_quat_t *dst, const nav_quat_t *src);
NAV_PRECISION nav_quat_norm(const nav_quat_t *q);
int nav_quat_normalize(nav_quat_t *q);
void nav_quat_conjugate(nav_quat_t *result, const nav_quat_t *q);
void nav_quat_multiply(nav_quat_t *result, const nav_quat_t *p,
                        const nav_quat_t *q);
void nav_quat_rotate_vector(nav_vector3_t *result, const nav_quat_t *q,
                             const nav_vector3_t *v);
void nav_quat_exp(nav_quat_t *result, const nav_vector3_t *v);
void nav_quat_log(nav_vector3_t *result, const nav_quat_t *q);
void nav_quat_slerp(nav_quat_t *result, const nav_quat_t *q0,
                     const nav_quat_t *q1, NAV_PRECISION t);
void nav_quat_inverse(nav_quat_t *result, const nav_quat_t *q);
void nav_quat_to_dcm(nav_dcm_t *dcm, const nav_quat_t *q);
void nav_dcm_to_quat(nav_quat_t *q, const nav_dcm_t *dcm);
void nav_quat_to_euler(nav_euler_t *euler, const nav_quat_t *q);
void nav_euler_to_quat(nav_quat_t *q, const nav_euler_t *euler);
void nav_dcm_to_euler(nav_euler_t *euler, const nav_dcm_t *dcm);
void nav_dcm_multiply(nav_dcm_t *result, const nav_dcm_t *a, const nav_dcm_t *b);
void nav_dcm_transpose(nav_dcm_t *result, const nav_dcm_t *dcm);
void nav_dcm_rotate_vector(nav_vector3_t *result, const nav_dcm_t *dcm,
                            const nav_vector3_t *v);
void nav_skew_symmetric(NAV_PRECISION m[9], const nav_vector3_t *v);
void nav_quat_kinematics(nav_quat_t *q, const nav_vector3_t *omega,
                          NAV_PRECISION dt);
int nav_triad(nav_dcm_t *dcm,
              const nav_vector3_t *body1, const nav_vector3_t *body2,
              const nav_vector3_t *ref1,  const nav_vector3_t *ref2);
void nav_att_error_to_quat(nav_quat_t *q_correction, const nav_vector3_t *psi);
void nav_quat_apply_correction(nav_quat_t *q, const nav_vector3_t *psi);

#endif /* NAV_ROTATION_H */
