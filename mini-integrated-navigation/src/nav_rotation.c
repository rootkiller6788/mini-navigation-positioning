/**
 * @file nav_rotation.c
 * @brief Rotation Representations Implementation
 */

#include "nav_rotation.h"
#include <math.h>
#include <string.h>

void nav_quat_identity(nav_quat_t *q) {
    if (!q) return;
    q->w = 1.0; q->x = 0.0; q->y = 0.0; q->z = 0.0;
}

void nav_quat_set(nav_quat_t *q, NAV_PRECISION w, NAV_PRECISION x,
                  NAV_PRECISION y, NAV_PRECISION z) {
    if (!q) return;
    q->w = w; q->x = x; q->y = y; q->z = z;
}

void nav_quat_copy(nav_quat_t *dst, const nav_quat_t *src) {
    if (!dst || !src) return;
    memcpy(dst, src, sizeof(nav_quat_t));
}

NAV_PRECISION nav_quat_norm(const nav_quat_t *q) {
    if (!q) return 0.0;
    return sqrt(q->w*q->w + q->x*q->x + q->y*q->y + q->z*q->z);
}

int nav_quat_normalize(nav_quat_t *q) {
    if (!q) return -1;
    NAV_PRECISION n = nav_quat_norm(q);
    if (n < 1e-15) return -1;
    NAV_PRECISION inv = 1.0 / n;
    q->w *= inv; q->x *= inv; q->y *= inv; q->z *= inv;
    return 0;
}

void nav_quat_conjugate(nav_quat_t *result, const nav_quat_t *q) {
    if (!result || !q) return;
    result->w =  q->w;
    result->x = -q->x;
    result->y = -q->y;
    result->z = -q->z;
}

void nav_quat_multiply(nav_quat_t *result, const nav_quat_t *p,
                        const nav_quat_t *q) {
    if (!result || !p || !q) return;
    NAV_PRECISION pw=p->w, px=p->x, py=p->y, pz=p->z;
    NAV_PRECISION qw=q->w, qx=q->x, qy=q->y, qz=q->z;
    result->w = pw*qw - px*qx - py*qy - pz*qz;
    result->x = pw*qx + px*qw + py*qz - pz*qy;
    result->y = pw*qy - px*qz + py*qw + pz*qx;
    result->z = pw*qz + px*qy - py*qx + pz*qw;
}

void nav_quat_rotate_vector(nav_vector3_t *result, const nav_quat_t *q,
                             const nav_vector3_t *v) {
    if (!result || !q || !v) return;
    NAV_PRECISION qw=q->w, qx=q->x, qy=q->y, qz=q->z;
    NAV_PRECISION vx=v->x, vy=v->y, vz=v->z;
    NAV_PRECISION tw = -qx*vx - qy*vy - qz*vz;
    NAV_PRECISION tx =  qw*vx + qy*vz - qz*vy;
    NAV_PRECISION ty =  qw*vy + qz*vx - qx*vz;
    NAV_PRECISION tz =  qw*vz + qx*vy - qy*vx;
    result->x = tw*(-qx) + tx*qw + ty*(-qz) - tz*(-qy);
    result->y = tw*(-qy) - tx*(-qz) + ty*qw + tz*(-qx);
    result->z = tw*(-qz) + tx*(-qy) - ty*(-qx) + tz*qw;
}

void nav_quat_exp(nav_quat_t *result, const nav_vector3_t *v) {
    if (!result || !v) return;
    NAV_PRECISION norm = sqrt(v->x*v->x + v->y*v->y + v->z*v->z);
    if (norm < 1e-15) { nav_quat_identity(result); return; }
    NAV_PRECISION half = 0.5 * norm;
    NAV_PRECISION s = sin(half) / norm;
    result->w = cos(half);
    result->x = v->x * s; result->y = v->y * s; result->z = v->z * s;
}

void nav_quat_log(nav_vector3_t *result, const nav_quat_t *q) {
    if (!result || !q) return;
    NAV_PRECISION w = q->w;
    if (w > 1.0) w = 1.0;
    if (w < -1.0) w = -1.0;
    NAV_PRECISION theta = acos(w);
    if (fabs(theta) < 1e-15) {
        result->x = 0.0; result->y = 0.0; result->z = 0.0;
        return;
    }
    NAV_PRECISION s = theta / sin(theta);
    result->x = q->x * s; result->y = q->y * s; result->z = q->z * s;
}

void nav_quat_slerp(nav_quat_t *result, const nav_quat_t *q0,
                     const nav_quat_t *q1, NAV_PRECISION t) {
    if (!result || !q0 || !q1) return;
    NAV_PRECISION dot = q0->w*q1->w + q0->x*q1->x + q0->y*q1->y + q0->z*q1->z;
    nav_quat_t q1_mod;
    if (dot < 0.0) {
        q1_mod.w = -q1->w; q1_mod.x = -q1->x;
        q1_mod.y = -q1->y; q1_mod.z = -q1->z;
        dot = -dot;
    } else { memcpy(&q1_mod, q1, sizeof(nav_quat_t)); }
    if (dot > 1.0) dot = 1.0;
    NAV_PRECISION theta = acos(dot);
    if (fabs(theta) < 1e-10) {
        NAV_PRECISION t1 = 1.0 - t;
        result->w = t1*q0->w + t*q1_mod.w;
        result->x = t1*q0->x + t*q1_mod.x;
        result->y = t1*q0->y + t*q1_mod.y;
        result->z = t1*q0->z + t*q1_mod.z;
        nav_quat_normalize(result);
        return;
    }
    NAV_PRECISION st = sin(theta);
    NAV_PRECISION s0 = sin((1.0-t)*theta)/st;
    NAV_PRECISION s1 = sin(t*theta)/st;
    result->w = s0*q0->w + s1*q1_mod.w;
    result->x = s0*q0->x + s1*q1_mod.x;
    result->y = s0*q0->y + s1*q1_mod.y;
    result->z = s0*q0->z + s1*q1_mod.z;
}

void nav_quat_inverse(nav_quat_t *result, const nav_quat_t *q) {
    if (!result || !q) return;
    NAV_PRECISION n2 = q->w*q->w + q->x*q->x + q->y*q->y + q->z*q->z;
    if (n2 < 1e-30) { memset(result, 0, sizeof(nav_quat_t)); return; }
    NAV_PRECISION inv = 1.0/n2;
    result->w = q->w*inv; result->x = -q->x*inv;
    result->y = -q->y*inv; result->z = -q->z*inv;
}

void nav_quat_to_dcm(nav_dcm_t *dcm, const nav_quat_t *q) {
    if (!dcm || !q) return;
    NAV_PRECISION w=q->w, x=q->x, y=q->y, z=q->z;
    NAV_PRECISION xx=x*x, yy=y*y, zz=z*z;
    NAV_PRECISION xy=x*y, xz=x*z, yz=y*z;
    NAV_PRECISION wx=w*x, wy=w*y, wz=w*z;
    dcm->m[0] = 1.0 - 2.0*(yy+zz);
    dcm->m[1] = 2.0*(xy - wz);
    dcm->m[2] = 2.0*(xz + wy);
    dcm->m[3] = 2.0*(xy + wz);
    dcm->m[4] = 1.0 - 2.0*(xx+zz);
    dcm->m[5] = 2.0*(yz - wx);
    dcm->m[6] = 2.0*(xz - wy);
    dcm->m[7] = 2.0*(yz + wx);
    dcm->m[8] = 1.0 - 2.0*(xx+yy);
}

void nav_dcm_to_quat(nav_quat_t *q, const nav_dcm_t *dcm) {
    if (!q || !dcm) return;
    NAV_PRECISION m00=dcm->m[0], m11=dcm->m[4], m22=dcm->m[8];
    NAV_PRECISION trace = m00+m11+m22;
    if (trace > 0.0) {
        NAV_PRECISION s = sqrt(trace+1.0)*2.0;
        q->w = 0.25*s;
        q->x = (dcm->m[7]-dcm->m[5])/s;
        q->y = (dcm->m[2]-dcm->m[6])/s;
        q->z = (dcm->m[3]-dcm->m[1])/s;
    } else if (m00 > m11 && m00 > m22) {
        NAV_PRECISION s = sqrt(1.0+m00-m11-m22)*2.0;
        q->w = (dcm->m[7]-dcm->m[5])/s; q->x = 0.25*s;
        q->y = (dcm->m[1]+dcm->m[3])/s;
        q->z = (dcm->m[2]+dcm->m[6])/s;
    } else if (m11 > m22) {
        NAV_PRECISION s = sqrt(1.0+m11-m00-m22)*2.0;
        q->w = (dcm->m[2]-dcm->m[6])/s;
        q->x = (dcm->m[1]+dcm->m[3])/s; q->y = 0.25*s;
        q->z = (dcm->m[5]+dcm->m[7])/s;
    } else {
        NAV_PRECISION s = sqrt(1.0+m22-m00-m11)*2.0;
        q->w = (dcm->m[3]-dcm->m[1])/s;
        q->x = (dcm->m[2]+dcm->m[6])/s;
        q->y = (dcm->m[5]+dcm->m[7])/s; q->z = 0.25*s;
    }
}

void nav_quat_to_euler(nav_euler_t *euler, const nav_quat_t *q) {
    if (!euler || !q) return;
    NAV_PRECISION w=q->w, x=q->x, y=q->y, z=q->z;
    NAV_PRECISION sp = 2.0*(w*y - z*x);
    if (sp > 1.0) sp = 1.0;
    if (sp < -1.0) sp = -1.0;
    euler->pitch = asin(sp);
    if (fabs(sp) > 0.9999) {
        euler->roll = 0.0;
        euler->yaw = atan2(-2.0*(w*x - y*z), 1.0-2.0*(x*x+y*y));
        return;
    }
    euler->roll = atan2(2.0*(w*x+y*z), 1.0-2.0*(x*x+y*y));
    euler->yaw  = atan2(2.0*(w*z+x*y), 1.0-2.0*(y*y+z*z));
}

void nav_euler_to_quat(nav_quat_t *q, const nav_euler_t *euler) {
    if (!q || !euler) return;
    NAV_PRECISION cr=cos(euler->roll*0.5), sr=sin(euler->roll*0.5);
    NAV_PRECISION cp=cos(euler->pitch*0.5), sp=sin(euler->pitch*0.5);
    NAV_PRECISION cy=cos(euler->yaw*0.5), sy=sin(euler->yaw*0.5);
    q->w = cr*cp*cy + sr*sp*sy;
    q->x = sr*cp*cy - cr*sp*sy;
    q->y = cr*sp*cy + sr*cp*sy;
    q->z = cr*cp*sy - sr*sp*cy;
}

void nav_dcm_to_euler(nav_euler_t *euler, const nav_dcm_t *dcm) {
    if (!euler || !dcm) return;
    euler->pitch = asin(-dcm->m[6]);
    if (fabs(dcm->m[6]) > 0.9999) {
        euler->roll = 0.0;
        euler->yaw = atan2(-dcm->m[1], dcm->m[0]);
        return;
    }
    euler->roll = atan2(dcm->m[7], dcm->m[8]);
    euler->yaw  = atan2(dcm->m[3], dcm->m[0]);
}

void nav_dcm_multiply(nav_dcm_t *result, const nav_dcm_t *a, const nav_dcm_t *b) {
    if (!result || !a || !b) return;
    const NAV_PRECISION *am=a->m, *bm=b->m;
    NAV_PRECISION *rm=result->m;
    for (int r=0; r<3; r++)
        for (int c=0; c<3; c++)
            rm[r*3+c] = am[r*3]*bm[c] + am[r*3+1]*bm[3+c] + am[r*3+2]*bm[6+c];
}

void nav_dcm_transpose(nav_dcm_t *result, const nav_dcm_t *dcm) {
    if (!result || !dcm) return;
    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
            result->m[j*3+i] = dcm->m[i*3+j];
}

void nav_dcm_rotate_vector(nav_vector3_t *result, const nav_dcm_t *dcm,
                            const nav_vector3_t *v) {
    if (!result || !dcm || !v) return;
    result->x = dcm->m[0]*v->x + dcm->m[1]*v->y + dcm->m[2]*v->z;
    result->y = dcm->m[3]*v->x + dcm->m[4]*v->y + dcm->m[5]*v->z;
    result->z = dcm->m[6]*v->x + dcm->m[7]*v->y + dcm->m[8]*v->z;
}

void nav_skew_symmetric(NAV_PRECISION m[9], const nav_vector3_t *v) {
    if (!m || !v) return;
    m[0]=0.0;  m[1]=-v->z; m[2]=v->y;
    m[3]=v->z; m[4]=0.0;   m[5]=-v->x;
    m[6]=-v->y; m[7]=v->x;  m[8]=0.0;
}

void nav_quat_kinematics(nav_quat_t *q, const nav_vector3_t *omega,
                          NAV_PRECISION dt) {
    if (!q || !omega) return;
    NAV_PRECISION qw=q->w, qx=q->x, qy=q->y, qz=q->z;
    NAV_PRECISION ox=omega->x, oy=omega->y, oz=omega->z;
    NAV_PRECISION hdt = 0.5*dt;
    q->w += hdt*(-qx*ox - qy*oy - qz*oz);
    q->x += hdt*( qw*ox + qy*oz - qz*oy);
    q->y += hdt*( qw*oy + qz*ox - qx*oz);
    q->z += hdt*( qw*oz + qx*oy - qy*ox);
    nav_quat_normalize(q);
}

int nav_triad(nav_dcm_t *dcm,
              const nav_vector3_t *b1, const nav_vector3_t *b2,
              const nav_vector3_t *r1, const nav_vector3_t *r2) {
    if (!dcm || !b1 || !b2 || !r1 || !r2) return -1;
    nav_vector3_t tb1, tb2, tb3, tr1, tr2, tr3;
    NAV_PRECISION n;
    /* Body triad */
    n = sqrt(b1->x*b1->x + b1->y*b1->y + b1->z*b1->z);
    if (n < 1e-15) return -1;
    tb1.x=b1->x/n; tb1.y=b1->y/n; tb1.z=b1->z/n;
    tb2.x=b1->y*b2->z - b1->z*b2->y;
    tb2.y=b1->z*b2->x - b1->x*b2->z;
    tb2.z=b1->x*b2->y - b1->y*b2->x;
    n = sqrt(tb2.x*tb2.x + tb2.y*tb2.y + tb2.z*tb2.z);
    if (n < 1e-15) return -1;
    tb2.x/=n; tb2.y/=n; tb2.z/=n;
    tb3.x=tb1.y*tb2.z - tb1.z*tb2.y;
    tb3.y=tb1.z*tb2.x - tb1.x*tb2.z;
    tb3.z=tb1.x*tb2.y - tb1.y*tb2.x;
    /* Reference triad */
    n = sqrt(r1->x*r1->x + r1->y*r1->y + r1->z*r1->z);
    if (n < 1e-15) return -1;
    tr1.x=r1->x/n; tr1.y=r1->y/n; tr1.z=r1->z/n;
    tr2.x=r1->y*r2->z - r1->z*r2->y;
    tr2.y=r1->z*r2->x - r1->x*r2->z;
    tr2.z=r1->x*r2->y - r1->y*r2->x;
    n = sqrt(tr2.x*tr2.x + tr2.y*tr2.y + tr2.z*tr2.z);
    if (n < 1e-15) return -1;
    tr2.x/=n; tr2.y/=n; tr2.z/=n;
    tr3.x=tr1.y*tr2.z - tr1.z*tr2.y;
    tr3.y=tr1.z*tr2.x - tr1.x*tr2.z;
    tr3.z=tr1.x*tr2.y - tr1.y*tr2.x;
    /* DCM = T_body * M_ref^T */
    NAV_PRECISION T[9] = {tb1.x,tb2.x,tb3.x, tb1.y,tb2.y,tb3.y, tb1.z,tb2.z,tb3.z};
    NAV_PRECISION M[9] = {tr1.x,tr2.x,tr3.x, tr1.y,tr2.y,tr3.y, tr1.z,tr2.z,tr3.z};
    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++) {
            dcm->m[i*3+j] = 0.0;
            for (int k=0; k<3; k++)
                dcm->m[i*3+j] += T[i*3+k] * M[j*3+k];
        }
    return 0;
}

void nav_att_error_to_quat(nav_quat_t *q_corr, const nav_vector3_t *psi) {
    if (!q_corr || !psi) return;
    nav_quat_set(q_corr, 1.0, 0.5*psi->x, 0.5*psi->y, 0.5*psi->z);
    nav_quat_normalize(q_corr);
}

void nav_quat_apply_correction(nav_quat_t *q, const nav_vector3_t *psi) {
    if (!q || !psi) return;
    nav_quat_t dq, tmp;
    nav_att_error_to_quat(&dq, psi);
    nav_quat_multiply(&tmp, &dq, q);
    memcpy(q, &tmp, sizeof(nav_quat_t));
}
