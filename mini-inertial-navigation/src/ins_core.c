/**
 * @file    ins_core.c
 * @brief   Inertial Navigation System core functions
 *
 * Knowledge Coverage:
 *   L1 (Definitions): Coordinate frames, WGS84 parameters, vector/matrix types
 *   L2 (Core Concepts): Earth/gravity models, coordinate transforms
 *   L3 (Math Structures): 3D vector algebra, rotation matrices, cross product
 *   L4 (Fundamental Laws): Somigliana gravity formula, Coriolis effect
 *
 * Reference:
 *   Titterton & Weston (2004), "Strapdown Inertial Navigation Technology", 2nd ed.
 *   Groves (2013), "Principles of GNSS, Inertial, and Multisensor Integrated
 *     Navigation Systems", 2nd ed., Artech House.
 *   Hofmann-Wellenhof & Moritz (2006), "Physical Geodesy", 2nd ed.
 *
 * Course Mapping:
 *   MIT 6.630 - Electromagnetic Waves (coordinate systems)
 *   Stanford EE267 - Virtual Reality (3D orientation, inertial sensing)
 *   Berkeley EE117 - Electromagnetic Fields (vector calculus)
 *   ETH 227-0455 - EM Waves (coordinate transforms)
 *   Illinois ECE 451 - EM Fields (vector analysis)
 */

#include "ins_core.h"
#include <string.h>
#include <stdio.h>
#include <math.h>

/* =========================================================================
 * L2: Gravity Model — Somigliana Formula (WGS84)
 * ========================================================================= */

double ins_gravity_wgs84(double lat, double alt) {
    double sin_lat = sin(lat);
    double sin2_lat = sin_lat * sin_lat;
    double g0;

    /* Somigliana closed formula: g at ellipsoid surface */
    double numerator = 1.0 + INS_K_SOMIGLIANA * sin2_lat;
    double denominator = sqrt(1.0 - INS_WGS84_E2 * sin2_lat);
    g0 = INS_GRAVITY_EQUATOR * numerator / denominator;

    /* Height correction (WGS84 free-air reduction) */
    double m_val = INS_WGS84_WE * INS_WGS84_WE * INS_WGS84_A * INS_WGS84_A
                   * INS_WGS84_B / INS_WGS84_GM;
    double h_over_a = alt / INS_WGS84_A;
    double corr = 1.0 - 2.0 * h_over_a * (1.0 + INS_WGS84_F + m_val
                   - 2.0 * INS_WGS84_F * sin2_lat)
                   + 3.0 * h_over_a * h_over_a;

    return g0 * corr;
}

void ins_gravity_ned(double lat, double alt, ins_vec3_t *g_ned) {
    if (!g_ned) return;
    double g_mag = ins_gravity_wgs84(lat, alt);
    ins_vec3_set(g_ned, 0.0, 0.0, g_mag);
}

/* =========================================================================
 * L2: Coordinate Transform — Geodetic <-> ECEF
 * ========================================================================= */

void ins_geodetic_to_ecef(const ins_geodetic_t *geo, ins_ecef_t *ecef) {
    if (!geo || !ecef) return;

    double sin_lat = sin(geo->lat);
    double cos_lat = cos(geo->lat);
    double sin_lon = sin(geo->lon);
    double cos_lon = cos(geo->lon);

    double N = INS_WGS84_A / sqrt(1.0 - INS_WGS84_E2 * sin_lat * sin_lat);

    ecef->x = (N + geo->alt) * cos_lat * cos_lon;
    ecef->y = (N + geo->alt) * cos_lat * sin_lon;
    ecef->z = (N * (1.0 - INS_WGS84_E2) + geo->alt) * sin_lat;
}

void ins_ecef_to_geodetic(const ins_ecef_t *ecef, ins_geodetic_t *geo) {
    if (!ecef || !geo) return;

    double p = sqrt(ecef->x * ecef->x + ecef->y * ecef->y);
    double z = ecef->z;
    double e2 = INS_WGS84_E2;
    double a = INS_WGS84_A;
    double b = INS_WGS84_B;
    double ep2 = INS_WGS84_EP2;

    double lon = atan2(ecef->y, ecef->x);
    double lat;

    /* Zhu (1994) closed-form refinement */
    double beta = atan2(z, p);
    double u = atan2(b * tan(beta), a);
    double sin_u = sin(u);
    double cos_u = cos(u);
    lat = atan2(z + ep2 * b * sin_u * sin_u * sin_u,
                p - e2 * a * cos_u * cos_u * cos_u);

    double sin_lat = sin(lat);
    double N = a / sqrt(1.0 - e2 * sin_lat * sin_lat);
    double alt = p / cos(lat) - N;

    /* Handle poles */
    if (p < 1e-12) {
        lat = (z > 0) ? M_PI / 2.0 : -M_PI / 2.0;
        alt = fabs(z) - b;
    }

    geo->lat = lat;
    geo->lon = lon;
    geo->alt = alt;
}

/* =========================================================================
 * L2: Rotation Matrix — ECEF to NED
 *
 * C_ned^ecef = R_y(-lat - pi/2) * R_z(lon)
 *
 * R_y(theta) at theta = -(lat+pi/2):
 *   [ -sin_lat  0  -cos_lat ]
 *   [     0     1      0    ]
 *   [  cos_lat  0  -sin_lat ]
 *
 * R_z(lon):
 *   [ cos_lon  -sin_lon  0 ]
 *   [ sin_lon   cos_lon  0 ]
 *   [    0         0     1 ]
 * ========================================================================= */

void ins_ecef_to_ned_dcm(double lat, double lon, ins_mat3_t *C) {
    if (!C) return;

    double slat = sin(lat), clat = cos(lat);
    double slon = sin(lon), clon = cos(lon);

    /* C = Ry(-lat-pi/2) * Rz(lon) */
    /* Row 0 of Ry: [-slat, 0, -clat] */
    C->m[0] = -slat * clon;           /* (-slat)*clon + 0*slon + (-clat)*0 */
    C->m[1] = -slat * (-slon);        /* (-slat)*(-slon) + 0*clon + (-clat)*0 */
    C->m[2] = -clat;                  /* (-slat)*0 + 0*0 + (-clat)*1 */
    /* Row 1 of Ry: [0, 1, 0] */
    C->m[3] = slon;                   /* 0*clon + 1*slon + 0*0 */
    C->m[4] = clon;                   /* 0*(-slon) + 1*clon + 0*0 */
    C->m[5] = 0.0;                    /* 0*0 + 1*0 + 0*1 */
    /* Row 2 of Ry: [clat, 0, -slat] */
    C->m[6] = clat * clon;            /* clat*clon + 0*slon + (-slat)*0 */
    C->m[7] = clat * (-slon);         /* clat*(-slon) + 0*clon + (-slat)*0 */
    C->m[8] = -slat;                  /* clat*0 + 0*0 + (-slat)*1 */
}

/* =========================================================================
 * L2: Transport Rate and Earth Rate in NED Frame
 * ========================================================================= */

void ins_transport_rate(const ins_vec3_t *v_ned, double lat, double alt,
                        ins_vec3_t *omega) {
    if (!v_ned || !omega) return;

    double M = ins_meridian_radius(lat);
    double N = ins_prime_vertical_radius(lat);
    double cos_lat = cos(lat);

    omega->x = v_ned->y / (N + alt);
    omega->y = -v_ned->x / (M + alt);
    if (fabs(cos_lat) < 1e-12) {
        omega->z = 0.0;
    } else {
        omega->z = -v_ned->y * tan(lat) / (N + alt);
    }
}

void ins_earth_rate_ned(double lat, ins_vec3_t *omega) {
    if (!omega) return;
    omega->x = INS_WGS84_WE * cos(lat);
    omega->y = 0.0;
    omega->z = -INS_WGS84_WE * sin(lat);
}

/* =========================================================================
 * L2: Vector Operations (ins_vec3_t)
 * ========================================================================= */

void ins_vec3_zero(ins_vec3_t *v) {
    if (!v) return;
    v->x = v->y = v->z = 0.0;
}

void ins_vec3_set(ins_vec3_t *v, double x, double y, double z) {
    if (!v) return;
    v->x = x; v->y = y; v->z = z;
}

void ins_vec3_copy(const ins_vec3_t *src, ins_vec3_t *dst) {
    if (!src || !dst) return;
    dst->x = src->x; dst->y = src->y; dst->z = src->z;
}

void ins_vec3_add(const ins_vec3_t *a, const ins_vec3_t *b, ins_vec3_t *c) {
    if (!a || !b || !c) return;
    c->x = a->x + b->x;
    c->y = a->y + b->y;
    c->z = a->z + b->z;
}

void ins_vec3_sub(const ins_vec3_t *a, const ins_vec3_t *b, ins_vec3_t *c) {
    if (!a || !b || !c) return;
    c->x = a->x - b->x;
    c->y = a->y - b->y;
    c->z = a->z - b->z;
}

void ins_vec3_scale(const ins_vec3_t *a, double s, ins_vec3_t *b) {
    if (!a || !b) return;
    b->x = a->x * s;
    b->y = a->y * s;
    b->z = a->z * s;
}

double ins_vec3_dot(const ins_vec3_t *a, const ins_vec3_t *b) {
    if (!a || !b) return 0.0;
    return a->x * b->x + a->y * b->y + a->z * b->z;
}

void ins_vec3_cross(const ins_vec3_t *a, const ins_vec3_t *b, ins_vec3_t *c) {
    if (!a || !b || !c) return;
    c->x = a->y * b->z - a->z * b->y;
    c->y = a->z * b->x - a->x * b->z;
    c->z = a->x * b->y - a->y * b->x;
}

double ins_vec3_norm(const ins_vec3_t *v) {
    if (!v) return 0.0;
    return sqrt(v->x * v->x + v->y * v->y + v->z * v->z);
}

int ins_vec3_normalize(ins_vec3_t *v) {
    if (!v) return -1;
    double n = ins_vec3_norm(v);
    if (n < 1e-15) return -1;
    double inv_n = 1.0 / n;
    v->x *= inv_n; v->y *= inv_n; v->z *= inv_n;
    return 0;
}

double ins_vec3_angle(const ins_vec3_t *a, const ins_vec3_t *b) {
    if (!a || !b) return 0.0;
    double dot = ins_vec3_dot(a, b);
    double na = ins_vec3_norm(a);
    double nb = ins_vec3_norm(b);
    if (na < 1e-15 || nb < 1e-15) return 0.0;
    double cos_theta = dot / (na * nb);
    if (cos_theta > 1.0) cos_theta = 1.0;
    if (cos_theta < -1.0) cos_theta = -1.0;
    return acos(cos_theta);
}

/* =========================================================================
 * L2: Matrix Operations (ins_mat3_t, row-major)
 * ========================================================================= */

void ins_vec3_skew(const ins_vec3_t *v, ins_mat3_t *S) {
    if (!v || !S) return;
    S->m[0] = 0.0;     S->m[1] = -v->z;   S->m[2] = v->y;
    S->m[3] = v->z;    S->m[4] = 0.0;     S->m[5] = -v->x;
    S->m[6] = -v->y;   S->m[7] = v->x;    S->m[8] = 0.0;
}

void ins_mat3_mul_vec(const ins_mat3_t *A, const ins_vec3_t *x, ins_vec3_t *y) {
    if (!A || !x || !y) return;
    y->x = A->m[0] * x->x + A->m[1] * x->y + A->m[2] * x->z;
    y->y = A->m[3] * x->x + A->m[4] * x->y + A->m[5] * x->z;
    y->z = A->m[6] * x->x + A->m[7] * x->y + A->m[8] * x->z;
}

void ins_mat3_mul(const ins_mat3_t *A, const ins_mat3_t *B, ins_mat3_t *C) {
    if (!A || !B || !C) return;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double sum = 0.0;
            for (int k = 0; k < 3; k++) {
                sum += A->m[i * 3 + k] * B->m[k * 3 + j];
            }
            C->m[i * 3 + j] = sum;
        }
    }
}

void ins_mat3_transpose(const ins_mat3_t *A, ins_mat3_t *AT) {
    if (!A || !AT) return;
    AT->m[0] = A->m[0]; AT->m[1] = A->m[3]; AT->m[2] = A->m[6];
    AT->m[3] = A->m[1]; AT->m[4] = A->m[4]; AT->m[5] = A->m[7];
    AT->m[6] = A->m[2]; AT->m[7] = A->m[5]; AT->m[8] = A->m[8];
}

void ins_mat3_identity(ins_mat3_t *M) {
    if (!M) return;
    memset(M->m, 0, sizeof(M->m));
    M->m[0] = M->m[4] = M->m[8] = 1.0;
}

/* =========================================================================
 * L2: Meridian and Prime Vertical Radii of Curvature
 * ========================================================================= */

double ins_meridian_radius(double lat) {
    double sin_lat = sin(lat);
    double sin2_lat = sin_lat * sin_lat;
    double denom = 1.0 - INS_WGS84_E2 * sin2_lat;
    return INS_WGS84_A * (1.0 - INS_WGS84_E2) / (denom * sqrt(denom));
}

double ins_prime_vertical_radius(double lat) {
    double sin_lat = sin(lat);
    return INS_WGS84_A / sqrt(1.0 - INS_WGS84_E2 * sin_lat * sin_lat);
}

/* =========================================================================
 * L3: Coriolis Acceleration Computation
 *
 * In the rotating Earth reference frame (NED), the Coriolis
 * pseudo-acceleration acting on a moving body is:
 *
 *   a_coriolis = -(2 * omega_ie^n + omega_en^n) x v^n
 *
 * where omega_ie^n = Earth rotation rate vector in NED
 *       omega_en^n = transport rate (NED frame motion over ellipsoid)
 *       v^n = velocity in NED frame
 *
 * This term is subtracted in the velocity integration to account
 * for the rotating reference frame.
 * ========================================================================= */

void ins_coriolis_accel(const ins_vec3_t *v_ned, double lat, double alt,
                        ins_vec3_t *a_coriolis) {
    if (!v_ned || !a_coriolis) return;

    ins_vec3_t w_ie, w_en, w_total;
    ins_earth_rate_ned(lat, &w_ie);
    ins_transport_rate(v_ned, lat, alt, &w_en);

    w_total.x = 2.0 * w_ie.x + w_en.x;
    w_total.y = 2.0 * w_ie.y + w_en.y;
    w_total.z = 2.0 * w_ie.z + w_en.z;

    ins_vec3_cross(v_ned, &w_total, a_coriolis);
}
