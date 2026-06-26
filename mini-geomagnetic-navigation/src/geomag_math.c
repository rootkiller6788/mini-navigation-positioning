/**
 * geomag_math.c -- Mathematical Utilities Implementation
 *
 * L3: 3x3 matrix operations, vector algebra, quaternion rotation,
 *     spherical geometry (great-circle), interpolation, numerical methods.
 *
 * Reference:
 *   Kuipers, "Quaternions and Rotation Sequences" (1999)
 *   Titterton & Weston, "Strapdown Inertial Navigation Technology" (2004)
 *   Press et al., "Numerical Recipes in C" (2007)
 */

#include "geomag_math.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ========================================================================
 * L3: 3x3 matrix multiply: C = A * B (all row-major)
 * C[i][j] = sum_k A[i][k] * B[k][j]
 * ======================================================================== */
void mat3x3_mult(const double A[9], const double B[9], double C[9]) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            C[i * 3 + j] = A[i * 3 + 0] * B[0 * 3 + j]
                         + A[i * 3 + 1] * B[1 * 3 + j]
                         + A[i * 3 + 2] * B[2 * 3 + j];
        }
    }
}

/* L3: Matrix transpose AT = A^T */
void mat3x3_transpose(const double A[9], double AT[9]) {
    AT[0] = A[0]; AT[1] = A[3]; AT[2] = A[6];
    AT[3] = A[1]; AT[4] = A[4]; AT[5] = A[7];
    AT[6] = A[2]; AT[7] = A[5]; AT[8] = A[8];
}

/* L3: Matrix-vector multiply: y = A * x */
void mat3x3_vec_mult(const double A[9], const double x[3], double y[3]) {
    y[0] = A[0] * x[0] + A[1] * x[1] + A[2] * x[2];
    y[1] = A[3] * x[0] + A[4] * x[1] + A[5] * x[2];
    y[2] = A[6] * x[0] + A[7] * x[1] + A[8] * x[2];
}

/* L3: Determinant */
double mat3x3_det(const double A[9]) {
    return A[0] * (A[4] * A[8] - A[5] * A[7])
         - A[1] * (A[3] * A[8] - A[5] * A[6])
         + A[2] * (A[3] * A[7] - A[4] * A[6]);
}

/* L3: Matrix inverse (Gauss-Jordan with partial pivoting) */
int mat3x3_inverse(const double A[9], double Ainv[9]) {
    /* Augmented matrix [A | I] */
    double aug[3][6];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            aug[i][j] = A[i * 3 + j];
            aug[i][j + 3] = (i == j) ? 1.0 : 0.0;
        }
    }

    for (int col = 0; col < 3; col++) {
        /* Find pivot */
        int pivot = col;
        double maxv = fabs(aug[col][col]);
        for (int row = col + 1; row < 3; row++) {
            if (fabs(aug[row][col]) > maxv) {
                maxv = fabs(aug[row][col]);
                pivot = row;
            }
        }
        if (maxv < 1e-15) return -1; /* singular */

        /* Swap rows */
        if (pivot != col) {
            for (int j = 0; j < 6; j++) {
                double tmp = aug[col][j];
                aug[col][j] = aug[pivot][j];
                aug[pivot][j] = tmp;
            }
        }

        /* Normalize pivot row */
        double piv_val = aug[col][col];
        for (int j = 0; j < 6; j++) aug[col][j] /= piv_val;

        /* Eliminate other rows */
        for (int row = 0; row < 3; row++) {
            if (row != col) {
                double factor = aug[row][col];
                for (int j = 0; j < 6; j++)
                    aug[row][j] -= factor * aug[col][j];
            }
        }
    }

    /* Extract inverse */
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            Ainv[i * 3 + j] = aug[i][j + 3];
    return 0;
}

/* L3: Identity matrix */
void mat3x3_identity(double A[9]) {
    A[0] = 1.0; A[1] = 0.0; A[2] = 0.0;
    A[3] = 0.0; A[4] = 1.0; A[5] = 0.0;
    A[6] = 0.0; A[7] = 0.0; A[8] = 1.0;
}

/* ========================================================================
 * L3: Vector operations
 * ======================================================================== */

double vec3_dot(const double a[3], const double b[3]) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

void vec3_cross(const double a[3], const double b[3], double c[3]) {
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

double vec3_norm(const double v[3]) {
    return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

int vec3_normalize(double v[3]) {
    double n = vec3_norm(v);
    if (n < 1e-15) return -1;
    v[0] /= n; v[1] /= n; v[2] /= n;
    return 0;
}

/* ========================================================================
 * L3: Quaternion algebra (Hamilton convention)
 *
 * q = w + xi + yj + zk
 * Multiplication: pq = (p_w*q_w - p_v·q_v, p_w*q_v + q_w*p_v + p_v x q_v)
 * ======================================================================== */

void quat_mult(const Quaternion *q, const Quaternion *p, Quaternion *r) {
    r->w = q->w * p->w - q->x * p->x - q->y * p->y - q->z * p->z;
    r->x = q->w * p->x + q->x * p->w + q->y * p->z - q->z * p->y;
    r->y = q->w * p->y - q->x * p->z + q->y * p->w + q->z * p->x;
    r->z = q->w * p->z + q->x * p->y - q->y * p->x + q->z * p->w;
}

void quat_conjugate(const Quaternion *q, Quaternion *qc) {
    qc->w =  q->w;
    qc->x = -q->x;
    qc->y = -q->y;
    qc->z = -q->z;
}

int quat_normalize(Quaternion *q) {
    double n = sqrt(q->w * q->w + q->x * q->x + q->y * q->y + q->z * q->z);
    if (n < 1e-15) return -1;
    q->w /= n; q->x /= n; q->y /= n; q->z /= n;
    return 0;
}

/* L3: Quaternion to Direction Cosine Matrix */
void quat_to_dcm(const Quaternion *q, double R[9]) {
    double w = q->w, x = q->x, y = q->y, z = q->z;
    double x2 = x * x, y2 = y * y, z2 = z * z;
    double wx = w * x, wy = w * y, wz = w * z;
    double xy = x * y, xz = x * z, yz = y * z;

    R[0] = 1.0 - 2.0 * (y2 + z2);
    R[1] = 2.0 * (xy - wz);
    R[2] = 2.0 * (xz + wy);

    R[3] = 2.0 * (xy + wz);
    R[4] = 1.0 - 2.0 * (x2 + z2);
    R[5] = 2.0 * (yz - wx);

    R[6] = 2.0 * (xz - wy);
    R[7] = 2.0 * (yz + wx);
    R[8] = 1.0 - 2.0 * (x2 + y2);
}

/* L3: Euler angles (ZYX intrinsic) to quaternion */
void euler_to_quat(double roll, double pitch, double yaw, Quaternion *q) {
    double cr = cos(roll * 0.5), sr = sin(roll * 0.5);
    double cp = cos(pitch * 0.5), sp = sin(pitch * 0.5);
    double cy = cos(yaw * 0.5), sy = sin(yaw * 0.5);

    q->w = cr * cp * cy + sr * sp * sy;
    q->x = sr * cp * cy - cr * sp * sy;
    q->y = cr * sp * cy + sr * cp * sy;
    q->z = cr * cp * sy - sr * sp * cy;
}

/* L3: Quaternion to Euler angles */
void quat_to_euler(const Quaternion *q, double *roll, double *pitch, double *yaw) {
    double w = q->w, x = q->x, y = q->y, z = q->z;

    /* Roll (phi) */
    double sinr_cosp = 2.0 * (w * x + y * z);
    double cosr_cosp = 1.0 - 2.0 * (x * x + y * y);
    *roll = atan2(sinr_cosp, cosr_cosp);

    /* Pitch (theta) */
    double sinp = 2.0 * (w * y - z * x);
    if (fabs(sinp) >= 1.0)
        *pitch = copysign(M_PI / 2.0, sinp);
    else
        *pitch = asin(sinp);

    /* Yaw (psi) */
    double siny_cosp = 2.0 * (w * z + x * y);
    double cosy_cosp = 1.0 - 2.0 * (y * y + z * z);
    *yaw = atan2(siny_cosp, cosy_cosp);
}

/* L3: Rotate a vector by quaternion: v' = q * [0,v] * q* */
void quat_rotate_vector(const Quaternion *q, const double v[3], double vp[3]) {
    Quaternion qv = { 0.0, v[0], v[1], v[2] };
    Quaternion qc;
    quat_conjugate(q, &qc);
    Quaternion tmp, result;
    quat_mult(&qv, &qc, &tmp);
    quat_mult(q, &tmp, &result);
    vp[0] = result.x;
    vp[1] = result.y;
    vp[2] = result.z;
}

/* L3: Spherical linear interpolation (SLERP) */
void quat_slerp(const Quaternion *q0, const Quaternion *q1, double t,
                Quaternion *q_out) {
    double dot = q0->w * q1->w + q0->x * q1->x + q0->y * q1->y + q0->z * q1->z;
    double sign = (dot < 0.0) ? -1.0 : 1.0;
    double d = fabs(dot);

    double theta = acos(d);
    double sin_theta = sin(theta);

    if (sin_theta < 1e-10) {
        /* Near-parallel: linear interpolation */
        double t0 = 1.0 - t, t1 = t * sign;
        q_out->w = t0 * q0->w + t1 * q1->w;
        q_out->x = t0 * q0->x + t1 * q1->x;
        q_out->y = t0 * q0->y + t1 * q1->y;
        q_out->z = t0 * q0->z + t1 * q1->z;
    } else {
        double s0 = sin((1.0 - t) * theta) / sin_theta;
        double s1 = sin(t * theta) / sin_theta * sign;
        q_out->w = s0 * q0->w + s1 * q1->w;
        q_out->x = s0 * q0->x + s1 * q1->x;
        q_out->y = s0 * q0->y + s1 * q1->y;
        q_out->z = s0 * q0->z + s1 * q1->z;
    }
}

/* ========================================================================
 * L3: Spherical geometry (great-circle calculations)
 * ======================================================================== */

/* L3: Haversine great-circle distance */
double great_circle_distance(const GeodeticCoord *p1, const GeodeticCoord *p2) {
    double lat1 = p1->lat * DEG2RAD, lon1 = p1->lon * DEG2RAD;
    double lat2 = p2->lat * DEG2RAD, lon2 = p2->lon * DEG2RAD;

    double dlat = lat2 - lat1;
    double dlon = lon2 - lon1;

    double a = sin(dlat * 0.5) * sin(dlat * 0.5)
             + cos(lat1) * cos(lat2) * sin(dlon * 0.5) * sin(dlon * 0.5);
    double c = 2.0 * atan2(sqrt(a), sqrt(1.0 - a));

    return GEOMAG_EARTH_RADIUS * c;
}

/* L3: Great-circle initial bearing */
double great_circle_bearing(const GeodeticCoord *p1, const GeodeticCoord *p2) {
    double lat1 = p1->lat * DEG2RAD, lon1 = p1->lon * DEG2RAD;
    double lat2 = p2->lat * DEG2RAD, lon2 = p2->lon * DEG2RAD;
    double dlon = lon2 - lon1;

    double y = sin(dlon) * cos(lat2);
    double x = cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(dlon);
    return atan2(y, x) * RAD2DEG;
}

/* L3: Destination from great-circle course */
void great_circle_destination(const GeodeticCoord *start, double bearing,
                               double distance, GeodeticCoord *dest) {
    double lat1 = start->lat * DEG2RAD, lon1 = start->lon * DEG2RAD;
    double brg = bearing * DEG2RAD;
    double dR = distance / GEOMAG_EARTH_RADIUS;

    double lat2 = asin(sin(lat1) * cos(dR) + cos(lat1) * sin(dR) * cos(brg));
    double lon2 = lon1 + atan2(sin(brg) * sin(dR) * cos(lat1),
                                cos(dR) - sin(lat1) * sin(lat2));

    dest->lat = lat2 * RAD2DEG;
    dest->lon = lon2 * RAD2DEG;
    dest->alt = start->alt;
}

/* ========================================================================
 * L5: Bilinear interpolation on magnetic map
 * ======================================================================== */

int bilinear_interpolate_2d(const double *data, int nrows, int ncols,
                             double row_start, double col_start,
                             double row_step, double col_step,
                             double row_val, double col_val, double *result) {
    if (!data || nrows < 2 || ncols < 2 || !result) return -1;

    double row_frac = (row_val - row_start) / row_step;
    double col_frac = (col_val - col_start) / col_step;

    int i0 = (int)floor(row_frac);
    int j0 = (int)floor(col_frac);

    if (i0 < 0) i0 = 0;
    if (j0 < 0) j0 = 0;
    if (i0 >= nrows - 1) i0 = nrows - 2;
    if (j0 >= ncols - 1) j0 = ncols - 2;

    double fx = row_frac - i0;
    double fy = col_frac - j0;

    double v00 = data[i0 * ncols + j0];
    double v10 = data[(i0 + 1) * ncols + j0];
    double v01 = data[i0 * ncols + j0 + 1];
    double v11 = data[(i0 + 1) * ncols + j0 + 1];

    *result = (1.0 - fx) * (1.0 - fy) * v00
            + fx * (1.0 - fy) * v10
            + (1.0 - fx) * fy * v01
            + fx * fy * v11;

    return 0;
}

int magmap_bilinear_interpolate(const MagneticMap *map, const GeodeticCoord *loc,
                                 double *field) {
    if (!map || !loc || !field || !map->total_field) return -1;

    return bilinear_interpolate_2d(map->total_field, map->nlat, map->nlon,
                                    map->origin.lat, map->origin.lon,
                                    map->grid_spacing_lat, map->grid_spacing_lon,
                                    loc->lat, loc->lon, field);
}

/* ========================================================================
 * L5: Golden-section search for 1D minimum
 *
 * phi = (1+sqrt(5))/2 = 1.618... (golden ratio)
 * Maintains bracketing triple: a < x1 < x2 < b
 * f(x1) and f(x2) evaluated; discard worst outer segment.
 *
 * Complexity: O(log((b-a)/tol)).
 * ======================================================================== */

double golden_section_search(double (*f)(double, void*), void *ctx,
                              double a, double b, double tol,
                              double *xmin, int *iter) {
    const double phi = 1.618033988749895;
    const double inv_phi = 1.0 / phi;

    double x1 = b - (b - a) * inv_phi;
    double x2 = a + (b - a) * inv_phi;
    double f1 = f(x1, ctx);
    double f2 = f(x2, ctx);
    int niter = 0;

    while (fabs(b - a) > tol) {
        if (f1 < f2) {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = b - (b - a) * inv_phi;
            f1 = f(x1, ctx);
        } else {
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + (b - a) * inv_phi;
            f2 = f(x2, ctx);
        }
        niter++;
        if (niter > 100) break;
    }

    if (f1 < f2) {
        *xmin = x1;
    } else {
        *xmin = x2;
    }
    *iter = niter;
    return f(*xmin, ctx);
}

/* ========================================================================
 * L5: 2D gradient descent with numerical gradient
 * ======================================================================== */

int gradient_descent_2d(double (*f)(double, double, void*), void *ctx,
                         double x0, double y0, double step, double tol,
                         int maxiter, double *x_opt, double *y_opt, int *iter) {
    double x = x0, y = y0;
    double f_curr = f(x, y, ctx);

    for (int k = 0; k < maxiter; k++) {
        double eps = 1e-6;
        double dfdx = (f(x + eps, y, ctx) - f(x - eps, y, ctx)) / (2.0 * eps);
        double dfdy = (f(x, y + eps, ctx) - f(x, y - eps, ctx)) / (2.0 * eps);

        double grad_norm = sqrt(dfdx * dfdx + dfdy * dfdy);
        if (grad_norm < tol) {
            *iter = k;
            *x_opt = x; *y_opt = y;
            return 0;
        }

        x -= step * dfdx / grad_norm;
        y -= step * dfdy / grad_norm;

        double f_new = f(x, y, ctx);
        if (f_new > f_curr) {
            step *= 0.5;
        } else {
            f_curr = f_new;
        }

        *iter = k + 1;
    }

    *x_opt = x; *y_opt = y;
    return -1;
}

/* ========================================================================
 * L5: Magnetic heading computation
 * ======================================================================== */

double mag_to_true_heading(double mag_heading, double declination) {
    double true_heading = mag_heading + declination;
    return wrap360(true_heading);
}

double angle_diff_deg(double a, double b) {
    double diff = a - b;
    return wrap180(diff);
}

double wrap180(double angle_deg) {
    double a = fmod(angle_deg, 360.0);
    if (a > 180.0) a -= 360.0;
    if (a <= -180.0) a += 360.0;
    return a;
}

double wrap360(double angle_deg) {
    double a = fmod(angle_deg, 360.0);
    if (a < 0.0) a += 360.0;
    return a;
}
