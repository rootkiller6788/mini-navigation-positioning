/**
 * @file indoor_positioning.c
 * @brief Core indoor positioning functions
 *
 * Implements: coordinate transforms, trilateration, TDOA multilateration,
 * weighted centroid, distance/RSSI conversion, quaternion utilities.
 */

#include <math.h>
#include <string.h>
#include <stdio.h>
#include "../include/indoor_positioning.h"

/* ============================================================================
 * L4 - RSSI Path Loss Model (Friis adapted for indoor)
 * ============================================================================ */

double rssi_to_distance(double rssi, const path_loss_model_t *model) {
    if (!model || model->path_loss_exp <= 0.0) {
        return -1.0;
    }
    /* d = d0 * 10^((RSSI_0 - rssi) / (10 * n))
     * d0 = 1.0 meter (reference distance) */
    double exponent = (model->rssi_at_1m - rssi) / (10.0 * model->path_loss_exp);
    double distance = pow(10.0, exponent);

    /* Guard against extreme values from bad RSSI */
    if (distance < 0.01) distance = 0.01;
    if (distance > 1000.0) distance = 1000.0;

    return distance;
}

double distance_to_rssi(double distance, const path_loss_model_t *model) {
    if (!model || distance <= 0.0) {
        return -200.0;  /* Sentinel for invalid */
    }
    /* RSSI(d) = RSSI_0 - 10 * n * log10(d / d0), d0 = 1.0 */
    if (distance < 0.01) distance = 0.01;
    double rssi = model->rssi_at_1m - 10.0 * model->path_loss_exp * log10(distance);
    return rssi;
}

/* ============================================================================
 * L5 - Linearized Least-Squares Trilateration (2D)
 * ============================================================================ */

int trilateration_2d(const position2d_t *anchor_positions,
                     const double *distances,
                     int n_anchors,
                     position2d_t *result) {
    if (!anchor_positions || !distances || !result) {
        return -1;
    }
    if (n_anchors < 3) {
        return -1;  /* Need at least 3 anchors for 2D */
    }

    /* Use anchor 0 as reference for linearization.
     * Equation i: (x-xi)^2 + (y-yi)^2 = di^2
     * Equation 0: (x-x0)^2 + (y-y0)^2 = d0^2
     * Subtract eq0 from eq_i, then solve for x,y:
     *   2*x*(x0 - xi) + 2*y*(y0 - yi) = di^2 - d0^2 + x0^2 - xi^2 + y0^2 - yi^2
     * This gives us N-1 linear equations A*x = b
     */

    double x0 = anchor_positions[0].x;
    double y0 = anchor_positions[0].y;
    double d0 = distances[0];

    /* Build A^T * A and A^T * b */
    double AtA[2][2] = {{0, 0}, {0, 0}};
    double Atb[2] = {0, 0};

    for (int i = 1; i < n_anchors; i++) {
        double xi = anchor_positions[i].x;
        double yi = anchor_positions[i].y;
        double di = distances[i];

        double a0 = 2.0 * (x0 - xi);  /* coefficient for x */
        double a1 = 2.0 * (y0 - yi);  /* coefficient for y */
        double bi = di*di - d0*d0 + x0*x0 - xi*xi + y0*y0 - yi*yi;

        AtA[0][0] += a0 * a0;
        AtA[0][1] += a0 * a1;
        AtA[1][0] += a0 * a1;
        AtA[1][1] += a1 * a1;
        Atb[0] += a0 * bi;
        Atb[1] += a1 * bi;
    }

    /* Solve the 2x2 system: AtA * pos = Atb */
    double det = AtA[0][0] * AtA[1][1] - AtA[0][1] * AtA[1][0];
    if (fabs(det) < 1e-15) {
        return -1;  /* Singular geometry */
    }

    result->x = (Atb[0] * AtA[1][1] - Atb[1] * AtA[0][1]) / det;
    result->y = (AtA[0][0] * Atb[1] - AtA[0][1] * Atb[0]) / det;

    /* Sanity check: result should be within reasonable range */
    if (fabs(result->x) > 1e6 || fabs(result->y) > 1e6) {
        return -1;
    }

    return 0;
}

/* ============================================================================
 * L5 - Non-linear Least-Squares Trilateration (3D, Gauss-Newton)
 * ============================================================================ */

int trilateration_3d(const position3d_t *anchor_positions,
                     const double *distances,
                     int n_anchors,
                     const position3d_t *initial_guess,
                     position3d_t *result,
                     int max_iterations,
                     double tolerance) {
    if (!anchor_positions || !distances || !result || !initial_guess) {
        return -1;
    }
    if (n_anchors < 4) {
        return -1;
    }
    if (max_iterations <= 0) max_iterations = 20;
    if (tolerance <= 0.0) tolerance = 1e-6;

    double x = initial_guess->x;
    double y = initial_guess->y;
    double z = initial_guess->z;

    for (int iter = 0; iter < max_iterations; iter++) {
        /* Build Jacobian H (Nx3) and residual vector r (N) */
        double HtH[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        double Htr[3] = {0, 0, 0};
        double max_correction = 0.0;

        for (int i = 0; i < n_anchors; i++) {
            double dx = x - anchor_positions[i].x;
            double dy = y - anchor_positions[i].y;
            double dz = z - anchor_positions[i].z;
            double dist_pred = sqrt(dx*dx + dy*dy + dz*dz);

            if (dist_pred < 1e-10) {
                dist_pred = 1e-10;  /* Avoid singularity */
            }

            /* Jacobian row: [dx/d, dy/d, dz/d] */
            double jx = dx / dist_pred;
            double jy = dy / dist_pred;
            double jz = dz / dist_pred;

            /* Residual: measured - predicted */
            double r = distances[i] - dist_pred;

            HtH[0][0] += jx * jx;
            HtH[0][1] += jx * jy;
            HtH[0][2] += jx * jz;
            HtH[1][0] += jy * jx;
            HtH[1][1] += jy * jy;
            HtH[1][2] += jy * jz;
            HtH[2][0] += jz * jx;
            HtH[2][1] += jz * jy;
            HtH[2][2] += jz * jz;

            Htr[0] += jx * r;
            Htr[1] += jy * r;
            Htr[2] += jz * r;

            if (fabs(r) > max_correction) max_correction = fabs(r);
        }

        /* Check convergence */
        if (max_correction < tolerance) {
            result->x = x;
            result->y = y;
            result->z = z;
            return 0;
        }

        /* Solve HtH * dx = Htr for correction dx using Cramer's rule */
        double det = HtH[0][0] * (HtH[1][1] * HtH[2][2] - HtH[1][2] * HtH[2][1])
                   - HtH[0][1] * (HtH[1][0] * HtH[2][2] - HtH[1][2] * HtH[2][0])
                   + HtH[0][2] * (HtH[1][0] * HtH[2][1] - HtH[1][1] * HtH[2][0]);

        if (fabs(det) < 1e-15) {
            break;  /* Singular, stop iterating */
        }

        /* dx = HtH^{-1} * Htr */
        double dx = (Htr[0] * (HtH[1][1]*HtH[2][2] - HtH[1][2]*HtH[2][1])
                   - HtH[0][1] * (Htr[1]*HtH[2][2] - HtH[1][2]*Htr[2])
                   + HtH[0][2] * (Htr[1]*HtH[2][1] - HtH[1][1]*Htr[2])) / det;
        double dy = (HtH[0][0] * (Htr[1]*HtH[2][2] - HtH[1][2]*Htr[2])
                   - Htr[0] * (HtH[1][0]*HtH[2][2] - HtH[1][2]*HtH[2][0])
                   + HtH[0][2] * (HtH[1][0]*Htr[2] - Htr[1]*HtH[2][0])) / det;
        double dz = (HtH[0][0] * (HtH[1][1]*Htr[2] - Htr[1]*HtH[2][1])
                   - HtH[0][1] * (HtH[1][0]*Htr[2] - Htr[1]*HtH[2][0])
                   + Htr[0] * (HtH[1][0]*HtH[2][1] - HtH[1][1]*HtH[2][0])) / det;

        x += dx;
        y += dy;
        z += dz;
    }

    /* Return best result even if not fully converged */
    result->x = x;
    result->y = y;
    result->z = z;
    return (max_iterations > 0) ? 0 : -1;
}

/* ============================================================================
 * L5 - Chan's TDOA Multilateration Algorithm
 * ============================================================================ */

int tdoa_multilateration(const position3d_t *anchor_positions,
                         const double *tdoa_measurements,
                         int n_anchors,
                         double speed_of_light,
                         position3d_t *result) {
    if (!anchor_positions || !tdoa_measurements || !result || n_anchors < 4) {
        return -1;
    }
    if (speed_of_light <= 0.0) {
        speed_of_light = SPEED_OF_LIGHT_MPS;
    }

    /* Convert TDOA to range differences: r_i1 = c * tdoa_i1
     * Reference is anchor[0] */
    double r_diff[32];  /* Max anchors in header */
    int M = n_anchors - 1;  /* Number of TDOA measurements */

    for (int i = 0; i < M && i < 32; i++) {
        r_diff[i] = speed_of_light * tdoa_measurements[i];
    }

    /* Anchor coordinates relative to reference anchor */
    double x1 = anchor_positions[0].x;
    double y1 = anchor_positions[0].y;
    double z1 = anchor_positions[0].z;

    /* Build linear system for Chan's first step.
     * Ga * za = h, where za = [x-x1, y-y1, z-z1, R1]^T
     * R1 = distance from user to reference anchor.
     */
    double Ga[128][4];  /* (M x 4) matrix, M <= 31 typically */
    double h_vec[128];
    double W[128];  /* Weights (diagonal of weight matrix) */

    for (int i = 0; i < M; i++) {
        int idx = i + 1;
        double xi = anchor_positions[idx].x;
        double yi = anchor_positions[idx].y;
        double zi = anchor_positions[idx].z;

        double Ki = xi*xi + yi*yi + zi*zi;
        double K1 = x1*x1 + y1*y1 + z1*z1;

        Ga[i][0] = -2.0 * (xi - x1);
        Ga[i][1] = -2.0 * (yi - y1);
        Ga[i][2] = -2.0 * (zi - z1);
        Ga[i][3] = -2.0 * r_diff[i];

        /* For WLS: compute weight using distance estimate */
        double dist_est = sqrt((xi-x1)*(xi-x1) + (yi-y1)*(yi-y1) + (zi-z1)*(zi-z1));
        if (dist_est < 1e-6) dist_est = 1e-6;
        W[i] = 1.0 / (dist_est * dist_est);

        h_vec[i] = r_diff[i] * r_diff[i] - Ki + K1;
    }

    /* Solve weighted least squares: za = (Ga^T * W * Ga)^{-1} * Ga^T * W * h
     * For simplicity, we solve with equal weights and then iterate */
    double GaTGa[4][4] = {{0}};
    double GaTh[4] = {0};

    for (int i = 0; i < M; i++) {
        for (int row = 0; row < 4; row++) {
            double wi_Ga = W[i] * Ga[i][row];
            GaTh[row] += wi_Ga * h_vec[i];
            for (int col = 0; col < 4; col++) {
                GaTGa[row][col] += wi_Ga * Ga[i][row] * Ga[i][col];
            }
        }
    }

    /* Solve 4x4 system via Gaussian elimination with partial pivot */
    double A[4][5];  /* Augmented matrix */
    for (int r = 0; r < 4; r++) {
        for (int c = 0; c < 4; c++) {
            A[r][c] = GaTGa[r][c];
        }
        A[r][4] = GaTh[r];
    }

    /* Gaussian elimination */
    for (int col = 0; col < 4; col++) {
        /* Find pivot */
        int max_row = col;
        double max_val = fabs(A[col][col]);
        for (int row = col + 1; row < 4; row++) {
            if (fabs(A[row][col]) > max_val) {
                max_val = fabs(A[row][col]);
                max_row = row;
            }
        }
        /* Swap rows */
        if (max_row != col) {
            for (int c = col; c < 5; c++) {
                double tmp = A[col][c];
                A[col][c] = A[max_row][c];
                A[max_row][c] = tmp;
            }
        }
        if (fabs(A[col][col]) < 1e-12) return -1;  /* Singular */
        /* Eliminate below */
        for (int row = col + 1; row < 4; row++) {
            double factor = A[row][col] / A[col][col];
            for (int c = col; c < 5; c++) {
                A[row][c] -= factor * A[col][c];
            }
        }
    }
    /* Back substitution */
    double za[4];
    for (int r = 3; r >= 0; r--) {
        double sum = A[r][4];
        for (int c = r + 1; c < 4; c++) {
            sum -= A[r][c] * za[c];
        }
        za[r] = sum / A[r][r];
    }

    /* First-step solution */
    double x_rel = za[0];  /* x - x1 */
    double y_rel = za[1];  /* y - y1 */
    double z_rel = za[2];  /* z - z1 */

    result->x = x1 + x_rel;
    result->y = y1 + y_rel;
    result->z = z1 + z_rel;

    /* Sanity check */
    if (isnan(result->x) || isnan(result->y) || isnan(result->z) ||
        isinf(result->x) || isinf(result->y) || isinf(result->z)) {
        return -1;
    }

    return 0;
}

/* ============================================================================
 * L5 - Weighted Centroid Localization
 * ============================================================================ */

int weighted_centroid_2d(const position2d_t *anchor_positions,
                         const double *distances,
                         int n_anchors,
                         position2d_t *result) {
    if (!anchor_positions || !distances || !result || n_anchors < 1) {
        return -1;
    }

    double weight_sum = 0.0;
    double wx = 0.0;
    double wy = 0.0;

    for (int i = 0; i < n_anchors; i++) {
        /* Weight = 1/distance to give more weight to closer anchors */
        double w = (distances[i] > 1e-6) ? 1.0 / distances[i] : 1e6;
        /* Cap extreme weights from very small distances */
        if (w > 1e6) w = 1e6;
        wx += w * anchor_positions[i].x;
        wy += w * anchor_positions[i].y;
        weight_sum += w;
    }

    if (weight_sum < 1e-12) {
        /* All distances are zero; return centroid */
        for (int i = 0; i < n_anchors; i++) {
            wx += anchor_positions[i].x;
            wy += anchor_positions[i].y;
        }
        result->x = wx / n_anchors;
        result->y = wy / n_anchors;
    } else {
        result->x = wx / weight_sum;
        result->y = wy / weight_sum;
    }

    return 0;
}

/* ============================================================================
 * L3 - Geodetic to ENU Coordinate Transformation
 * ============================================================================ */

/* WGS-84 ellipsoid parameters */
#define WGS84_A 6378137.0           /* Semi-major axis (meters) */
#define WGS84_F (1.0 / 298.257223563) /* Flattening */
#define WGS84_E2 (2.0 * WGS84_F - WGS84_F * WGS84_F) /* First eccentricity squared */

void geodetic_to_enu(double lat, double lon, double alt,
                     double ref_lat, double ref_lon, double ref_alt,
                     position3d_t *enu) {
    if (!enu) return;

    /* Compute ECEF coordinates from geodetic */
    double sin_lat = sin(lat);
    double cos_lat = cos(lat);
    double sin_lon = sin(lon);
    double cos_lon = cos(lon);

    double N = WGS84_A / sqrt(1.0 - WGS84_E2 * sin_lat * sin_lat);
    double x_ecef = (N + alt) * cos_lat * cos_lon;
    double y_ecef = (N + alt) * cos_lat * sin_lon;
    double z_ecef = (N * (1.0 - WGS84_E2) + alt) * sin_lat;

    /* Compute reference ECEF coordinates */
    double sin_ref_lat = sin(ref_lat);
    double cos_ref_lat = cos(ref_lat);
    double sin_ref_lon = sin(ref_lon);
    double cos_ref_lon = cos(ref_lon);

    double N_ref = WGS84_A / sqrt(1.0 - WGS84_E2 * sin_ref_lat * sin_ref_lat);
    double x_ref = (N_ref + ref_alt) * cos_ref_lat * cos_ref_lon;
    double y_ref = (N_ref + ref_alt) * cos_ref_lat * sin_ref_lon;
    double z_ref = (N_ref * (1.0 - WGS84_E2) + ref_alt) * sin_ref_lat;

    /* Vector from reference to point in ECEF */
    double dx = x_ecef - x_ref;
    double dy = y_ecef - y_ref;
    double dz = z_ecef - z_ref;

    /* Rotate ECEF difference to ENU:
     * E = -sin(lon) * dx + cos(lon) * dy
     * N = -sin(lat)*cos(lon) * dx - sin(lat)*sin(lon) * dy + cos(lat) * dz
     * U = cos(lat)*cos(lon) * dx + cos(lat)*sin(lon) * dy + sin(lat) * dz
     */
    enu->x = -sin_ref_lon * dx + cos_ref_lon * dy;
    enu->y = -sin_ref_lat * cos_ref_lon * dx - sin_ref_lat * sin_ref_lon * dy
             + cos_ref_lat * dz;
    enu->z = cos_ref_lat * cos_ref_lon * dx + cos_ref_lat * sin_ref_lon * dy
             + sin_ref_lat * dz;
}

/* ============================================================================
 * Quaternion Math Utilities
 * ============================================================================ */

void quaternion_to_euler(const quaternion_t *q,
                         double *roll, double *pitch, double *yaw) {
    if (!q || !roll || !pitch || !yaw) return;

    double w = q->w, x = q->x, y = q->y, z = q->z;

    /* Roll (x-axis rotation) */
    double sinr_cosp = 2.0 * (w * x + y * z);
    double cosr_cosp = 1.0 - 2.0 * (x * x + y * y);
    *roll = atan2(sinr_cosp, cosr_cosp);

    /* Pitch (y-axis rotation) */
    double sinp = 2.0 * (w * y - z * x);
    if (fabs(sinp) >= 1.0) {
        *pitch = copysign(M_PI / 2.0, sinp);
    } else {
        *pitch = asin(sinp);
    }

    /* Yaw (z-axis rotation) */
    double siny_cosp = 2.0 * (w * z + x * y);
    double cosy_cosp = 1.0 - 2.0 * (y * y + z * z);
    *yaw = atan2(siny_cosp, cosy_cosp);
}

void euler_to_quaternion(double roll, double pitch, double yaw, quaternion_t *q) {
    if (!q) return;

    double cr = cos(roll * 0.5);
    double sr = sin(roll * 0.5);
    double cp = cos(pitch * 0.5);
    double sp = sin(pitch * 0.5);
    double cy = cos(yaw * 0.5);
    double sy = sin(yaw * 0.5);

    q->w = cr * cp * cy + sr * sp * sy;
    q->x = sr * cp * cy - cr * sp * sy;
    q->y = cr * sp * cy + sr * cp * sy;
    q->z = cr * cp * sy - sr * sp * cy;

    quaternion_normalize(q);
}

void quaternion_rotate_vector(const double v[3], const quaternion_t *q,
                              double result[3]) {
    if (!v || !q || !result) return;

    /* v' = q * v * q* where v is treated as pure quaternion (0, vx, vy, vz) */
    double w = q->w, x = q->x, y = q->y, z = q->z;

    /* q * v */
    double qv_w = -x * v[0] - y * v[1] - z * v[2];
    double qv_x =  w * v[0] + y * v[2] - z * v[1];
    double qv_y =  w * v[1] + z * v[0] - x * v[2];
    double qv_z =  w * v[2] + x * v[1] - y * v[0];

    /* (q*v) * q* */
    result[0] = qv_w * (-x) + qv_x * w + qv_y * (-z) - qv_z * (-y);
    result[1] = qv_w * (-y) + qv_y * w + qv_z * (-x) - qv_x * (-z);
    result[2] = qv_w * (-z) + qv_z * w + qv_x * (-y) - qv_y * (-x);
}
