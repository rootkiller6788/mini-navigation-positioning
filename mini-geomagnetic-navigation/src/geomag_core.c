/**
 * geomag_core.c -- Core implementations
 *
 * Coordinate transformations (geodetic <-> ECEF), magnetic elements
 * computation, magnetometer calibration.
 *
 * Reference: NIMA TR8350.2 (WGS84), Bowring (1976)
 */

#include "geomag_core.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ========================================================================
 * L3: Geodetic to ECEF conversion
 *
 * Formula (WGS84 ellipsoid):
 *   N(phi) = a / sqrt(1 - e^2 * sin^2(phi))    [prime vertical radius]
 *   X = (N + h) * cos(phi) * cos(lambda)
 *   Y = (N + h) * cos(phi) * sin(lambda)
 *   Z = (N * (1 - e^2) + h) * sin(phi)
 *
 * where phi = geodetic latitude [rad], lambda = longitude [rad].
 *
 * Complexity: O(1), purely algebraic.
 * ======================================================================== */
void geodetic_to_ecef(const GeodeticCoord *geo, ECEFCoord *ecef)
{
    double lat_rad = geo->lat * DEG2RAD;
    double lon_rad = geo->lon * DEG2RAD;
    double sin_lat = sin(lat_rad);
    double cos_lat = cos(lat_rad);
    double sin_lon = sin(lon_rad);
    double cos_lon = cos(lon_rad);

    double N = WGS84_A / sqrt(1.0 - WGS84_E2 * sin_lat * sin_lat);

    ecef->x = (N + geo->alt) * cos_lat * cos_lon;
    ecef->y = (N + geo->alt) * cos_lat * sin_lon;
    ecef->z = (N * (1.0 - WGS84_E2) + geo->alt) * sin_lat;
}

/* ========================================================================
 * L3: ECEF to Geodetic conversion (Bowring 1976 iterative method)
 *
 * Algorithm:
 *   1. Longitude is exact: lambda = atan2(Y, X)
 *   2. Iterate for latitude and height:
 *      p = sqrt(X^2 + Y^2)
 *      phi_0 = atan2(Z, p * (1 - e^2))
 *      Repeat:
 *        N_i = a / sqrt(1 - e^2 * sin^2(phi_i))
 *        h_i = p/cos(phi_i) - N_i
 *        phi_{i+1} = atan2(Z, p * (1 - e^2 * N_i/(N_i + h_i)))
 *      Until |phi_{i+1} - phi_i| < epsilon
 *
 * Converges in 2-3 iterations for h < 1e6 m (below LEO).
 *
 * Reference: Bowring, B.R., "Transformation from Spatial to Geographical
 *   Coordinates", Survey Review, Vol. 23, No. 181, 1976, pp. 323-327.
 * Complexity: O(k) where k = number of iterations (typically 2-5).
 * ======================================================================== */
void ecef_to_geodetic(const ECEFCoord *ecef, GeodeticCoord *geo)
{
    double p = sqrt(ecef->x * ecef->x + ecef->y * ecef->y);
    double lon_rad;

    /* Longitude: exact from atan2 */
    if (p < 1e-12) {
        /* At the poles, longitude is arbitrary; set to 0 */
        lon_rad = 0.0;
    } else {
        lon_rad = atan2(ecef->y, ecef->x);
    }
    geo->lon = lon_rad * RAD2DEG;

    /* Iterative Bowring method for latitude and height */
    double lat_rad;
    if (p < 1e-12) {
        /* At poles */
        lat_rad = (ecef->z > 0.0) ? M_PI / 2.0 : -M_PI / 2.0;
        geo->alt = fabs(ecef->z) - WGS84_B;
    } else {
        /* First approximation */
        lat_rad = atan2(ecef->z, p * (1.0 - WGS84_E2));

        double sin_lat, N, h;
        for (int iter = 0; iter < 10; iter++) {
            sin_lat = sin(lat_rad);
            N = WGS84_A / sqrt(1.0 - WGS84_E2 * sin_lat * sin_lat);
            h = p / cos(lat_rad) - N;

            double lat_new = atan2(ecef->z / p,
                                    1.0 - WGS84_E2 * N / (N + h));

            if (fabs(lat_new - lat_rad) < 1e-12)
                break;
            lat_rad = lat_new;
        }
        sin_lat = sin(lat_rad);
        N = WGS84_A / sqrt(1.0 - WGS84_E2 * sin_lat * sin_lat);
        geo->alt = p / cos(lat_rad) - N;
    }
    geo->lat = lat_rad * RAD2DEG;
}

/* ========================================================================
 * L3: ECEF to NED rotation matrix
 *
 * R_ned_ecef(lat, lon) =
 *   [ -sin(phi)*cos(lambda),  -sin(phi)*sin(lambda),  cos(phi) ]
 *   [ -sin(lambda),            cos(lambda),           0        ]
 *   [ -cos(phi)*cos(lambda),  -cos(phi)*sin(lambda),  -sin(phi) ]
 *
 * where phi = geodetic latitude, lambda = longitude.
 *
 * This matrix rotates a vector expressed in ECEF to NED:
 *   v_ned = R_ned_ecef * v_ecef
 *
 * Properties: R is orthonormal (R^T = R^{-1}), det(R) = 1.
 * ======================================================================== */
void ecef_to_ned_rotation(double lat, double lon, double R[9])
{
    double sin_phi = sin(lat);
    double cos_phi = cos(lat);
    double sin_lambda = sin(lon);
    double cos_lambda = cos(lon);

    R[0] = -sin_phi * cos_lambda;  R[1] = -sin_phi * sin_lambda;  R[2] =  cos_phi;
    R[3] = -sin_lambda;             R[4] =  cos_lambda;             R[5] =  0.0;
    R[6] = -cos_phi * cos_lambda;  R[7] = -cos_phi * sin_lambda;  R[8] = -sin_phi;
}

/* ========================================================================
 * L3: Transform magnetic vector from ECEF to NED frame
 *
 * B_ned = R_ned_ecef * B_ecef
 * ======================================================================== */
void mag_ecef_to_ned(const MagVector *B_ecef, double lat, double lon,
                     MagVector *B_ned)
{
    double R[9];
    ecef_to_ned_rotation(lat, lon, R);

    double Bx = B_ecef->bx, By = B_ecef->by, Bz = B_ecef->bz;
    B_ned->bx = R[0] * Bx + R[1] * By + R[2] * Bz;
    B_ned->by = R[3] * Bx + R[4] * By + R[5] * Bz;
    B_ned->bz = R[6] * Bx + R[7] * By + R[8] * Bz;
}

/* ========================================================================
 * L2: Compute magnetic elements from NED field vector
 *
 * X = B_north [nT]     (North component)
 * Y = B_east  [nT]      (East component)
 * Z = B_down  [nT]      (Down component)
 *
 * Horizontal: H = sqrt(X^2 + Y^2)
 * Total:      F = sqrt(X^2 + Y^2 + Z^2)
 * Declination: D = atan2(Y, X)  [deg, +East of true North]
 * Inclination: I = atan2(Z, H)  [deg, +Downward]
 *
 * Edge case: at magnetic poles (H ~ 0), inclination approaches +/-90 deg
 *            and declination is undefined (set to 0).
 *
 * Reference: NOAA, "Magnetic Field Calculators", geomag.bgs.ac.uk.
 * Complexity: O(1).
 * ======================================================================== */
void compute_magnetic_elements(const MagVector *B, MagneticElements *elem)
{
    double X = B->bx;
    double Y = B->by;
    double Z = B->bz;

    double H = sqrt(X * X + Y * Y);
    double F = sqrt(X * X + Y * Y + Z * Z);

    elem->north_component = X;
    elem->east_component  = Y;
    elem->vertical        = Z;
    elem->horizontal      = H;
    elem->total_intensity = F;

    if (H > 1e-3) {
        elem->declination = atan2(Y, X) * RAD2DEG;
        elem->inclination = atan2(Z, H) * RAD2DEG;
    } else {
        /* Near magnetic pole -- declination undefined */
        elem->declination = 0.0;
        if (Z > 0.0)
            elem->inclination = 90.0;
        else if (Z < 0.0)
            elem->inclination = -90.0;
        else
            elem->inclination = 0.0;
    }
}

/* ========================================================================
 * L2: Magnetic field magnitude
 *
 * |B| = sqrt(Bx^2 + By^2 + Bz^2)
 *
 * Complexity: O(1).
 * ======================================================================== */
double mag_magnitude(const MagVector *B)
{
    return sqrt(B->bx * B->bx + B->by * B->by + B->bz * B->bz);
}

/* ========================================================================
 * L2: Magnetic field horizontal intensity
 *
 * H = sqrt(B_north^2 + B_east^2)
 *
 * Complexity: O(1).
 * ======================================================================== */
double mag_horizontal(const MagVector *B)
{
    return sqrt(B->bx * B->bx + B->by * B->by);
}

/* ========================================================================
 * L5: Magnetometer calibration -- hard-iron and soft-iron correction
 *
 * Measurement model:
 *   B_raw = S * M * (B_true + B_hard) + w
 *
 * Calibration (inverse):
 *   B_cal = M^{-1} * S^{-1} * (B_raw - B_hard)
 *
 * Supplied scale_inv[9] is the combined M^{-1} * S^{-1} matrix.
 * bias[3] is the hard-iron offset B_hard.
 *
 * Algorithm:
 *   1. Subtract hard-iron bias: d = B_raw - bias
 *   2. Apply inverse scale/misalignment matrix: B_cal = scale_inv * d
 *
 * Reference: Gebre-Egziabher et al., "Calibration of Strapdown
 *   Magnetometers in Magnetic Field Domain", ASCE J. Aerospace
 *   Engineering, Vol. 19, No. 2, 2006, pp. 87-102.
 *
 * Complexity: O(1) (constant 3x3 matrix-vector multiply).
 * ======================================================================== */
void magnetometer_calibrate(const double raw[3], const double bias[3],
                            const double scale_inv[9], double cal[3])
{
    /* Subtract hard-iron bias */
    double d[3];
    d[0] = raw[0] - bias[0];
    d[1] = raw[1] - bias[1];
    d[2] = raw[2] - bias[2];

    /* Apply inverse soft-iron + scale correction */
    cal[0] = scale_inv[0] * d[0] + scale_inv[1] * d[1] + scale_inv[2] * d[2];
    cal[1] = scale_inv[3] * d[0] + scale_inv[4] * d[1] + scale_inv[5] * d[2];
    cal[2] = scale_inv[6] * d[0] + scale_inv[7] * d[1] + scale_inv[8] * d[2];
}
