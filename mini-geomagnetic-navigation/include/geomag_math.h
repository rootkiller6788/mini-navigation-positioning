/**
 * geomag_math.h -- Mathematical Utilities for Geomagnetic Navigation
 *
 * L3: Mathematical structures for coordinate transforms, rotations,
 *     spherical geometry, and numerical methods specific to
 *     geomagnetic field computations.
 *
 * Topics:
 *   - 3x3 matrix operations (multiply, transpose, invert)
 *   - Vector cross/dot products
 *   - Quaternion algebra for attitude representation
 *   - Spherical trigonometry (great-circle distance, bearing)
 *   - Rotation matrices (Euler angles, axis-angle, DCM)
 *   - Interpolation (bilinear on magnetic maps)
 *   - Numerical root-finding (for magnetic pole location)
 *
 * References:
 *   Kuipers, "Quaternions and Rotation Sequences" (1999)
 *   Titterton & Weston, "Strapdown Inertial Navigation Technology" (2004)
 *   Press et al., "Numerical Recipes in C" (2007)
 */

#ifndef GEOMAG_MATH_H
#define GEOMAG_MATH_H

#include "geomag_core.h"

/* ========================================================================
 * L3: 3x3 matrix multiply: C = A * B (all row-major)
 * ======================================================================== */
void mat3x3_mult(const double A[9], const double B[9], double C[9]);

/**
 * L3: 3x3 matrix transpose: AT = A^T
 */
void mat3x3_transpose(const double A[9], double AT[9]);

/**
 * L3: 3x3 matrix-vector multiply: y = A * x
 */
void mat3x3_vec_mult(const double A[9], const double x[3], double y[3]);

/**
 * L3: 3x3 matrix determinant.
 *
 * det(A) = a11*(a22*a33-a23*a32) - a12*(a21*a33-a23*a31) + a13*(a21*a32-a22*a31)
 */
double mat3x3_det(const double A[9]);

/**
 * L3: 3x3 matrix inverse (Gauss-Jordan elimination).
 *
 * Computes A^{-1}. Returns 0 on success, -1 if singular (|det| < 1e-15).
 */
int mat3x3_inverse(const double A[9], double Ainv[9]);

/**
 * L3: 3x3 identity matrix.
 */
void mat3x3_identity(double A[9]);

/**
 * L3: Vector dot product: a·b
 */
double vec3_dot(const double a[3], const double b[3]);

/**
 * L3: Vector cross product: c = a x b
 */
void vec3_cross(const double a[3], const double b[3], double c[3]);

/**
 * L3: Vector norm (L2): sqrt(x^2+y^2+z^2)
 */
double vec3_norm(const double v[3]);

/**
 * L3: Normalize vector to unit length. Returns 0 on success, -1 if zero vector.
 */
int vec3_normalize(double v[3]);

/* ========================================================================
 * L3: Quaternion representation (Hamilton convention)
 *
 * q = q_w + q_x*i + q_y*j + q_z*k
 * i*j = k, j*k = i, k*i = j
 * i^2 = j^2 = k^2 = i*j*k = -1
 *
 * Unit quaternion for rotation:
 *   q = [cos(theta/2), sin(theta/2)*axis_x, sin(theta/2)*axis_y, sin(theta/2)*axis_z]
 * ======================================================================== */
typedef struct {
    double w, x, y, z;
} Quaternion;

/**
 * L3: Quaternion multiplication (Hamilton product): r = q * p
 */
void quat_mult(const Quaternion *q, const Quaternion *p, Quaternion *r);

/**
 * L3: Quaternion conjugate: q* = [w, -x, -y, -z]
 */
void quat_conjugate(const Quaternion *q, Quaternion *qc);

/**
 * L3: Quaternion normalization (returns unit quaternion). Returns 0 on success.
 */
int quat_normalize(Quaternion *q);

/**
 * L3: Convert quaternion to 3x3 rotation matrix (row-major).
 *
 * R = [ 1-2(y^2+z^2),    2(xy-wz),      2(xz+wy)    ]
 *     [   2(xy+wz),    1-2(x^2+z^2),    2(yz-wx)    ]
 *     [   2(xz-wy),      2(yz+wx),    1-2(x^2+y^2)  ]
 */
void quat_to_dcm(const Quaternion *q, double R[9]);

/**
 * L3: Convert Euler angles (roll, pitch, yaw) to quaternion.
 *
 * Order: ZYX intrinsic = yaw (psi) about Z, pitch (theta) about Y, roll (phi) about X.
 *
 * q = [c(phi/2)c(theta/2)c(psi/2) + s(phi/2)s(theta/2)s(psi/2),
 *      s(phi/2)c(theta/2)c(psi/2) - c(phi/2)s(theta/2)s(psi/2),
 *      c(phi/2)s(theta/2)c(psi/2) + s(phi/2)c(theta/2)s(psi/2),
 *      c(phi/2)c(theta/2)s(psi/2) - s(phi/2)s(theta/2)c(psi/2)]
 */
void euler_to_quat(double roll, double pitch, double yaw, Quaternion *q);

/**
 * L3: Convert quaternion to Euler angles (roll, pitch, yaw).
 *
 * roll  = atan2(2(wx+yz), 1-2(x^2+y^2))
 * pitch = asin(2(wy-zx))
 * yaw   = atan2(2(wz+xy), 1-2(y^2+z^2))
 */
void quat_to_euler(const Quaternion *q, double *roll, double *pitch, double *yaw);

/**
 * L3: Rotate a vector by a quaternion: v' = q * [0,v] * q*
 *
 * Equivalent to: v' = R(q) * v where R(q) is the DCM.
 */
void quat_rotate_vector(const Quaternion *q, const double v[3], double vp[3]);

/**
 * L3: Spherical linear interpolation (SLERP) between two unit quaternions.
 *
 * slerp(q0, q1, t) = q0*(q0^{-1}*q1)^t
 *                  = [sin((1-t)*Omega)*q0 + sin(t*Omega)*q1] / sin(Omega)
 * where Omega = acos(q0·q1)
 *
 * @param q0, q1  Input unit quaternions
 * @param t        Interpolation parameter [0,1]
 * @param q_out    Output interpolated quaternion
 */
void quat_slerp(const Quaternion *q0, const Quaternion *q1, double t, Quaternion *q_out);

/* ========================================================================
 * L3: Spherical geometry on Earth surface
 * ======================================================================== */

/**
 * L3: Great-circle distance between two points on sphere (Haversine formula).
 *
 * a = sin^2(dlat/2) + cos(lat1)*cos(lat2)*sin^2(dlon/2)
 * c = 2*atan2(sqrt(a), sqrt(1-a))
 * d = R * c
 *
 * Uses spherical Earth approximation (R = GEOMAG_EARTH_RADIUS).
 *
 * @param p1, p2  Geodetic coordinates (altitude ignored)
 * @return Distance [m]
 */
double great_circle_distance(const GeodeticCoord *p1, const GeodeticCoord *p2);

/**
 * L3: Initial bearing (azimuth) from p1 to p2 on great circle.
 *
 * theta = atan2(sin(dlon)*cos(lat2),
 *               cos(lat1)*sin(lat2) - sin(lat1)*cos(lat2)*cos(dlon))
 *
 * @param p1, p2  Geodetic coordinates
 * @return Bearing [deg], 0=North, 90=East
 */
double great_circle_bearing(const GeodeticCoord *p1, const GeodeticCoord *p2);

/**
 * L3: Destination point given start, bearing, and distance.
 *
 * lat2 = asin(sin(lat1)*cos(d/R) + cos(lat1)*sin(d/R)*cos(theta))
 * lon2 = lon1 + atan2(sin(theta)*sin(d/R)*cos(lat1),
 *                     cos(d/R) - sin(lat1)*sin(lat2))
 *
 * @param start    Starting coordinates
 * @param bearing  Bearing [deg]
 * @param distance Distance [m]
 * @param dest     Output destination coordinates
 */
void great_circle_destination(const GeodeticCoord *start, double bearing,
                               double distance, GeodeticCoord *dest);

/* ========================================================================
 * L5: Interpolation for magnetic maps
 * ======================================================================== */

/**
 * L5: Bilinear interpolation on a regular lat-lon grid.
 *
 * Given grid values at (i,j), (i+1,j), (i,j+1), (i+1,j+1),
 * interpolate at fractional position (fx, fy) within the cell.
 *
 * v = (1-fx)*(1-fy)*v(i,j) + fx*(1-fy)*v(i+1,j)
 *   + (1-fx)*fy*v(i,j+1) + fx*fy*v(i+1,j+1)
 *
 * Edge cases: clamps to grid boundaries.
 *
 * @param map    Magnetic map grid
 * @param loc    Query location
 * @param field  Output: interpolated value from total_field array
 * @return 0 on success, -1 if point outside grid
 */
int magmap_bilinear_interpolate(const MagneticMap *map, const GeodeticCoord *loc,
                                 double *field);

/**
 * L5: 2D linear interpolation on a regular grid (generic).
 *
 * @param data        Flattened grid data (row-major)
 * @param nrows, ncols Grid dimensions
 * @param row_start, col_start Starting values for each axis
 * @param row_step, col_step Grid spacing
 * @param row_val, col_val Query position
 * @param result       Output interpolated value
 * @return 0 on success, -1 if out of bounds
 */
int bilinear_interpolate_2d(const double *data, int nrows, int ncols,
                             double row_start, double col_start,
                             double row_step, double col_step,
                             double row_val, double col_val, double *result);

/* ========================================================================
 * L5: Numerical methods
 * ======================================================================== */

/**
 * L5: Golden-section search for 1D minimum of a function.
 *
 * Finds minimum of f(x) in [a,b] to tolerance tol.
 * Uses callback function pointer.
 *
 * Complexity: O(log((b-a)/tol)).
 *
 * @param f     Objective function (user-provided)
 * @param a, b  Search interval bounds
 * @param tol   Tolerance for convergence
 * @param xmin  Output: argmin location
 * @param iter  Output: number of iterations used
 * @return Minimum function value
 */
double golden_section_search(double (*f)(double, void*), void *ctx,
                              double a, double b, double tol,
                              double *xmin, int *iter);

/**
 * L5: 2D gradient descent for magnetic field optimization.
 *
 * Minimizes f(x,y) using numerical gradient approximation.
 *
 * @param f      Objective function f(x,y,ctx)
 * @param ctx    User context pointer
 * @param x0, y0 Initial guess
 * @param step   Step size
 * @param tol    Convergence tolerance
 * @param maxiter Maximum iterations
 * @param x_opt, y_opt Output optimal point
 * @param iter   Output: iterations used
 * @return 0 if converged, -1 if max iterations exceeded
 */
int gradient_descent_2d(double (*f)(double, double, void*), void *ctx,
                         double x0, double y0, double step, double tol,
                         int maxiter, double *x_opt, double *y_opt, int *iter);

/**
 * L5: Magnetic heading compensation using declination.
 *
 * True heading = Magnetic heading + Declination
 * (Declination positive East, following navigation convention)
 *
 * @param mag_heading  Magnetic heading [deg, 0=mag North]
 * @param declination  Local declination [deg, +East]
 * @return True heading [deg, 0=true North]
 */
double mag_to_true_heading(double mag_heading, double declination);

/**
 * L3: Angular difference (shortest signed difference).
 *
 * delta = wrap180(target - source)
 *
 * @param a, b Angles [deg]
 * @return Signed difference [deg], range (-180, 180]
 */
double angle_diff_deg(double a, double b);

/**
 * L3: Wrap angle to [-180, 180] degrees.
 */
double wrap180(double angle_deg);

/**
 * L3: Wrap angle to [0, 360] degrees.
 */
double wrap360(double angle_deg);

#endif /* GEOMAG_MATH_H */
