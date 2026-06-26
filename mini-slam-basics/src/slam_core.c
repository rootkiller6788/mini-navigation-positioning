/**
 * @file    slam_core.c
 * @brief   Core SLAM utilities: pose composition, transforms, linear algebra
 *
 * Implements foundational operations shared across all SLAM backends.
 * Each function corresponds to an independent knowledge point.
 *
 * Reference: Thrun, Burgard & Fox (2005) "Probabilistic Robotics".
 */

#include "slam_types.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

/* =========================================================================
 * L1: Angle normalization — essential for SE(2) operations
 * ========================================================================= */

/**
 * @brief Normalize angle to [−π, π)
 *
 * Uses floor-based wrap for numerical stability near ±π.
 * Avoids repeated addition/subtraction that accumulates error.
 *
 * Knowledge: Angular wrapping is fundamental to SE(2) operations.
 * All pose operations must normalize angles to prevent drift.
 */
double slam_normalize_angle(double theta) {
    if (theta >= -M_PI && theta < M_PI) return theta;
    double two_pi = 2.0 * M_PI;
    theta = fmod(theta, two_pi);
    if (theta >= M_PI)  theta -= two_pi;
    if (theta < -M_PI)  theta += two_pi;
    /* Handle edge case of exactly -PI → +PI */
    if (theta <= -M_PI + 1e-12) theta = M_PI;
    return theta;
}

/**
 * @brief Angular difference: normalize(a − b)
 *
 * Correctly handles wrap-around. Example: diff(3.1, −3.1) = 6.2 − 2π ≈ −0.083
 * This is crucial for computing heading errors in SLAM.
 */
double slam_angle_diff(double a, double b) {
    return slam_normalize_angle(a - b);
}

/* =========================================================================
 * L1: SE(2) Pose Composition (⊕ operator)
 * ========================================================================= */

/**
 * @brief Compose two SE(2) poses: c = a ⊕ b
 *
 * This is the group operation on SE(2):
 *   c.x   = a.x + b.x·cos(a.θ) − b.y·sin(a.θ)
 *   c.y   = a.y + b.x·sin(a.θ) + b.y·cos(a.θ)
 *   c.θ   = a.θ + b.θ   (normalized)
 *
 * Usage:
 *   - World_from_robot ⊕ robot_from_sensor = world_from_sensor
 *   - Odometry chaining: pose_{t+1} = pose_t ⊕ odom_t
 *
 * The operation is associative but NOT commutative.
 */
void slam_pose_compose(const slam_pose2d_t *a,
                       const slam_pose2d_t *b,
                       slam_pose2d_t *c) {
    double ca = cos(a->theta), sa = sin(a->theta);
    c->x     = a->x + b->x * ca - b->y * sa;
    c->y     = a->y + b->x * sa + b->y * ca;
    c->theta = slam_normalize_angle(a->theta + b->theta);
}

/**
 * @brief Inverse of SE(2) pose: inv(a) ⊕ a = (0,0,0)
 *
 * The inverse of (x, y, θ) is:
 *   inv.x   = −x·cos(θ) − y·sin(θ)
 *   inv.y   =  x·sin(θ) − y·cos(θ)
 *   inv.θ   = −θ
 *
 * Proof: inv(a) ⊕ a = (−x·cos(θ)−y·sin(θ) + x·cos(−θ)−y·sin(−θ), ...)
 *                     = (0, 0, 0)  [since cos(−θ)=cos(θ), sin(−θ)=−sin(θ)]
 */
void slam_pose_inverse(const slam_pose2d_t *a,
                       slam_pose2d_t *inv) {
    double ca = cos(a->theta), sa = sin(a->theta);
    inv->x     = -a->x * ca - a->y * sa;
    inv->y     =  a->x * sa - a->y * ca;
    inv->theta = slam_normalize_angle(-a->theta);
}

/**
 * @brief Relative pose between two world-frame poses: rel = b ⊖ a
 *
 * Defined as: rel = inv(a) ⊕ b
 * This gives the pose of b expressed in a's coordinate frame.
 *
 * Used for computing odometry between two robot poses.
 */
void slam_pose_relative(const slam_pose2d_t *a,
                        const slam_pose2d_t *b,
                        slam_pose2d_t *rel) {
    double dx = b->x - a->x;
    double dy = b->y - a->y;
    double ca = cos(a->theta), sa = sin(a->theta);
    rel->x     =  dx * ca + dy * sa;
    rel->y     = -dx * sa + dy * ca;
    rel->theta = slam_normalize_angle(b->theta - a->theta);
}

/**
 * @brief Transform a point from robot to world frame
 *
 *   p_world = R(θ)·p_robot + [x, y]^T
 *
 * where R(θ) = [cos(θ), −sin(θ); sin(θ), cos(θ)]
 */
void slam_transform_point(const slam_pose2d_t *pose,
                          double px, double py,
                          double *wx, double *wy) {
    double c = cos(pose->theta), s = sin(pose->theta);
    *wx = pose->x + px * c - py * s;
    *wy = pose->y + px * s + py * c;
}

/**
 * @brief Transform a point from world to robot frame (inverse transform)
 *
 *   p_robot = R(−θ)·(p_world − [x, y]^T)
 */
void slam_transform_point_inv(const slam_pose2d_t *pose,
                               double wx, double wy,
                               double *px, double *py) {
    double dx = wx - pose->x;
    double dy = wy - pose->y;
    double c = cos(pose->theta), s = sin(pose->theta);
    *px =  dx * c + dy * s;
    *py = -dx * s + dy * c;
}

/* =========================================================================
 * L3: Rotation matrix construction
 * ========================================================================= */

/**
 * @brief Build 2D rotation matrix from angle θ
 *
 * R(θ) = [cos(θ), −sin(θ);
 *         sin(θ),  cos(θ)]
 *
 * R(θ) ∈ SO(2), a Lie group with Lie algebra so(2).
 * R(θ₁)·R(θ₂) = R(θ₁ + θ₂)  (the group is abelian for 2D).
 */
void slam_rotation_matrix_2d(double theta, double R[4]) {
    double c = cos(theta), s = sin(theta);
    R[0] = c;  R[1] = -s;
    R[2] = s;  R[3] =  c;
}

/**
 * @brief Derivative of rotation matrix w.r.t. θ
 *
 * dR/dθ = [−sin(θ), −cos(θ);
 *          cos(θ), −sin(θ)]
 *
 * Used in SE(2) Jacobians for graph optimization.
 */
void slam_rotation_derivative_2d(double theta, double dR[4]) {
    double c = cos(theta), s = sin(theta);
    dR[0] = -s;  dR[1] = -c;
    dR[2] =  c;  dR[3] = -s;
}

/**
 * @brief SE(2) transformation as 3×3 homogeneous matrix
 *
 * T = [R(θ), [x; y];
 *      0, 0,    1   ]
 *
 * T ∈ SE(2) ⊂ GL(3, ℝ). This representation makes composition
 * a simple matrix multiplication.
 */
void slam_se2_to_matrix(const slam_pose2d_t *pose, double T[9]) {
    double c = cos(pose->theta), s = sin(pose->theta);
    T[0] = c;  T[1] = -s; T[2] = pose->x;
    T[3] = s;  T[4] =  c; T[5] = pose->y;
    T[6] = 0;  T[7] =  0; T[8] = 1.0;
}

/**
 * @brief SE(2) matrix multiplication: C = A·B
 *
 * Homogeneous matrix multiply (3×3)·(3×3) for SE(2) composition.
 * Used to verify that the Lie group formulation is consistent.
 */
void slam_se2_matmul(const double A[9], const double B[9], double C[9]) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double sum = 0.0;
            for (int k = 0; k < 3; k++) {
                sum += A[i*3 + k] * B[k*3 + j];
            }
            C[i*3 + j] = sum;
        }
    }
}

/* =========================================================================
 * L3: Basic Linear Algebra Operations
 * ========================================================================= */

/**
 * @brief Matrix-vector multiply: y = A·x (row-major storage)
 *
 * y_i = Σ_j A_{i,j} · x_j
 *
 * Complexity: O(m·n). Used throughout SLAM for Jacobian application,
 * covariance propagation, and innovation computation.
 */
void slam_matvec_mul(const double *A, const double *x,
                     int m, int n, double *y) {
    for (int i = 0; i < m; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            sum += A[i * n + j] * x[j];
        }
        y[i] = sum;
    }
}

/**
 * @brief Matrix-matrix multiply: C = A·B (row-major)
 *
 * C_{i,j} = Σ_k A_{i,k} · B_{k,j}
 *
 * All matrices stored row-major. Dimensions: A(m×k), B(k×n), C(m×n).
 */
void slam_matmul(const double *A, const double *B,
                 int m, int k, int n, double *C) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int l = 0; l < k; l++) {
                sum += A[i * k + l] * B[l * n + j];
            }
            C[i * n + j] = sum;
        }
    }
}

/**
 * @brief Matrix transpose multiply: C = A^T·B
 *
 * C_{i,j} = Σ_k A_{k,i} · B_{k,j}
 */
void slam_matmul_AT_B(const double *A, const double *B,
                       int k, int m, int n, double *C) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int l = 0; l < k; l++) {
                sum += A[l * m + i] * B[l * n + j];
            }
            C[i * n + j] = sum;
        }
    }
}

/**
 * @brief Matrix multiply by transpose: C = A·B^T
 *
 * C_{i,j} = Σ_k A_{i,k} · B_{j,k}
 */
void slam_matmul_A_BT(const double *A, const double *B,
                       int m, int k, int n, double *C) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int l = 0; l < k; l++) {
                sum += A[i * k + l] * B[j * k + l];
            }
            C[i * n + j] = sum;
        }
    }
}

/**
 * @brief Compute determinant of 2×2 matrix
 *
 * det([a b; c d]) = a·d − b·c
 */
double slam_det2x2(const double A[4]) {
    return A[0]*A[3] - A[1]*A[2];
}

/**
 * @brief Compute inverse of 2×2 matrix
 *
 * A^{-1} = (1/det)·[ d, −b; −c, a ]
 * Returns 0 if singular (|det| < 1e-12).
 */
int slam_inv2x2(const double A[4], double Ainv[4]) {
    double det = slam_det2x2(A);
    if (fabs(det) < 1e-12) return 0;
    double invdet = 1.0 / det;
    Ainv[0] =  A[3] * invdet;
    Ainv[1] = -A[1] * invdet;
    Ainv[2] = -A[2] * invdet;
    Ainv[3] =  A[0] * invdet;
    return 1;
}

/**
 * @brief Cholesky decomposition L·L^T = A for 3×3 SPD matrix
 *
 * Cholesky-Banachiewicz algorithm:
 *   For i=0..n-1:
 *     L[i,i] = √(A[i,i] − Σ_{k=0}^{i-1} L[i,k]²)
 *     For j=i+1..n-1:
 *       L[j,i] = (A[j,i] − Σ_{k=0}^{i-1} L[j,k]·L[i,k]) / L[i,i]
 *
 * Returns 0 if not positive-definite (numerical check).
 * Used in 3×3 innovation covariance inversion for FastSLAM.
 */
int slam_cholesky_3x3(const double A[9], double L[9]) {
    memset(L, 0, 9 * sizeof(double));

    for (int i = 0; i < 3; i++) {
        double sum = 0.0;
        for (int k = 0; k < i; k++) {
            sum += L[i*3 + k] * L[i*3 + k];
        }
        double diag = A[i*3 + i] - sum;
        if (diag <= 1e-12) return 0;
        L[i*3 + i] = sqrt(diag);

        for (int j = i + 1; j < 3; j++) {
            sum = 0.0;
            for (int k = 0; k < i; k++) {
                sum += L[j*3 + k] * L[i*3 + k];
            }
            L[j*3 + i] = (A[j*3 + i] - sum) / L[i*3 + i];
        }
    }
    return 1;
}

/**
 * @brief Solve L·y = b (forward substitution), L is lower-triangular
 *
 * y_i = (b_i − Σ_{j=0}^{i−1} L_{i,j}·y_j) / L_{i,i}
 */
void slam_forward_sub_3x3(const double L[9], const double b[3], double y[3]) {
    y[0] = b[0] / L[0];
    y[1] = (b[1] - L[3]*y[0]) / L[4];
    y[2] = (b[2] - L[6]*y[0] - L[7]*y[1]) / L[8];
}

/**
 * @brief Solve L^T·x = y (back substitution), L^T is upper-triangular
 *
 * x_i = (y_i − Σ_{j=i+1}^{n-1} L^T_{i,j}·x_j) / L^T_{i,i}
 */
void slam_back_sub_3x3(const double L[9], const double y[3], double x[3]) {
    x[2] = y[2] / L[8];
    x[1] = (y[1] - L[7]*x[2]) / L[4];
    x[0] = (y[0] - L[3]*x[1] - L[6]*x[2]) / L[0];
}

/**
 * @brief Solve A·x = b for 3×3 using Cholesky: A = L·L^T
 *
 * 1. Cholesky: A = L·L^T
 * 2. Forward:  L·y = b   → y
 * 3. Back:     L^T·x = y → x
 *
 * Returns 0 if A is not SPD.
 */
int slam_solve_cholesky_3x3(const double A[9], const double b[3], double x[3]) {
    double L[9];
    if (!slam_cholesky_3x3(A, L)) return 0;
    double y[3];
    slam_forward_sub_3x3(L, b, y);
    slam_back_sub_3x3(L, y, x);
    return 1;
}

/* =========================================================================
 * L3: Covariance Operations
 * ========================================================================= */

/**
 * @brief Initialize identity matrix
 */
void slam_eye(double *M, int n) {
    memset(M, 0, n * n * sizeof(double));
    for (int i = 0; i < n; i++) {
        M[i * n + i] = 1.0;
    }
}

/**
 * @brief Copy matrix (row-major)
 */
void slam_mat_copy(const double *src, double *dst, int n) {
    memcpy(dst, src, n * n * sizeof(double));
}

/**
 * @brief Matrix addition: C = A + B
 */
void slam_mat_add(const double *A, const double *B, double *C, int n) {
    for (int i = 0; i < n * n; i++) C[i] = A[i] + B[i];
}

/**
 * @brief Scalar multiply matrix: A = s·A
 */
void slam_mat_scale(double *A, double s, int n) {
    for (int i = 0; i < n * n; i++) A[i] *= s;
}

/**
 * @brief Compute Σ_new = J·Σ·J^T + Q  (common EKF pattern)
 *
 * Used for covariance prediction and observation innovation.
 * Σ(3×3), J(m×3), Q(m×m) row-major.
 *
 * Complexity: O(m²·3) = O(m²) for small m.
 */
void slam_covariance_propagate(const double *Sigma, const double *J,
                                const double *Q, int m, double *S) {
    /* S = J·Σ */
    slam_matmul(J, Sigma, m, 3, 3, S);
    /* S = (J·Σ)·J^T */
    double temp[9]; /* m × m max for our use is 3x3 */
    slam_matmul_A_BT(S, J, m, 3, m, temp);
    /* S = J·Σ·J^T + Q */
    for (int i = 0; i < m * m; i++) {
        S[i] = temp[i] + Q[i];
    }
}

/* =========================================================================
 * L3: Information Form Operations
 * ========================================================================= */

/**
 * @brief Convert covariance to information matrix: Ω = Σ^{-1}
 *
 * The information form (canonical parameterization) is dual to the
 * covariance (moments) parameterization. Key properties:
 *   - Marginalization in moments = conditioning in information
 *   - Conditioning in moments = marginalization in information
 *
 * For SEIF (Sparse Extended Information Filter), the information
 * matrix is approximately sparse, enabling O(N) updates vs O(N²)
 * for EKF-SLAM.
 *
 * Returns 0 if singular.
 */
int slam_cov_to_info(const double *Sigma, double *Omega, int n) {
    if (n == 3) {
        /* Use direct 3×3 inverse via cofactors */
        double det = Sigma[0]*(Sigma[4]*Sigma[8] - Sigma[5]*Sigma[7])
                   - Sigma[1]*(Sigma[3]*Sigma[8] - Sigma[5]*Sigma[6])
                   + Sigma[2]*(Sigma[3]*Sigma[7] - Sigma[4]*Sigma[6]);
        if (fabs(det) < 1e-12) return 0;
        double invdet = 1.0 / det;

        /* Cofactor matrix (adjugate) */
        Omega[0] = (Sigma[4]*Sigma[8] - Sigma[5]*Sigma[7]) * invdet;
        Omega[1] = (Sigma[2]*Sigma[7] - Sigma[1]*Sigma[8]) * invdet;
        Omega[2] = (Sigma[1]*Sigma[5] - Sigma[2]*Sigma[4]) * invdet;
        Omega[3] = (Sigma[5]*Sigma[6] - Sigma[3]*Sigma[8]) * invdet;
        Omega[4] = (Sigma[0]*Sigma[8] - Sigma[2]*Sigma[6]) * invdet;
        Omega[5] = (Sigma[3]*Sigma[2] - Sigma[0]*Sigma[5]) * invdet;
        Omega[6] = (Sigma[3]*Sigma[7] - Sigma[4]*Sigma[6]) * invdet;
        Omega[7] = (Sigma[6]*Sigma[1] - Sigma[0]*Sigma[7]) * invdet;
        Omega[8] = (Sigma[0]*Sigma[4] - Sigma[1]*Sigma[3]) * invdet;
        return 1;
    }
    /* For larger n, use LU decomposition — not implemented */
    return 0;
}

/**
 * @brief Information vector update: ξ = Ω·μ
 *
 * In the information form, the Gaussian is parameterized by
 * (ξ, Ω) rather than (μ, Σ). This is useful for graph SLAM
 * where information updates are additive.
 */
void slam_info_vector(const double *Omega, const double *mu,
                      double *xi, int n) {
    slam_matvec_mul(Omega, mu, n, n, xi);
}

/* =========================================================================
 * L2: Map Management
 * ========================================================================= */

/**
 * @brief Initialize a 2D feature map
 */
int slam_map2d_init(slam_map2d_t *map, int capacity) {
    assert(map != NULL);
    assert(capacity > 0);
    map->capacity = capacity;
    map->count = 0;
    map->landmarks = (slam_landmark2d_t *)calloc(capacity,
                                                  sizeof(slam_landmark2d_t));
    if (!map->landmarks) return SLAM_ERR_MEMORY;
    map->last_update = 0;
    return SLAM_OK;
}

/**
 * @brief Free map memory
 */
void slam_map2d_free(slam_map2d_t *map) {
    if (map && map->landmarks) {
        free(map->landmarks);
        map->landmarks = NULL;
        map->count = 0;
        map->capacity = 0;
    }
}

/**
 * @brief Add landmark to map, returns its index
 */
int slam_map2d_add_landmark(slam_map2d_t *map,
                             const slam_landmark2d_t *lm,
                             int *idx) {
    if (map->count >= map->capacity) return SLAM_ERR_MEMORY;
    map->landmarks[map->count] = *lm;
    map->landmarks[map->count].id = map->count;
    *idx = map->count;
    map->count++;
    return SLAM_OK;
}

/* =========================================================================
 * L4: Determinant computation for monotonicity verification
 * ========================================================================= */

/**
 * @brief Compute determinant of a 3×3 matrix (for robot covariance)
 *
 * det = a·(e·i − f·h) − b·(d·i − f·g) + c·(d·h − e·g)
 */
double slam_det3x3(const double M[9]) {
    return M[0]*(M[4]*M[8] - M[5]*M[7])
         - M[1]*(M[3]*M[8] - M[5]*M[6])
         + M[2]*(M[3]*M[7] - M[4]*M[6]);
}

/**
 * @brief Determinant of leading k×k submatrix of an n×n matrix
 *
 * Used for SLAM monotonicity verification: the determinant of
 * any map sub-block should be non-increasing.
 */
double slam_subdet(const double *M, int n, int k) {
    assert(k <= n);
    /* For k=1, it's just M[0,0] */
    if (k == 1) return M[0];
    /* For small k, compute explicitly via cofactor expansion */
    /* This implementation handles k up to 4 for landmark sub-blocks */
    if (k == 2) {
        return M[0]*M[n+1] - M[1]*M[n];
    }
    if (k == 3) {
        double M3[9] = {
            M[0],      M[1],      M[2],
            M[n],      M[n+1],    M[n+2],
            M[2*n],    M[2*n+1],  M[2*n+2]
        };
        return slam_det3x3(M3);
    }
    /* k=4: 4×4 sub-determinant */
    if (k == 4) {
        double a00=M[0], a01=M[1], a02=M[2], a03=M[3];
        double a10=M[n], a11=M[n+1], a12=M[n+2], a13=M[n+3];
        double a20=M[2*n], a21=M[2*n+1], a22=M[2*n+2], a23=M[2*n+3];
        double a30=M[3*n], a31=M[3*n+1], a32=M[3*n+2], a33=M[3*n+3];

        return a00 * (a11*(a22*a33 - a23*a32) - a12*(a21*a33 - a23*a31) + a13*(a21*a32 - a22*a31))
             - a01 * (a10*(a22*a33 - a23*a32) - a12*(a20*a33 - a23*a30) + a13*(a20*a32 - a22*a30))
             + a02 * (a10*(a21*a33 - a23*a31) - a11*(a20*a33 - a23*a30) + a13*(a20*a31 - a21*a30))
             - a03 * (a10*(a21*a32 - a22*a31) - a11*(a20*a32 - a22*a30) + a12*(a20*a31 - a21*a30));
    }
    return 0.0;
}

/* =========================================================================
 * L2: SLAM System Initialization
 * ========================================================================= */

/**
 * @brief Initialize SLAM configuration with sensible defaults
 */
void slam_config_default(slam_config_t *cfg) {
    assert(cfg != NULL);
    cfg->backend              = SLAM_BACKEND_EKF;
    cfg->sensor_type          = SLAM_SENSOR_RANGE_BEARING;
    cfg->motion_type          = SLAM_MOTION_VELOCITY;
    cfg->da_method            = SLAM_DA_MAHALANOBIS_GATE;
    cfg->max_landmarks        = 100;
    cfg->max_particles        = 50;
    cfg->sigma_r              = 0.1;
    cfg->sigma_b              = 0.02;
    cfg->sigma_v              = 0.1;
    cfg->sigma_omega          = 0.05;
    cfg->mahalanobis_gate     = 5.991;
    cfg->max_range            = 30.0;
    cfg->max_graph_iterations = 50;
    cfg->convergence_thresh   = 1e-5;
    cfg->lm_lambda_init       = 1e-3;
    cfg->loop_closure_radius  = 5.0;
    cfg->loop_closure_min_steps = 20;
    cfg->use_robust_kernel    = true;
    cfg->robust_kernel_delta  = 1.0;
}

/**
 * @brief Initialize SLAM metrics to zero
 */
void slam_metrics_init(slam_metrics_t *m) {
    memset(m, 0, sizeof(slam_metrics_t));
}

/* =========================================================================
 * Self-test
 * ========================================================================= */

#ifdef SLAM_CORE_SELFTEST
#include <stdio.h>
int main(void) {
    printf("=== slam_core self-test ===\n");

    /* Angle normalization */
    assert(fabs(slam_normalize_angle(3.5)) < M_PI);
    assert(fabs(slam_normalize_angle(-3.5)) < M_PI);
    assert(fabs(slam_normalize_angle(0.0)) < 1e-12);
    assert(fabs(slam_normalize_angle(2.0 * M_PI)) < 1e-12);
    printf("  angle: OK\n");

    /* SE(2) composition and inverse */
    slam_pose2d_t a = {1.0, 2.0, M_PI/4};
    slam_pose2d_t b = {0.5, 0.3, M_PI/6};
    slam_pose2d_t c, inv_a, back;
    slam_pose_compose(&a, &b, &c);
    slam_pose_inverse(&a, &inv_a);
    slam_pose_compose(&inv_a, &a, &back);
    assert(fabs(back.x) < 1e-12 && fabs(back.y) < 1e-12
           && fabs(back.theta) < 1e-12);
    printf("  SE(2) compose/inverse: OK\n");

    /* Point transform round-trip */
    double wx, wy, px2, py2;
    slam_transform_point(&a, 1.0, 0.5, &wx, &wy);
    slam_transform_point_inv(&a, wx, wy, &px2, &py2);
    assert(fabs(px2 - 1.0) < 1e-12 && fabs(py2 - 0.5) < 1e-12);
    printf("  transform round-trip: OK\n");

    /* 2×2 inverse */
    double M[4] = {4, 7, 2, 6}, Minv[4];
    assert(slam_inv2x2(M, Minv));
    /* verify: M·Minv ≈ I */
    double I_check[4];
    slam_matmul(M, Minv, 2, 2, 2, I_check);
    assert(fabs(I_check[0]-1.0) < 1e-10 && fabs(I_check[3]-1.0) < 1e-10);
    printf("  2x2 inverse: OK\n");

    /* Cholesky 3×3 */
    double A[9] = {4, 12, -16, 12, 37, -43, -16, -43, 98};
    double L[9], x[3], b[3] = {1, 2, 3};
    int ok = slam_solve_cholesky_3x3(A, b, x);
    assert(ok);
    /* Verify: A·x ≈ b */
    double Ax[3];
    slam_matvec_mul(A, x, 3, 3, Ax);
    assert(fabs(Ax[0]-b[0]) < 1e-10 && fabs(Ax[1]-b[1]) < 1e-10
           && fabs(Ax[2]-b[2]) < 1e-10);
    printf("  Cholesky 3x3: OK\n");

    /* Determinant */
    double M3[9] = {6, 1, 1, 4, -2, 5, 2, 8, 7};
    /* det = 6·(−14−40) − 1·(28−10) + 1·(32+4) = 6·(−54) − 1·18 + 1·36 = −324 − 18 + 36 = −306 */
    double det = slam_det3x3(M3);
    assert(fabs(det - (-306.0)) < 1e-10);
    printf("  det3x3: OK\n");

    printf("ALL TESTS PASSED\n");
    return 0;
}
#endif
