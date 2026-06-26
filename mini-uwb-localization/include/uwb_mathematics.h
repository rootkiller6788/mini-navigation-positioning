/**
 * mini-uwb-localization: Linear Algebra and Numerical Utilities
 *
 * Self-contained vector/matrix operations for UWB positioning
 * without external BLAS/LAPACK dependency. Provides least squares
 * solvers, small-matrix inversion, Cholesky decomposition, and
 * singular value decomposition for D < 10.
 *
 * Reference: Golub & Van Loan (2013) "Matrix Computations", 4th Ed
 * Reference: Press et al. (2007) "Numerical Recipes", 3rd Ed
 *
 * Knowledge Coverage: L3 Mathematical Structures (Linear Algebra)
 *                      L5 Algorithms (LS, Cholesky, SVD)
 */

#ifndef UWB_MATHEMATICS_H
#define UWB_MATHEMATICS_H

#include <math.h>
#include <string.h>

#define UWB_MATH_MAX_DIM  12
#define UWB_MATH_EPSILON  1e-12

/*
 * Vector operations
 */

/* v = a * v (scale in-place) */
void vec_scale(double *v, int n, double a);

/* v = v1 + v2 */
void vec_add(double *v, const double *v1, const double *v2, int n);

/* v = v1 - v2 */
void vec_sub(double *v, const double *v1, const double *v2, int n);

/* dot = v1^T * v2 */
double vec_dot(const double *v1, const double *v2, int n);

/* ||v||_2 */
double vec_norm2(const double *v, int n);

/* v = v / ||v||_2 */
void vec_normalize(double *v, int n);

/* result = A * x (matrix-vector multiply, A is m x n row-major) */
void mat_vec_mul(double *result, const double *A, const double *x, int m, int n);

/*
 * Matrix operations (row-major storage)
 */

/* C = A * B, A is m x k, B is k x n, C is m x n */
void mat_mul(double *C, const double *A, const double *B, int m, int k, int n);

/* B = A^T, A is m x n */
void mat_transpose(double *B, const double *A, int m, int n);

/* A = A + alpha * I (add scaled identity to square matrix) */
void mat_add_identity(double *A, int n, double alpha);

/* Set matrix A to identity of size n */
void mat_set_identity(double *A, int n);

/* Copy matrix B = A */
void mat_copy(double *B, const double *A, int m, int n);

/* Frobenius norm = sqrt(sum(A_ij^2)) */
double mat_norm_frobenius(const double *A, int m, int n);

/*
 * Linear system solvers
 */

/*
 * Solve Ax = b for n x n system using Gaussian elimination with
 * partial pivoting. A and b are modified in-place.
 * Complexity: O(n^3). Stable for small n <= UWB_MATH_MAX_DIM.
 * Returns 0 if matrix is singular.
 */
int solve_linear_gauss(double *A, double *b, int n);

/*
 * Solve Ax = b using Cholesky decomposition (A must be symmetric positive definite).
 * A = L * L^T, then forward-substitute L*y = b, back-substitute L^T*x = y.
 * Complexity: O(n^3/3). Returns 0 if A is not SPD.
 */
int solve_linear_cholesky(double *A, double *b, int n);

/*
 * Compute A^-1 using Gaussian elimination with full pivoting.
 * A_inv is n x n output. A is modified in-place.
 * Returns 0 if singular.
 */
int mat_inverse(double *A_inv, double *A, int n);

/*
 * Solve least squares: min ||A*x - b||_2.
 * A is m x n, m >= n (over-determined).
 * Solution: x = (A^T * A)^(-1) * A^T * b via normal equations.
 * Returns residual norm.
 */
double solve_least_squares(const double *A, const double *b, int m, int n,
                           double *x);

/*
 * Weighted least squares: min ||W^(1/2) * (A*x - b)||_2.
 * W is the diagonal weight matrix, weights[i] = 1/sigma_i^2.
 */
double solve_weighted_least_squares(const double *A, const double *b,
                                    const double *weights, int m, int n,
                                    double *x);

/*
 * Singular Value Decomposition for m x n matrix (m >= n).
 * A = U * S * V^T
 * U is m x n, S is n (diagonal), V is n x n.
 * Uses Golub-Reinsch algorithm with householder bidiagonalization.
 */
int svd_decompose(double *A, int m, int n, double *U, double *S, double *V);

/*
 * Compute condition number of matrix A.
 * cond(A) = sigma_max / sigma_min (via SVD).
 */
double mat_condition_number(const double *A, int m, int n);

/*
 * Compute determinant of n x n matrix.
 * Uses LU decomposition with partial pivoting.
 */
double mat_determinant(double *A, int n);

/*
 * Compute trace of n x n matrix.
 */
double mat_trace(const double *A, int n);

/*
 * Check if matrix is symmetric positive definite.
 * Uses Cholesky attempt: succeeds if and only if SPD.
 */
int mat_is_spd(const double *A, int n);

/*
 * Solve upper triangular system U*x = b (back-substitution).
 * U is n x n upper triangular. x overwrites b.
 */
void solve_upper_triangular(const double *U, double *b, int n);

/*
 * Solve lower triangular system L*x = b (forward-substitution).
 * L is n x n lower triangular. x overwrites b.
 */
void solve_lower_triangular(const double *L, double *b, int n);

/*
 * QR decomposition via modified Gram-Schmidt.
 * A is m x n, Q is m x n (orthonormal columns), R is n x n (upper triangular).
 */
void qr_decompose_mgs(const double *A, int m, int n, double *Q, double *R);

/*
 * Compute eigenvalues of symmetric 2x2 or 3x3 matrix using
 * analytical formulas (no iteration needed).
 * Returns number of real eigenvalues.
 */
int symm_eigenvalues_2x2(const double *A, double *eigenvalues);
int symm_eigenvalues_3x3(const double *A, double *eigenvalues);

/*
 * Householder transformation: compute vector v such that
 * (I - 2*v*v^T/(v^T*v)) * x = [||x||, 0, ..., 0]^T
 */
void householder_vector(const double *x, int n, double *v, double *beta);

/*
 * Wilkinson shift for symmetric tridiagonal matrix.
 * Returns eigenvalue estimate closest to T[n-1][n-1].
 */
double wilkinson_shift(const double *diag, const double *offdiag, int n);

/* === Numerical stability utilities === */

/* Stable sqrt(a^2 + b^2) avoiding overflow */
double hypot2(double a, double b);

/* Stable 2x2 matrix inversion */
int mat_inverse_2x2(double *A_inv, const double *A);

/* Stable 3x3 matrix inversion using cofactor expansion */
int mat_inverse_3x3(double *A_inv, const double *A);

#endif /* UWB_MATHEMATICS_H */
