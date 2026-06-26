/**
 * mini-uwb-localization: UWB Positioning Algorithms
 *
 * Implements multilateration positioning algorithms for converting UWB
 * range measurements into position estimates. Supports linear least squares,
 * nonlinear optimization (Gauss-Newton, Levenberg-Marquardt), and
 * TDoA hyperbolic positioning (Chan, Taylor series).
 *
 * Reference: Zekavat & Buehrer (2019) "Handbook of Position Location"
 * Reference: Chan & Ho (1994) IEEE TSP
 * Reference: Foy (1976) IEEE AES
 *
 * Knowledge Coverage: L4 Fundamental Laws (GDOP, CRLB, FIM)
 *                      L5 Algorithms (LS, WLS, GN, LM, Chan)
 *                      L6 Canonical Problems (Trilateration, TDoA)
 */

#ifndef UWB_POSITIONING_H
#define UWB_POSITIONING_H

#include "uwb_types.h"

typedef enum {
    POS_SOLVER_LINEAR_LS    = 0,
    POS_SOLVER_WEIGHTED_LS  = 1,
    POS_SOLVER_GAUSS_NEWTON = 2,
    POS_SOLVER_LEVENBERG_MARQUARDT = 3,
    POS_SOLVER_CHAN_TDOA    = 4,
    POS_SOLVER_TAYLOR_TDOA  = 5
} uwb_positioning_solver_t;

typedef enum {
    POS_DIM_2D = 2,
    POS_DIM_3D = 3
} uwb_positioning_dim_t;

typedef struct {
    uwb_pos3d_t position;
    uwb_covariance_t covariance;
    double residual_norm;
    double gdop;
    double hdop;
    double vdop;
    double pdop;
    int num_iterations;
    int converged;
    double fisher_trace;
    uwb_positioning_solver_t solver_used;
} uwb_positioning_result_t;

typedef struct {
    uwb_ranging_meas_t *measurements;
    int count;
    uwb_pos3d_t initial_guess;
    double *weights;
} uwb_positioning_input_t;

/*
 * 2D Trilateration (closed-form, 3 anchors)
 * Subtracts eq1 from eq2, eq3 yields linear Ax=b.
 * A = 2*[x2-x1 y2-y1; x3-x1 y3-y1]
 */
int trilateration_2d(const uwb_pos3d_t *anchors, const double *ranges,
                     uwb_pos2d_t *result);

/*
 * 3D Trilateration (closed-form, 4 anchors)
 */
int trilateration_3d(const uwb_pos3d_t *anchors, const double *ranges,
                     uwb_pos3d_t *result);

/*
 * Linear Least Squares Multilateration.
 * For N anchors: min ||A*p - b||_2, p = (A^T A)^-1 A^T b.
 * Complexity: O(N*D^2 + D^3) where D=2 or 3.
 */
int multilateration_linear_ls(const uwb_pos3d_t *anchors, const double *ranges,
                              int num_anchors, uwb_positioning_dim_t dim,
                              uwb_positioning_result_t *result);

/*
 * Weighted Least Squares: p = (A^T W A)^-1 A^T W b, w_i = 1/sigma_i^2.
 */
int multilateration_weighted_ls(const uwb_pos3d_t *anchors, const double *ranges,
                                const double *weights, int num_anchors,
                                uwb_positioning_dim_t dim,
                                uwb_positioning_result_t *result);

/*
 * Gauss-Newton iteration for non-linear least squares:
 * J_k = Jacobian, delta = (J_k^T J_k)^-1 J_k^T * residual
 */
int multilateration_gauss_newton(const uwb_pos3d_t *anchors, const double *ranges,
                                 int num_anchors, uwb_positioning_dim_t dim,
                                 uwb_positioning_result_t *result);

/*
 * Levenberg-Marquardt: adds damping lambda*I to J^T J.
 * lambda adaptively adjusted.
 */
int multilateration_levenberg_marquardt(const uwb_pos3d_t *anchors,
                                        const double *ranges, int num_anchors,
                                        uwb_positioning_dim_t dim,
                                        uwb_positioning_result_t *result);

/*
 * Chan TDoA algorithm: two-step WLS for hyperbolic positioning.
 * Reference: Chan & Ho (1994)
 */
int tdoa_positioning_chan(const uwb_pos3d_t *anchors, const double *tdoa_values,
                          int num_anchors, uwb_positioning_dim_t dim,
                          uwb_positioning_result_t *result);

/*
 * Taylor-series TDoA expansion.
 * Reference: Foy (1976)
 */
int tdoa_positioning_taylor(const uwb_pos3d_t *anchors, const double *tdoa_values,
                            int num_anchors, uwb_positioning_dim_t dim,
                            uwb_positioning_result_t *result);

/*
 * Compute all DOP metrics (GDOP, PDOP, HDOP, VDOP).
 * G = (H^T H)^-1, GDOP = sqrt(trace(G))
 */
void compute_dop_metrics(const uwb_pos3d_t *anchors, int num_anchors,
                         const uwb_pos3d_t *tag_pos,
                         uwb_positioning_result_t *result);

double evaluate_anchor_geometry(const uwb_pos3d_t *anchors, int num_anchors,
                                const uwb_pos3d_t *tag_pos);

/*
 * CRLB for position estimation: sqrt(trace(F^-1))
 * F = (1/sigma^2) * H^T H (Fisher Information Matrix)
 */
double crlb_position(const uwb_pos3d_t *anchors, int num_anchors,
                     const uwb_pos3d_t *tag_pos, double range_variance);

/*
 * Update error metrics incrementally (Welford online algorithm).
 * Reference: Welford (1962) Technometrics 4(3)
 */
void error_metrics_update(uwb_error_metrics_t *metrics,
                          const uwb_pos3d_t *estimated,
                          const uwb_pos3d_t *ground_truth,
                          uwb_positioning_dim_t dim);

/*
 * Compute CEP50, CEP90 from accumulated error samples.
 */
void error_metrics_finalize(uwb_error_metrics_t *metrics);

#endif /* UWB_POSITIONING_H */
