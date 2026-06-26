#ifndef GNSS_POSITION_H
#define GNSS_POSITION_H
#include "gnss_common.h"
#include "gnss_pseudorange.h"

#ifdef __cplusplus
extern "C" {
#endif

/* =========================================================================
 * L1: Position-Velocity-Time (PVT) solution structure
 * ========================================================================= */

/** Single GNSS PVT solution (epoch result) */
typedef struct {
    gnss_ecef_t  pos;          /* Receiver position in ECEF [m] */
    double       clock_bias;   /* Receiver clock bias [m] (c·dt_r) */
    double       clock_drift;  /* Receiver clock drift [m/s] */
    int          n_sats;       /* Number of satellites used */
    double       gdop;         /* Geometric Dilution of Precision */
    double       pdop;         /* Position DOP */
    double       hdop;         /* Horizontal DOP */
    double       vdop;         /* Vertical DOP */
    double       tdop;         /* Time DOP */
    double       rms_residual; /* RMS post-fit residual [m] */
    int          valid;        /* 1 = valid solution */
    int          iterations;   /* LS iteration count */
} gnss_pvt_solution_t;

/** Weighting strategy for weighted least squares */
typedef enum {
    GNSS_WGT_UNIFORM      = 0,  /* All satellites weight = 1 */
    GNSS_WGT_ELEVATION    = 1,  /* Weight = sin²(elevation) */
    GNSS_WGT_SNR          = 2,  /* Weight ∝ C/N₀ */
    GNSS_WGT_ELEV_SNR     = 3   /* Weight = sin²(elevation) × C/N₀ */
} gnss_weight_strategy_t;

/* -------------------------------------------------------------------------
 * L1: Satellite measurement bundle for PVT
 * ------------------------------------------------------------------------- */

/**
 * @brief Corrected pseudorange for one satellite
 *
 * After applying all corrections (clock, iono, tropo, Sagnac),
 * the corrected pseudorange relates to receiver position by:
 *   P_corr = ||sat_pos - rx_pos|| + b
 * where b = c·dt_r is the receiver clock bias in meters.
 */
typedef struct {
    gnss_ecef_t  sat_pos;      /* Satellite ECEF position [m] */
    double       pseudorange;  /* Corrected pseudorange [m] */
    double       elevation;    /* Elevation angle [rad] */
    double       azimuth;      /* Azimuth angle [rad] */
    double       weight;       /* Weight for WLS (1.0 for uniform) */
    int          used;         /* 1 = used in solution, 0 = excluded */
} gnss_sat_meas_t;

/* -------------------------------------------------------------------------
 * L2: Design Matrix (H-matrix)
 *
 * Row i for satellite i:
 *   h_i = [-e_x, -e_y, -e_z, 1]
 * where e = (sat_i - rx_0) / ||sat_i - rx_0|| is the LOS unit vector.
 *
 * The linearized observation equation:
 *   ΔP = H · Δx
 * where Δx = [Δx_rx, Δy_rx, Δz_rx, Δb]^T
 *       ΔP = [ΔP₁, ΔP₂, ..., ΔPₙ]^T
 * ------------------------------------------------------------------------- */

typedef struct {
    int      n_sats;           /* Number of satellites */
    int      n_params;         /* 4 for position+clock */
    double  *data;             /* Row-major [n_sats × 4] */
} gnss_design_matrix_t;

/* -------------------------------------------------------------------------
 * L3: Normal Equation Matrix (4×4)
 *
 * N = Hᵀ·W·H, b = Hᵀ·W·ΔP, Δx = N⁻¹·b
 * ------------------------------------------------------------------------- */

typedef struct {
    double N[4][4];            /* Normal matrix */
    double rhs[4];             /* Right-hand side */
} gnss_normal_eqn_t;

/* -------------------------------------------------------------------------
 * L3: Bancroft algorithm (direct closed-form solution)
 *
 * Bancroft, S. (1985). "An Algebraic Solution of the GPS Equations."
 * IEEE Trans. AES, 21(1), 56-59.
 *
 * Solves the nonlinear system directly using Lorentz inner product,
 * reducing to a quadratic equation. Does not require an initial
 * position guess. Used for initialization of iterative solvers.
 * ------------------------------------------------------------------------- */

typedef struct {
    gnss_ecef_t pos1;          /* First solution (near Earth) */
    gnss_ecef_t pos2;          /* Second solution (often far from Earth) */
    double      bias1;         /* Clock bias for solution 1 [m] */
    double      bias2;         /* Clock bias for solution 2 [m] */
    int         valid;         /* 1 = real solutions obtained */
} gnss_bancroft_result_t;

/* -------------------------------------------------------------------------
 * L3: Iterative Least Squares State
 * ------------------------------------------------------------------------- */

typedef struct {
    gnss_ecef_t  pos;          /* Current position estimate [m] */
    double       clock_bias;   /* Current clock bias estimate [m] */
    double       residual_rms; /* RMS residual [m] */
    int          iteration;    /* Iteration number */
    double       delta_norm;   /* ||Δx|| this iteration */
    int          converged;    /* 1 = converged */
} gnss_ls_state_t;

/* -------------------------------------------------------------------------
 * API: Position Solution Algorithms
 * ------------------------------------------------------------------------- */

/**
 * @brief Bancroft direct solution (closed-form)
 *
 * Uses Minkowski/Lorentz inner product to solve GPS pseudorange equations
 * algebraically. Returns two possible solutions; the one closer to Earth's
 * surface is the correct user position.
 *
 * Input: n_sats corrected pseudorange measurements.
 * Requires: n_sats ≥ 4.
 *
 * Complexity: O(n_sats) matrix algebra, no iteration.
 */
int gnss_bancroft_solve(const gnss_sat_meas_t *meas, int n_sats,
                         gnss_bancroft_result_t *result);

/**
 * @brief Iterative least squares PVT solution
 *
 * Gauss-Newton iteration on the nonlinear observation equations:
 *   1. Start with initial position guess (or Bancroft solution)
 *   2. Compute range residuals: Δρ_i = P_corr,i - ||sat_i - pos|| - b
 *   3. Build design matrix H (n_sats × 4)
 *   4. Solve Δx = (HᵀH)⁻¹·Hᵀ·Δρ
 *   5. Update: pos += Δx[0:3], b += Δx[3]
 *   6. Repeat until ||Δx|| < tol or max_iter
 *
 * @param meas        Corrected satellite measurements [n_sats]
 * @param n_sats      Number of satellites (≥4)
 * @param initial_pos Initial position guess (e.g., Bancroft or last fix)
 * @param max_iter    Maximum iterations (typ. 10)
 * @param tol         Convergence tolerance [m]
 * @param solution    Output PVT solution
 * @return 0 = success, <0 = failure
 *
 * Reference: IS-GPS-200, Tsui (2005) Ch.5
 */
int gnss_ls_position_solve(const gnss_sat_meas_t *meas, int n_sats,
                            gnss_ecef_t initial_pos,
                            int max_iter, double tol,
                            gnss_pvt_solution_t *solution);

/**
 * @brief Weighted least squares PVT solution
 *
 * Same as LS but with diagonal weight matrix W.
 * Solves Δx = (HᵀWH)⁻¹·HᵀW·Δρ.
 *
 * Weight schemes:
 *   - UNIFORM:    w_i = 1.0
 *   - ELEVATION:  w_i = sin²(el_i)  (de-weights low-elevation sats)
 *   - SNR:        w_i = 10^(C/N₀/10)
 *   - ELEV_SNR:   w_i = sin²(el_i) · 10^(C/N₀/10)
 */
int gnss_wls_position_solve(const gnss_sat_meas_t *meas, int n_sats,
                             gnss_ecef_t initial_pos,
                             gnss_weight_strategy_t wgt_strategy,
                             int max_iter, double tol,
                             gnss_pvt_solution_t *solution);

/**
 * @brief Velocity estimation from Doppler measurements
 *
 * v_rx = (HᵀH)⁻¹·Hᵀ·ṁ_corr
 * where ṁ_corr is the corrected range-rate for each satellite.
 *
 * Range rate equation:
 *   ṁ = -(Δf / f_carrier)·c  = e·(v_sat - v_rx) + clock_drift
 */
int gnss_velocity_solve(const gnss_sat_meas_t *meas, int n_sats,
                         const double doppler_hz[], const double carrier_freq,
                         const double sat_vel[][3],
                         gnss_ecef_t rx_pos, double *clock_drift,
                         double rx_vel[3]);

/* -------------------------------------------------------------------------
 * API: Design matrix and Normal equations
 * ------------------------------------------------------------------------- */

int  gnss_design_matrix_build(const gnss_sat_meas_t *meas, int n_sats,
                               gnss_ecef_t rx_guess,
                               gnss_design_matrix_t *H);
void gnss_design_matrix_free(gnss_design_matrix_t *H);

int gnss_normal_eqn_solve(const gnss_design_matrix_t *H,
                           const double residuals[], int n_sats,
                           const double weights[],
                           double delta_x[4]);

int gnss_compute_dop_from_H(const gnss_design_matrix_t *H,
                             double *gdop, double *pdop,
                             double *hdop, double *vdop, double *tdop);

/** @brief RAIM fault detection via residual-based method */
int gnss_raim_fde(const gnss_sat_meas_t *meas, int n_sats,
                   gnss_ecef_t rx_guess, double sigma_threshold,
                   int excluded_mask[]);

/** @brief Compute weighted RMS of post-fit residuals */
double gnss_postfit_rms(const gnss_sat_meas_t *meas, int n_sats,
                         gnss_ecef_t rx_pos, double clock_bias);

#ifdef __cplusplus
}
#endif
#endif /* GNSS_POSITION_H */
