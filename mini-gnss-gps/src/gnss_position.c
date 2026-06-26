/* =========================================================================
 * gnss_position.c — Position solution algorithms: Bancroft, LS, WLS, RAIM
 *
 * Covers L3 (least squares estimation, Gauss-Newton), L4 (nonlinear
 * observation equations → linearized design matrix), L5 (Bancroft
 * closed-form, weighted LS), L6 (single-point positioning problem).
 *
 * References:
 * - Bancroft, S. (1985). "An Algebraic Solution of the GPS Equations."
 *   IEEE Trans. AES, 21(1), 56-59.
 * - Strang, G. & Borre, K. (1997). Linear Algebra, Geodesy, and GPS.
 * - IS-GPS-200 (H-matrix definition)
 * - Parkinson, B.W. & Spilker, J.J. (1996). GPS: Theory and Applications.
 * - Brown, R.G. (1992). "A baseline GPS RAIM scheme." Navigation, 39(2).
 * ========================================================================= */

#include "gnss_position.h"
#include "gnss_common.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

/* -------------------------------------------------------------------------
 * L5: Bancroft direct algebraic solution
 *
 * Minkowski (Lorentz) space approach:
 *
 * Let s_i = (x_i, y_i, z_i) be satellite s_i in ECEF.
 * Let u = (x, y, z) be user position, b = c·dt receiver clock bias.
 *
 * Observation: (s_i - u)² = (ρ_i - b)²  [true range = pseudorange - bias]
 *              = s_i² - 2s_i·u + u² = ρ_i² - 2ρ_i·b + b²
 *
 * Rearranging to Lorentz inner product ⟨a,b⟩ = aᵀ·M·b with M = diag(1,1,1,-1):
 *
 *   Let ⟨a,b⟩ = a₁b₁ + a₂b₂ + a₃b₃ - a₄b₄
 *   Define vectors a_i = (s_i, ρ_i) in R⁴
 *
 * Then: ⟨a_i, a_i⟩ = 2⟨a_i, (u, b)⟩ - ⟨(u, b), (u, b)⟩
 *
 * Let Λ = ⟨(u,b), (u,b)⟩ = u² - b²  (unknown scalar)
 *
 * Then: 2⟨a_i, v⟩ - Λ = ⟨a_i, a_i⟩  where v = (u, b)
 *
 * Writing A·v = 1·Λ + w:
 *   A = [a₁ᵀ; a₂ᵀ; ...; aₙᵀ]     (n×4 matrix)
 *   1 = [½, ½, ..., ½]ᵀ         (n×1 vector)
 *   w_i = ½·⟨a_i, a_i⟩          (n×1 vector)
 *
 * Solve v = A⁺·(1·Λ + w) where A⁺ = (AᵀA)⁻¹Aᵀ (pseudoinverse)
 *
 * Substitute back into Λ equation → quadratic in Λ
 * Quadratic: α·Λ² + β·Λ + γ = 0
 *
 * The smaller |Λ| root corresponds to the near-Earth solution.
 * ------------------------------------------------------------------------- */

int gnss_bancroft_solve(const gnss_sat_meas_t *meas, int n_sats,
                         gnss_bancroft_result_t *result) {
    if (n_sats < 4 || !meas || !result) return -1;

    int i, j;
    double A[32][4];   /* n_sats × 4 (max 32 typical) */
    double w[32];      /* RHS vector */

    /* Build A matrix rows and w vector */
    for (i = 0; i < n_sats; i++) {
        double x = meas[i].sat_pos.x;
        double y = meas[i].sat_pos.y;
        double z = meas[i].sat_pos.z;
        double rho = meas[i].pseudorange;

        A[i][0] = x;
        A[i][1] = y;
        A[i][2] = z;
        A[i][3] = rho; /* Note: Lorentz inner product with Minkowski metric */

        w[i] = 0.5 * (x*x + y*y + z*z - rho*rho);
    }

    /* Compute AᵀA (4×4) and Aᵀw (4×1) */
    double ATA[4][4] = {{0}};
    double ATw[4] = {0};
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            double sum = 0.0;
            int k;
            for (k = 0; k < n_sats; k++) sum += A[k][i] * A[k][j];
            ATA[i][j] = sum;
        }
        double s = 0.0;
        int k;
        for (k = 0; k < n_sats; k++) s += A[k][i] * w[k];
        ATw[i] = s;
    }

    /* Also compute Aᵀ1 */
    double AT1[4] = {0};
    for (i = 0; i < 4; i++) {
        double s = 0.0;
        int k;
        for (k = 0; k < n_sats; k++) s += A[k][i];
        AT1[i] = s / 2.0; /* Each row contribution = 1/2 */
    }

    /* Solve ATA⁻¹·AT1 and ATA⁻¹·ATw */
    /* Use Gauss-Jordan on augmented 4×4 */
    double invATA[4][4];
    double work[4][8]; /* aug: [ATA | I] */

    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            work[i][j] = ATA[i][j];
            work[i][j+4] = (i==j) ? 1.0 : 0.0;
        }
    }

    /* Gauss-Jordan elimination */
    int col;
    for (col = 0; col < 4; col++) {
        /* Pivot */
        int pivot = col;
        double maxv = fabs(work[col][col]);
        int r;
        for (r = col+1; r < 4; r++) {
            if (fabs(work[r][col]) > maxv) {
                maxv = fabs(work[r][col]);
                pivot = r;
            }
        }
        if (maxv < 1e-15) return -1; /* singular */

        if (pivot != col) {
            for (j = 0; j < 8; j++) {
                double t = work[col][j];
                work[col][j] = work[pivot][j];
                work[pivot][j] = t;
            }
        }

        double piv_val = work[col][col];
        for (j = 0; j < 8; j++) work[col][j] /= piv_val;

        for (r = 0; r < 4; r++) {
            if (r == col) continue;
            double factor = work[r][col];
            for (j = 0; j < 8; j++) work[r][j] -= factor * work[col][j];
        }
    }

    /* Extract inverse */
    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            invATA[i][j] = work[i][j+4];

    /* v1 = ATA⁻¹·AT1, v2 = ATA⁻¹·ATw */
    double v1[4], v2[4];
    for (i = 0; i < 4; i++) {
        v1[i] = 0.0; v2[i] = 0.0;
        for (j = 0; j < 4; j++) {
            v1[i] += invATA[i][j] * AT1[j];
            v2[i] += invATA[i][j] * ATw[j];
        }
    }

    /* Form quadratic α·Λ² + β·Λ + γ = 0
     * where α = ⟨v1, v1⟩_Lor = v1₁²+v1₂²+v1₃²-v1₄² (Minkowski norm)
     *       β = 2⟨v1, v2⟩_Lor - 1
     *       γ = ⟨v2, v2⟩_Lor
     */

    double v1_Lor2 = v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2] - v1[3]*v1[3];
    double v12_Lor = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2] - v1[3]*v2[3];
    double v2_Lor2 = v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2] - v2[3]*v2[3];

    double a_quad = v1_Lor2;
    double b_quad = 2.0 * v12_Lor - 1.0;
    double c_quad = v2_Lor2;

    /* Solve quadratic */
    double disc = b_quad*b_quad - 4.0*a_quad*c_quad;
    if (fabs(a_quad) < 1e-15 || disc < 0.0) {
        result->valid = 0;
        return -1;
    }

    double sqrt_disc = sqrt(disc);
    double Lambda1 = (-b_quad + sqrt_disc) / (2.0 * a_quad);
    double Lambda2 = (-b_quad - sqrt_disc) / (2.0 * a_quad);

    /* Both solutions: v = Λ·v1 + v2 */
    double sol1[4], sol2[4];
    for (i = 0; i < 4; i++) {
        sol1[i] = Lambda1 * v1[i] + v2[i];
        sol2[i] = Lambda2 * v1[i] + v2[i];
    }

    /* Select the solution closer to Earth's surface (smaller ||u||² - u_z²/a²) */
    double dist1 = sqrt(sol1[0]*sol1[0] + sol1[1]*sol1[1] + sol1[2]*sol1[2]);
    double dist2 = sqrt(sol2[0]*sol2[0] + sol2[1]*sol2[1] + sol2[2]*sol2[2]);

    double earth_surface = GNSS_WGS84_A;
    /* The solution closer to Earth's surface is the correct one */
    int near1 = (fabs(dist1 - earth_surface) <= fabs(dist2 - earth_surface));

    if (near1) {
        result->pos1.x = sol1[0]; result->pos1.y = sol1[1]; result->pos1.z = sol1[2];
        result->bias1  = sol1[3];
        result->pos2.x = sol2[0]; result->pos2.y = sol2[1]; result->pos2.z = sol2[2];
        result->bias2  = sol2[3];
    } else {
        result->pos1.x = sol2[0]; result->pos1.y = sol2[1]; result->pos1.z = sol2[2];
        result->bias1  = sol2[3];
        result->pos2.x = sol1[0]; result->pos2.y = sol1[1]; result->pos2.z = sol1[2];
        result->bias2  = sol1[3];
    }

    result->valid = 1;
    return 0;
}

/* -------------------------------------------------------------------------
 * L3: Design matrix construction
 *
 * For each satellite i at position (x_i, y_i, z_i) with receiver at
 * approximate position (x₀, y₀, z₀):
 *
 *   r_i = √((x_i-x₀)² + (y_i-y₀)² + (z_i-z₀)²)
 *
 *   e_i = [(x_i-x₀)/r_i, (y_i-y₀)/r_i, (z_i-z₀)/r_i] = LOS unit vector
 *
 *   H[i,:] = [-e_i_x, -e_i_y, -e_i_z, 1.0]
 *
 * The -e_i orientation is because:
 *   ΔP_i = P_corr,i - (r_i + b₀) ≈ -e_i·Δx_rx + 1·Δb
 * ------------------------------------------------------------------------- */

int gnss_design_matrix_build(const gnss_sat_meas_t *meas, int n_sats,
                               gnss_ecef_t rx_guess,
                               gnss_design_matrix_t *H) {
    if (n_sats < 4 || !meas || !H) return -1;

    H->n_sats = n_sats;
    H->n_params = 4;
    H->data = (double*)malloc(n_sats * 4 * sizeof(double));
    if (!H->data) return -2;

    int i;
    for (i = 0; i < n_sats; i++) {
        double dx = meas[i].sat_pos.x - rx_guess.x;
        double dy = meas[i].sat_pos.y - rx_guess.y;
        double dz = meas[i].sat_pos.z - rx_guess.z;
        double range = sqrt(dx*dx + dy*dy + dz*dz);

        double inv_range = (range > 1e-10) ? (1.0 / range) : 0.0;
        double e_x = dx * inv_range;
        double e_y = dy * inv_range;
        double e_z = dz * inv_range;

        H->data[i*4 + 0] = -e_x;
        H->data[i*4 + 1] = -e_y;
        H->data[i*4 + 2] = -e_z;
        H->data[i*4 + 3] = 1.0;
    }

    return 0;
}

void gnss_design_matrix_free(gnss_design_matrix_t *H) {
    if (H && H->data) {
        free(H->data);
        H->data = NULL;
    }
}

/* -------------------------------------------------------------------------
 * L3: Solve normal equations Δx = (Hᵀ·W·H)⁻¹·Hᵀ·W·ΔP
 *
 * Uses Cholesky decomposition for the 4×4 symmetric positive-definite
 * normal matrix N = Hᵀ·W·H.
 * ------------------------------------------------------------------------- */

int gnss_normal_eqn_solve(const gnss_design_matrix_t *H,
                           const double residuals[], int n_sats,
                           const double weights[],
                           double delta_x[4]) {
    if (!H || !residuals || n_sats < 4) return -1;

    double N[4][4] = {{0}};
    double b[4] = {0};
    int i, j, k;

    /* N = Hᵀ·W·H, b = Hᵀ·W·r */
    for (i = 0; i < n_sats; i++) {
        double w = (weights) ? weights[i] : 1.0;
        double Hi[4];
        Hi[0] = H->data[i*4+0];
        Hi[1] = H->data[i*4+1];
        Hi[2] = H->data[i*4+2];
        Hi[3] = H->data[i*4+3];

        for (j = 0; j < 4; j++) {
            for (k = 0; k < 4; k++) {
                N[j][k] += w * Hi[j] * Hi[k];
            }
            b[j] += w * Hi[j] * residuals[i];
        }
    }

    /* Solve N·Δx = b via Cholesky (N is 4×4 SPD) */
    double L[4][4] = {{0}};
    for (i = 0; i < 4; i++) {
        for (j = 0; j <= i; j++) {
            double sum = N[i][j];
            for (k = 0; k < j; k++)
                sum -= L[i][k] * L[j][k];
            if (i == j) {
                if (sum <= 1e-30) return -1;
                L[i][j] = sqrt(sum);
            } else {
                L[i][j] = sum / L[j][j];
            }
        }
    }

    /* Forward substitution: L·y = b */
    double y[4];
    for (i = 0; i < 4; i++) {
        double sum = b[i];
        for (j = 0; j < i; j++) sum -= L[i][j] * y[j];
        y[i] = sum / L[i][i];
    }

    /* Back substitution: Lᵀ·Δx = y */
    for (i = 3; i >= 0; i--) {
        double sum = y[i];
        for (j = i+1; j < 4; j++) sum -= L[j][i] * delta_x[j];
        delta_x[i] = sum / L[i][i];
    }

    return 0;
}

/* -------------------------------------------------------------------------
 * L5: Iterative least squares PVT solution
 * ------------------------------------------------------------------------- */

int gnss_ls_position_solve(const gnss_sat_meas_t *meas, int n_sats,
                            gnss_ecef_t initial_pos,
                            int max_iter, double tol,
                            gnss_pvt_solution_t *solution) {
    if (n_sats < 4 || !meas || !solution) return -1;

    memset(solution, 0, sizeof(gnss_pvt_solution_t));
    gnss_ecef_t pos = initial_pos;
    double clock_bias_m = 0.0;

    int iter;
    for (iter = 0; iter < max_iter; iter++) {
        /* Compute residuals */
        double *residuals = (double*)malloc(n_sats * sizeof(double));
        if (!residuals) return -2;

        int i;
        for (i = 0; i < n_sats; i++) {
            double range = gnss_geometric_range(meas[i].sat_pos, pos);
            residuals[i] = meas[i].pseudorange - (range + clock_bias_m);
        }

        /* Build design matrix */
        gnss_design_matrix_t H;
        if (gnss_design_matrix_build(meas, n_sats, pos, &H) != 0) {
            free(residuals);
            return -3;
        }

        /* Solve */
        double delta_x[4];
        int solve_rc = gnss_normal_eqn_solve(&H, residuals, n_sats, NULL, delta_x);
        gnss_design_matrix_free(&H);

        if (solve_rc != 0) {
            free(residuals);
            return solve_rc;
        }

        /* Update */
        pos.x += delta_x[0];
        pos.y += delta_x[1];
        pos.z += delta_x[2];
        clock_bias_m += delta_x[3];

        double delta_norm = sqrt(delta_x[0]*delta_x[0] + delta_x[1]*delta_x[1]
                               + delta_x[2]*delta_x[2] + delta_x[3]*delta_x[3]);

        /* Compute RMS residual */
        double sum_sq = 0.0;
        for (i = 0; i < n_sats; i++) {
            double range = gnss_geometric_range(meas[i].sat_pos, pos);
            double r = meas[i].pseudorange - (range + clock_bias_m);
            sum_sq += r * r;
        }
        double rms = sqrt(sum_sq / (double)n_sats);
        free(residuals);

        /* Check convergence */
        if (delta_norm < tol) {
            solution->pos = pos;
            solution->clock_bias = clock_bias_m;
            solution->n_sats = n_sats;
            solution->rms_residual = rms;
            solution->valid = 1;
            solution->iterations = iter + 1;

            /* Compute DOP */
            gnss_design_matrix_t H2;
            if (gnss_design_matrix_build(meas, n_sats, pos, &H2) == 0) {
                gnss_compute_dop_from_H(&H2, &solution->gdop, &solution->pdop,
                                         &solution->hdop, &solution->vdop,
                                         &solution->tdop);
                gnss_design_matrix_free(&H2);
            }
            return 0;
        }
    }

    /* Did not converge within max_iter */
    solution->pos = pos;
    solution->clock_bias = clock_bias_m;
    solution->n_sats = n_sats;
    solution->valid = 0;
    solution->iterations = max_iter;
    return -4;
}

/* -------------------------------------------------------------------------
 * L5: Weighted least squares PVT solution
 * ------------------------------------------------------------------------- */

int gnss_wls_position_solve(const gnss_sat_meas_t *meas, int n_sats,
                             gnss_ecef_t initial_pos,
                             gnss_weight_strategy_t wgt_strategy,
                             int max_iter, double tol,
                             gnss_pvt_solution_t *solution) {
    if (n_sats < 4 || !meas || !solution) return -1;

    memset(solution, 0, sizeof(gnss_pvt_solution_t));
    gnss_ecef_t pos = initial_pos;
    double clock_bias_m = 0.0;

    /* Allocate weights */
    double *weights = (double*)malloc(n_sats * sizeof(double));
    if (!weights) return -2;
    int i;
    for (i = 0; i < n_sats; i++) weights[i] = 1.0;

    int iter;
    for (iter = 0; iter < max_iter; iter++) {
        /* Update weights based on current geometry */
        for (i = 0; i < n_sats; i++) {
            switch (wgt_strategy) {
            case GNSS_WGT_UNIFORM:
                weights[i] = 1.0;
                break;
            case GNSS_WGT_ELEVATION: {
                double el = meas[i].elevation;
                if (el < 0.0) el = 0.0;
                double sin_el = sin(el);
                weights[i] = sin_el * sin_el + 0.01; /* floor */
                break;
            }
            case GNSS_WGT_SNR: {
                double snr_linear = pow(10.0, meas[i].weight / 10.0);
                weights[i] = (snr_linear > 0.01) ? snr_linear : 0.01;
                break;
            }
            default:
                weights[i] = 1.0;
                break;
            }
        }

        /* Compute residuals */
        double *residuals = (double*)malloc(n_sats * sizeof(double));
        if (!residuals) { free(weights); return -2; }

        for (i = 0; i < n_sats; i++) {
            double range = gnss_geometric_range(meas[i].sat_pos, pos);
            residuals[i] = meas[i].pseudorange - (range + clock_bias_m);
        }

        gnss_design_matrix_t H;
        if (gnss_design_matrix_build(meas, n_sats, pos, &H) != 0) {
            free(residuals); free(weights); return -3;
        }

        double delta_x[4];
        int solve_rc = gnss_normal_eqn_solve(&H, residuals, n_sats, weights, delta_x);
        gnss_design_matrix_free(&H);

        if (solve_rc != 0) {
            free(residuals); free(weights); return solve_rc;
        }

        pos.x += delta_x[0]; pos.y += delta_x[1]; pos.z += delta_x[2];
        clock_bias_m += delta_x[3];

        double delta_norm = sqrt(delta_x[0]*delta_x[0] + delta_x[1]*delta_x[1]
                               + delta_x[2]*delta_x[2] + delta_x[3]*delta_x[3]);

        double sum_sq = 0.0;
        for (i = 0; i < n_sats; i++) {
            double range = gnss_geometric_range(meas[i].sat_pos, pos);
            double r = meas[i].pseudorange - (range + clock_bias_m);
            sum_sq += weights[i] * r * r;
        }
        double rms = sqrt(sum_sq / (double)n_sats);
        free(residuals);

        if (delta_norm < tol) {
            solution->pos = pos;
            solution->clock_bias = clock_bias_m;
            solution->n_sats = n_sats;
            solution->rms_residual = rms;
            solution->valid = 1;
            solution->iterations = iter + 1;

            gnss_design_matrix_t H2;
            if (gnss_design_matrix_build(meas, n_sats, pos, &H2) == 0) {
                gnss_compute_dop_from_H(&H2, &solution->gdop, &solution->pdop,
                                         &solution->hdop, &solution->vdop,
                                         &solution->tdop);
                gnss_design_matrix_free(&H2);
            }
            free(weights);
            return 0;
        }
    }

    free(weights);
    solution->pos = pos;
    solution->clock_bias = clock_bias_m;
    solution->n_sats = n_sats;
    solution->valid = 0;
    solution->iterations = max_iter;
    return -4;
}

/* -------------------------------------------------------------------------
 * L3: DOP computation from design matrix
 *
 * Q = (HᵀH)⁻¹
 * GDOP = √(q₁₁+q₂₂+q₃₃+q₄₄)
 * PDOP = √(q₁₁+q₂₂+q₃₃)
 * TDOP = √(q₄₄)
 * HDOP/VDOP require Q rotation to ENU frame
 * ------------------------------------------------------------------------- */

int gnss_compute_dop_from_H(const gnss_design_matrix_t *H,
                             double *gdop, double *pdop,
                             double *hdop, double *vdop, double *tdop) {
    if (!H || H->n_sats < 4) return -1;

    /* Compute HᵀH */
    double HTH[4][4] = {{0}};
    int i, j, k;
    for (i = 0; i < H->n_sats; i++) {
        for (j = 0; j < 4; j++) {
            for (k = 0; k < 4; k++) {
                HTH[j][k] += H->data[i*4+j] * H->data[i*4+k];
            }
        }
    }

    /* Invert HTH → Q using Cholesky + back-substitution */
    double L[4][4] = {{0}};
    for (i = 0; i < 4; i++) {
        for (j = 0; j <= i; j++) {
            double sum = HTH[i][j];
            for (k = 0; k < j; k++) sum -= L[i][k] * L[j][k];
            if (i == j) {
                if (sum <= 1e-30) return -1;
                L[i][j] = sqrt(sum);
            } else {
                L[i][j] = sum / L[j][j];
            }
        }
    }

    /* Solve L·Lᵀ·Q = I → Q columns by forward/back substitution */
    double Q[4][4];
    int col;
    for (col = 0; col < 4; col++) {
        double e[4] = {(col==0)?1.0:0.0, (col==1)?1.0:0.0,
                        (col==2)?1.0:0.0, (col==3)?1.0:0.0};
        double y[4];
        /* Forward: L·y = e */
        for (i = 0; i < 4; i++) {
            double sum = e[i];
            for (j = 0; j < i; j++) sum -= L[i][j] * y[j];
            y[i] = sum / L[i][i];
        }
        /* Backward: Lᵀ·q = y */
        for (i = 3; i >= 0; i--) {
            double sum = y[i];
            for (j = i+1; j < 4; j++) sum -= L[j][i] * Q[j][col];
            Q[i][col] = sum / L[i][i];
        }
    }

    if (gdop) *gdop = sqrt(Q[0][0] + Q[1][1] + Q[2][2] + Q[3][3]);
    if (pdop) *pdop = sqrt(Q[0][0] + Q[1][1] + Q[2][2]);
    if (tdop) *tdop = sqrt(Q[3][3]);

    /* For HDOP/VDOP, rotate 3×3 position sub-matrix to ENU first.
     * Use a simplified estimate: HDOP ≈ √(PDOP² - VDOP²), VDOP ≈ √Q₂₂ in ENU.
     * Since we don't have the LLA here, approximate:
     *   HDOP ≈ √(Q₁₁+Q₂₂)/2·1.5  (crude but demonstrates the concept) */
    if (hdop || vdop) {
        double x_local = Q[0][0]*0.4 + Q[1][1]*0.4 + Q[2][2]*0.2;
        double z_local = Q[0][0]*0.2 + Q[1][1]*0.2 + Q[2][2]*0.6;
        if (hdop) *hdop = sqrt(x_local);
        if (vdop) *vdop = sqrt(z_local);
    }

    return 0;
}

/* -------------------------------------------------------------------------
 * L6: RAIM Fault Detection and Exclusion (residual-based)
 * ------------------------------------------------------------------------- */

int gnss_raim_fde(const gnss_sat_meas_t *meas, int n_sats,
                   gnss_ecef_t rx_guess, double sigma_threshold,
                   int excluded_mask[]) {
    if (n_sats < 5 || !meas || !excluded_mask) return -1;

    /* RAIM needs at least 5 satellites (4 for solution + 1 for detection) */
    int i;
    for (i = 0; i < n_sats; i++) excluded_mask[i] = 0;

    /* Compute LS solution with all satellites */
    gnss_design_matrix_t H;
    gnss_design_matrix_build(meas, n_sats, rx_guess, &H);

    double *residuals = (double*)malloc(n_sats * sizeof(double));
    if (!residuals) { gnss_design_matrix_free(&H); return -2; }

    /* Compute HᵀH */
    double HTH[4][4] = {{0}};
    int j, k;
    for (i = 0; i < n_sats; i++) {
        for (j = 0; j < 4; j++)
            for (k = 0; k < 4; k++)
                HTH[j][k] += H.data[i*4+j] * H.data[i*4+k];
    }

    /* Approximate Q = (HᵀH)⁻¹ using Cholesky */
    double L[4][4] = {{0}};
    for (i = 0; i < 4; i++) {
        for (j = 0; j <= i; j++) {
            double sum = HTH[i][j];
            for (k = 0; k < j; k++) sum -= L[i][k] * L[j][k];
            if (i == j) {
                L[i][j] = (sum > 1e-30) ? sqrt(sum) : 0.0;
            } else {
                L[i][j] = (L[j][j] > 1e-30) ? (sum / L[j][j]) : 0.0;
            }
        }
    }

    /* Compute S = I - H(HᵀH)⁻¹Hᵀ (the residual sensitivity matrix diagonal) */
    for (i = 0; i < n_sats; i++) {
        /* Solve (HᵀH)·coeff = Hᵀ row i */
        double Hi[4];
        Hi[0] = H.data[i*4+0]; Hi[1] = H.data[i*4+1];
        Hi[2] = H.data[i*4+2]; Hi[3] = H.data[i*4+3];

        /* Forward: L·y = Hi */
        double y[4];
        int p;
        for (p = 0; p < 4; p++) {
            double sum = Hi[p];
            for (k = 0; k < p; k++) sum -= L[p][k] * y[k];
            y[p] = (L[p][p] > 1e-30) ? (sum / L[p][p]) : 0.0;
        }
        /* Backward: Lᵀ·coeff = y */
        double coeff[4];
        for (p = 3; p >= 0; p--) {
            double sum = y[p];
            for (k = p+1; k < 4; k++) sum -= L[k][p] * coeff[k];
            coeff[p] = (L[p][p] > 1e-30) ? (sum / L[p][p]) : 0.0;
        }

        double h_dot_coeff = Hi[0]*coeff[0]+Hi[1]*coeff[1]
                            +Hi[2]*coeff[2]+Hi[3]*coeff[3];
        double S_ii = 1.0 - h_dot_coeff;

        /* Standardized residual: r_i / (σ·√S_ii) */
        double std_res = (S_ii > 1e-10) ? fabs(residuals[i] / sqrt(S_ii)) : 0.0;
        if (std_res > sigma_threshold) {
            excluded_mask[i] = 1;
        }
    }

    free(residuals);
    gnss_design_matrix_free(&H);
    return 0;
}

/* -------------------------------------------------------------------------
 * L6: Post-fit RMS
 * ------------------------------------------------------------------------- */

double gnss_postfit_rms(const gnss_sat_meas_t *meas, int n_sats,
                         gnss_ecef_t rx_pos, double clock_bias) {
    int i;
    double sum_sq = 0.0;
    for (i = 0; i < n_sats; i++) {
        double range = gnss_geometric_range(meas[i].sat_pos, rx_pos);
        double residual = meas[i].pseudorange - (range + clock_bias);
        sum_sq += residual * residual;
    }
    return sqrt(sum_sq / (double)n_sats);
}
