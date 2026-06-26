/**
 * @file sensor_fusion.c
 * @brief Sensor fusion algorithms: Kalman, EKF, UKF, particle filter,
 *        complementary filter, matrix utilities
 */

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "../include/sensor_fusion.h"
#include "../include/indoor_positioning.h"

/* ============================================================================
 * L3 - Matrix Utility Functions
 * ============================================================================ */

void matrix_vec_mult(const double *A, const double *x, double *y,
                     int rows, int cols) {
    if (!A || !x || !y) return;
    for (int r = 0; r < rows; r++) {
        y[r] = 0.0;
        for (int c = 0; c < cols; c++) {
            y[r] += A[r * cols + c] * x[c];
        }
    }
}

void matrix_mult(const double *A, const double *B, double *C,
                 int m, int k, int n) {
    if (!A || !B || !C) return;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int p = 0; p < k; p++) {
                sum += A[i * k + p] * B[p * n + j];
            }
            C[i * n + j] = sum;
        }
    }
}

void matrix_transpose(const double *A, double *B, int rows, int cols) {
    if (!A || !B) return;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            B[j * rows + i] = A[i * cols + j];
        }
    }
}

/* ============================================================================
 * L3 - Cholesky Decomposition (Banachiewicz algorithm)
 * ============================================================================ */

int cholesky_decompose(double *A, int n) {
    if (!A || n <= 0) return -1;

    for (int j = 0; j < n; j++) {
        double sum = 0.0;
        for (int k = 0; k < j; k++) {
            sum += A[j * n + k] * A[j * n + k];
        }
        double diag = A[j * n + j] - sum;
        if (diag <= 1e-15) {
            return -1;  /* Not positive definite */
        }
        A[j * n + j] = sqrt(diag);

        for (int i = j + 1; i < n; i++) {
            sum = 0.0;
            for (int k = 0; k < j; k++) {
                sum += A[i * n + k] * A[j * n + k];
            }
            A[i * n + j] = (A[i * n + j] - sum) / A[j * n + j];
        }
    }
    return 0;
}

void cholesky_solve(const double *L, const double *b, double *x, int n) {
    if (!L || !b || !x || n <= 0) return;

    /* Forward substitution: L * y = b */
    double y[KF_STATE_MAX_DIM];
    for (int i = 0; i < n; i++) {
        double sum = b[i];
        for (int j = 0; j < i; j++) {
            sum -= L[i * n + j] * y[j];
        }
        y[i] = sum / L[i * n + i];
    }

    /* Back substitution: L^T * x = y */
    for (int i = n - 1; i >= 0; i--) {
        double sum = y[i];
        for (int j = i + 1; j < n; j++) {
            sum -= L[j * n + i] * x[j];
        }
        x[i] = sum / L[i * n + i];
    }
}

/* ============================================================================
 * L3 - Mahalanobis Distance
 * ============================================================================ */

double mahalanobis_distance(const double *x, const double *mean,
                            const double *S_inv, int n) {
    if (!x || !mean || !S_inv || n <= 0) return 0.0;
    double diff[KF_STATE_MAX_DIM];
    for (int i = 0; i < n && i < KF_STATE_MAX_DIM; i++) {
        diff[i] = x[i] - mean[i];
    }
    /* Compute diff^T * S^{-1} * diff */
    double temp[KF_STATE_MAX_DIM];
    matrix_vec_mult(S_inv, diff, temp, n, n);
    double d2 = 0.0;
    for (int i = 0; i < n; i++) {
        d2 += diff[i] * temp[i];
    }
    return d2;
}

/* ============================================================================
 * L5 - Linear Kalman Filter
 * ============================================================================ */

void kf_init(kalman_filter_t *kf, int n_state, int n_measure) {
    if (!kf) return;
    if (n_state > KF_STATE_MAX_DIM) n_state = KF_STATE_MAX_DIM;
    if (n_measure > KF_STATE_MAX_DIM) n_measure = KF_STATE_MAX_DIM;
    memset(kf, 0, sizeof(kalman_filter_t));
    kf->n_state = n_state;
    kf->n_measure = n_measure;
    /* Initialize F as identity */
    for (int i = 0; i < n_state; i++) kf->F[i][i] = 1.0;
    /* Initialize Q and R as identity */
    for (int i = 0; i < n_state; i++) kf->Q[i][i] = 1.0;
    for (int i = 0; i < n_measure; i++) kf->R[i][i] = 1.0;
    /* Initialize P with large uncertainty */
    for (int i = 0; i < n_state; i++) kf->P[i][i] = 1000.0;
}

void kf_set_F(kalman_filter_t *kf, const double *F) {
    if (!kf || !F) return;
    for (int r = 0; r < kf->n_state; r++)
        for (int c = 0; c < kf->n_state; c++)
            kf->F[r][c] = F[r * kf->n_state + c];
}

void kf_set_H(kalman_filter_t *kf, const double *H) {
    if (!kf || !H) return;
    for (int r = 0; r < kf->n_measure; r++)
        for (int c = 0; c < kf->n_state; c++)
            kf->H[r][c] = H[r * kf->n_state + c];
}

void kf_set_Q(kalman_filter_t *kf, const double *Q) {
    if (!kf || !Q) return;
    for (int r = 0; r < kf->n_state; r++)
        for (int c = 0; c < kf->n_state; c++)
            kf->Q[r][c] = Q[r * kf->n_state + c];
}

void kf_set_R(kalman_filter_t *kf, const double *R) {
    if (!kf || !R) return;
    for (int r = 0; r < kf->n_measure; r++)
        for (int c = 0; c < kf->n_measure; c++)
            kf->R[r][c] = R[r * kf->n_measure + c];
}

void kf_set_initial(kalman_filter_t *kf, const double *x0, const double *P0) {
    if (!kf || !x0 || !P0) return;
    for (int i = 0; i < kf->n_state; i++) kf->x[i] = x0[i];
    for (int r = 0; r < kf->n_state; r++)
        for (int c = 0; c < kf->n_state; c++)
            kf->P[r][c] = P0[r * kf->n_state + c];
    kf->initialized = 1;
}

/* Internal: solve linear system via Cholesky for symmetric positive definite matrix */
static int solve_sym_posdef(const double *A, const double *b, double *x, int n) {
    double L[KF_STATE_MAX_DIM * KF_STATE_MAX_DIM];
    memcpy(L, A, n * n * sizeof(double));
    if (cholesky_decompose(L, n) != 0) return -1;
    cholesky_solve(L, b, x, n);
    return 0;
}

void kf_predict(kalman_filter_t *kf) {
    if (!kf || !kf->initialized) return;
    int n = kf->n_state;

    /* x = F * x */
    double x_pred[KF_STATE_MAX_DIM];
    matrix_vec_mult((const double *)kf->F, kf->x, x_pred, n, n);
    memcpy(kf->x, x_pred, n * sizeof(double));

    /* P = F * P * F^T + Q */
    double FPFt[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM];
    double FP[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM];
    double Ft[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM];

    matrix_mult((const double *)kf->F, (const double *)kf->P,
                (double *)FP, n, n, n);
    matrix_transpose((const double *)kf->F, (double *)Ft, n, n);
    matrix_mult((const double *)FP, (const double *)Ft,
                (double *)FPFt, n, n, n);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            kf->P[i][j] = FPFt[i][j] + kf->Q[i][j];
}

void kf_update(kalman_filter_t *kf, const double *z) {
    if (!kf || !z || !kf->initialized) return;
    int n = kf->n_state;
    int m = kf->n_measure;

    /* Innovation: y = z - H*x */
    double y[KF_STATE_MAX_DIM];
    double Hx[KF_STATE_MAX_DIM];
    matrix_vec_mult((const double *)kf->H, kf->x, Hx, m, n);
    for (int i = 0; i < m; i++) y[i] = z[i] - Hx[i];

    /* Innovation covariance: S = H*P*H^T + R */
    double PHt[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM];
    double Ht[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM];
    double S[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM];

    matrix_transpose((const double *)kf->H, (double *)Ht, m, n);
    matrix_mult((const double *)kf->P, (const double *)Ht,
                (double *)PHt, n, n, m);
    matrix_mult((const double *)kf->H, (const double *)PHt,
                (double *)S, m, n, m);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
            S[i][j] += kf->R[i][j];

    /* Kalman gain: K = P*H^T * S^{-1} */
    /* Solve S*Kt = PHt for Kt (transpose of K), so Kt = S^{-1} * PHt */
    double Kt[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM];
    for (int col = 0; col < n; col++) {
        double rhs[KF_STATE_MAX_DIM];
        for (int i = 0; i < m; i++) {
            rhs[i] = PHt[col][i];  /* (PHt)^T = H*P */
        }
        double sol[KF_STATE_MAX_DIM];
        /* Need to solve S * sol = rhs. But we have PHt = P*H^T with dims n x m.
         * Let's use the correct formula: K = P*H^T * S^{-1}
         * So K is n x m. We compute S^{-1} * (H*P)^T .
         * Actually let me compute: K = PHt * S^{-1}
         * Solve S^T * K^T = PHt^T
         * Since S is symmetric, solve S * K_row = PHt_row */
        if (solve_sym_posdef((const double *)S, rhs, sol, m) == 0) {
            for (int i = 0; i < m; i++) Kt[col][i] = sol[i];
        } else {
            /* Fallback: use simplified gain K = P*Ht / R (scalar approx) */
            for (int i = 0; i < m; i++) Kt[col][i] = 0.0;
            /* Use diagonal S for simplicity in fallback */
            for (int i = 0; i < m; i++)
                if (S[i][i] > 1e-12) Kt[col][i] = PHt[col][i] / S[i][i];
        }
    }

    /* State update: x = x + K * y
     * K is n x m: K[i][j] = Kt[j][i] */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            kf->x[i] += Kt[i][j] * y[j];
        }
    }

    /* Covariance update: P = (I - K*H) * P */
    double KH[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM];
    double I_KH[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM];
    double P_new[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM];

    /* Compute K * H: n x n */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int p = 0; p < m; p++) {
                sum += Kt[i][p] * kf->H[p][j];
            }
            KH[i][j] = sum;
        }
    }

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            I_KH[i][j] = (i == j) ? 1.0 - KH[i][j] : -KH[i][j];

    matrix_mult((const double *)I_KH, (const double *)kf->P,
                (double *)P_new, n, n, n);

    /* Joseph form for numerical stability: P = (I-KH)*P*(I-KH)^T + K*R*K^T */
    double temp1[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM];
    double I_KHt[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM];
    matrix_transpose((const double *)I_KH, (double *)I_KHt, n, n);
    matrix_mult((const double *)P_new, (const double *)I_KHt,
                (double *)temp1, n, n, n);

    /* Add K*R*K^T */
    double KR[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM];
    double KRKt[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            double sum = 0.0;
            for (int p = 0; p < m; p++) sum += Kt[i][p] * kf->R[p][j];
            KR[i][j] = sum;
        }
    }
    /* Kt is n x m, so K^T is m x n = transpose of Kt
     * Need K^T: K^T[j][k] = Kt[k][j]
     * KRKt = (n x m) * (m x n) -> n x n */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int p = 0; p < m; p++) sum += KR[i][p] * Kt[j][p];  /* Kt[j][p] = K^T[p][j] */
            KRKt[i][j] = sum;
        }
    }

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            kf->P[i][j] = temp1[i][j] + KRKt[i][j];
}

void kf_get_state(const kalman_filter_t *kf, double *x) {
    if (!kf || !x) return;
    memcpy(x, kf->x, kf->n_state * sizeof(double));
}

void kf_get_covariance_diag(const kalman_filter_t *kf, double *diag) {
    if (!kf || !diag) return;
    for (int i = 0; i < kf->n_state; i++) diag[i] = kf->P[i][i];
}

/* ============================================================================
 * L5 - Extended Kalman Filter
 * ============================================================================ */

void ekf_init(ekf_t *ekf, int n_state, int n_measure) {
    if (!ekf) return;
    if (n_state > KF_STATE_MAX_DIM) n_state = KF_STATE_MAX_DIM;
    if (n_measure > KF_STATE_MAX_DIM) n_measure = KF_STATE_MAX_DIM;
    memset(ekf, 0, sizeof(ekf_t));
    ekf->n_state = n_state;
    ekf->n_measure = n_measure;
    for (int i = 0; i < n_state; i++) ekf->P[i][i] = 1000.0;
    for (int i = 0; i < n_state; i++) ekf->Q[i][i] = 0.01;
    for (int i = 0; i < n_measure; i++) ekf->R[i][i] = 1.0;
}

void ekf_predict(ekf_t *ekf, ekf_transition_fn f, ekf_jacobian_fn F_jacobian,
                 double dt, void *user_data) {
    if (!ekf || !f || !F_jacobian) return;
    int n = ekf->n_state;

    /* Compute predicted state: x_pred = f(x, dt) */
    double x_pred[KF_STATE_MAX_DIM];
    f(ekf->x, dt, user_data, x_pred);

    /* Compute Jacobian: F = ∂f/∂x at current state */
    double F_mat[KF_STATE_MAX_DIM * KF_STATE_MAX_DIM];
    F_jacobian(ekf->x, dt, user_data, F_mat);

    /* P = F * P * F^T + Q */
    double FP[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM];
    double Ft[KF_STATE_MAX_DIM * KF_STATE_MAX_DIM];
    double FPFt[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM];

    matrix_mult(F_mat, (const double *)ekf->P, (double *)FP, n, n, n);
    matrix_transpose(F_mat, Ft, n, n);
    matrix_mult((const double *)FP, Ft, (double *)FPFt, n, n, n);

    memcpy(ekf->x, x_pred, n * sizeof(double));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            ekf->P[i][j] = FPFt[i][j] + ekf->Q[i][j];
}

void ekf_update(ekf_t *ekf, const double *z,
                ekf_measurement_fn h, ekf_jacobian_fn H_jacobian,
                void *user_data) {
    if (!ekf || !z || !h || !H_jacobian) return;
    int n = ekf->n_state;
    int m = ekf->n_measure;

    /* Predicted measurement */
    double z_pred[KF_STATE_MAX_DIM];
    h(ekf->x, user_data, z_pred);

    /* Innovation */
    double y[KF_STATE_MAX_DIM];
    for (int i = 0; i < m; i++) y[i] = z[i] - z_pred[i];

    /* Measurement Jacobian */
    double H_mat[KF_STATE_MAX_DIM * KF_STATE_MAX_DIM];
    H_jacobian(ekf->x, 0.0, user_data, H_mat);

    /* Innovation covariance: S = H*P*H^T + R */
    double PHt[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM];
    double Ht[KF_STATE_MAX_DIM * KF_STATE_MAX_DIM];
    double S[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM];

    matrix_transpose(H_mat, Ht, m, n);
    matrix_mult((const double *)ekf->P, Ht, (double *)PHt, n, n, m);
    matrix_mult(H_mat, (const double *)PHt, (double *)S, m, n, m);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
            S[i][j] += ekf->R[i][j];

    /* Kalman gain: K = P*H^T * S^{-1} */
    double Kt[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM];
    for (int col = 0; col < n; col++) {
        double rhs[KF_STATE_MAX_DIM];
        for (int i = 0; i < m; i++) rhs[i] = PHt[col][i];
        double sol[KF_STATE_MAX_DIM];
        if (solve_sym_posdef((const double *)S, rhs, sol, m) == 0) {
            for (int i = 0; i < m; i++) Kt[col][i] = sol[i];
        } else {
            for (int i = 0; i < m; i++) Kt[col][i] = 0.0;
            for (int i = 0; i < m; i++)
                if (S[i][i] > 1e-12) Kt[col][i] = PHt[col][i] / S[i][i];
        }
    }

    /* State update */
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            ekf->x[i] += Kt[i][j] * y[j];

    /* Covariance: Joseph form */
    double KH[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int p = 0; p < m; p++) sum += Kt[i][p] * H_mat[p * n + j];
            KH[i][j] = sum;
        }

    double I_KH[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            I_KH[i][j] = (i == j) ? 1.0 - KH[i][j] : -KH[i][j];

    double temp[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM];
    matrix_mult((const double *)I_KH, (const double *)ekf->P, (double *)temp, n, n, n);

    double I_KHt[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM];
    matrix_transpose((const double *)I_KH, (double *)I_KHt, n, n);

    double P_joseph[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM];
    matrix_mult((const double *)temp, (const double *)I_KHt, (double *)P_joseph, n, n, n);

    /* K*R*K^T */
    double KR[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++) {
            double sum = 0.0;
            for (int p = 0; p < m; p++) sum += Kt[i][p] * ekf->R[p][j];
            KR[i][j] = sum;
        }

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int p = 0; p < m; p++) sum += KR[i][p] * Kt[j][p];
            ekf->P[i][j] = P_joseph[i][j] + sum;
        }

    ekf->initialized = 1;
}

void ekf_get_state(const ekf_t *ekf, double *x) {
    if (!ekf || !x) return;
    memcpy(x, ekf->x, ekf->n_state * sizeof(double));
}

/* ============================================================================
 * L6 - EKF Setup: Constant Velocity Ranging
 * ============================================================================ */

/* State: [x, y, vx, vy] */
void cv_transition_fn(const double *x, double dt, void *user_data,
                      double *x_pred) {
    (void)user_data;
    x_pred[0] = x[0] + x[2] * dt;  /* x = x + vx*dt */
    x_pred[1] = x[1] + x[3] * dt;  /* y = y + vy*dt */
    x_pred[2] = x[2];               /* vx constant */
    x_pred[3] = x[3];               /* vy constant */
}

void cv_jacobian_fn(const double *x, double dt, void *user_data,
                    double *F) {
    (void)x; (void)user_data;
    /* F = [[1,0,dt,0],[0,1,0,dt],[0,0,1,0],[0,0,0,1]] */
    memset(F, 0, 16 * sizeof(double));
    F[0] = 1.0; F[1] = 0.0; F[2] = dt;  F[3] = 0.0;
    F[4] = 0.0; F[5] = 1.0; F[6] = 0.0; F[7] = dt;
    F[8] = 0.0; F[9] = 0.0; F[10] = 1.0; F[11] = 0.0;
    F[12] = 0.0; F[13] = 0.0; F[14] = 0.0; F[15] = 1.0;
}

void ekf_range_measurement_fn(const double *x, void *user_data, double *z_pred) {
    ekf_range_measurement_data_t *d = (ekf_range_measurement_data_t *)user_data;
    double dx = x[0] - d->anchor_x;
    double dy = x[1] - d->anchor_y;
    z_pred[0] = sqrt(dx*dx + dy*dy);
}

void ekf_range_jacobian_fn(const double *x, double dt, void *user_data,
                           double *H) {
    (void)dt;
    ekf_range_measurement_data_t *d = (ekf_range_measurement_data_t *)user_data;
    double dx = x[0] - d->anchor_x;
    double dy = x[1] - d->anchor_y;
    double dist = sqrt(dx*dx + dy*dy);
    if (dist < 1e-10) dist = 1e-10;
    memset(H, 0, KF_STATE_MAX_DIM * KF_STATE_MAX_DIM * sizeof(double));
    H[0] = dx / dist;  /* ∂h/∂x */
    H[1] = dy / dist;  /* ∂h/∂y */
    /* ∂h/∂vx = 0, ∂h/∂vy = 0 */
}

void ekf_setup_constant_velocity_ranging(ekf_t *ekf, double dt,
                                         double process_noise_accel,
                                         double measurement_noise_dist,
                                         double initial_pos_x,
                                         double initial_pos_y) {
    if (!ekf) return;
    ekf_init(ekf, 4, 1);

    /* Process noise (discrete white noise acceleration model):
     * Q = G * sigma_a^2 * G^T
     * G = [dt^2/2, 0; 0, dt^2/2; dt, 0; 0, dt] */
    double q = process_noise_accel * process_noise_accel;
    double dt2 = dt * dt;
    double dt3_2 = dt2 * dt / 2.0;
    double dt4_4 = dt2 * dt2 / 4.0;

    double Q_flat[16] = {
        dt4_4 * q, 0,        dt3_2 * q, 0,
        0,        dt4_4 * q, 0,        dt3_2 * q,
        dt3_2 * q, 0,        dt2 * q,  0,
        0,        dt3_2 * q, 0,        dt2 * q
    };
    for (int r = 0; r < 4; r++)
        for (int c = 0; c < 4; c++)
            ekf->Q[r][c] = Q_flat[r * 4 + c];

    ekf->R[0][0] = measurement_noise_dist * measurement_noise_dist;

    ekf->x[0] = initial_pos_x;
    ekf->x[1] = initial_pos_y;

    double P0[16] = {
        100.0, 0, 0, 0,  0, 100.0, 0, 0,
        0, 0, 10.0, 0,   0, 0, 0, 10.0
    };
    for (int r = 0; r < 4; r++)
        for (int c = 0; c < 4; c++)
            ekf->P[r][c] = P0[r * 4 + c];

    ekf->initialized = 1;
}

/* ============================================================================
 * L8 Advanced - Unscented Kalman Filter
 * ============================================================================ */

void ukf_init(ukf_t *ukf, int n_state, int n_measure,
              double alpha, double beta, double kappa) {
    if (!ukf) return;
    if (n_state > KF_STATE_MAX_DIM) n_state = KF_STATE_MAX_DIM;
    if (n_measure > KF_STATE_MAX_DIM) n_measure = KF_STATE_MAX_DIM;
    memset(ukf, 0, sizeof(ukf_t));
    ukf->n_state = n_state;
    ukf->n_measure = n_measure;
    ukf->alpha = alpha;
    ukf->beta = beta;
    ukf->kappa = kappa;
    ukf->lambda = alpha * alpha * (n_state + kappa) - n_state;
    for (int i = 0; i < n_state; i++) ukf->P[i][i] = 1000.0;
    for (int i = 0; i < n_state; i++) ukf->Q[i][i] = 0.01;
    for (int i = 0; i < n_measure; i++) ukf->R[i][i] = 1.0;
}

void ukf_predict(ukf_t *ukf, ekf_transition_fn f, double dt, void *user_data) {
    if (!ukf || !f) return;
    int n = ukf->n_state;
    double lambda = ukf->lambda;

    /* Generate sigma points (2*n+1) */
    int n_sigma = 2 * n + 1;

    /* Cholesky decomposition of P */
    double sqrtP[KF_STATE_MAX_DIM * KF_STATE_MAX_DIM];
    memcpy(sqrtP, ukf->P, n * n * sizeof(double));
    if (cholesky_decompose(sqrtP, n) != 0) {
        /* Fallback: diagonal sqrt */
        memset(sqrtP, 0, n * n * sizeof(double));
        for (int i = 0; i < n; i++)
            sqrtP[i * n + i] = sqrt(fabs(ukf->P[i][i]));
    }
    /* Scale by sqrt(n + lambda) */
    double scale = sqrt(n + lambda);
    for (int i = 0; i < n * n; i++) sqrtP[i] *= scale;

    /* Sigma points */
    double sigma_x[UKF_SIGMA_PTS_MAX][KF_STATE_MAX_DIM];

    /* sigma_0 = x */
    for (int j = 0; j < n; j++) sigma_x[0][j] = ukf->x[j];

    /* sigma_i = x + sqrtP_column_i */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            sigma_x[i + 1][j] = ukf->x[j] + sqrtP[j * n + i];
        }
    }
    /* sigma_{i+n} = x - sqrtP_column_i */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            sigma_x[i + n + 1][j] = ukf->x[j] - sqrtP[j * n + i];
        }
    }

    /* Propagate sigma points */
    double sigma_x_pred[UKF_SIGMA_PTS_MAX][KF_STATE_MAX_DIM];
    for (int i = 0; i < n_sigma; i++) {
        f(sigma_x[i], dt, user_data, sigma_x_pred[i]);
    }

    /* Weights */
    double Wm[UKF_SIGMA_PTS_MAX];
    double Wc[UKF_SIGMA_PTS_MAX];
    Wm[0] = lambda / (n + lambda);
    Wc[0] = lambda / (n + lambda) + (1.0 - ukf->alpha * ukf->alpha + ukf->beta);
    for (int i = 1; i < n_sigma; i++) {
        Wm[i] = 0.5 / (n + lambda);
        Wc[i] = 0.5 / (n + lambda);
    }

    /* Predicted mean */
    memset(ukf->x, 0, n * sizeof(double));
    for (int i = 0; i < n_sigma; i++)
        for (int j = 0; j < n; j++)
            ukf->x[j] += Wm[i] * sigma_x_pred[i][j];

    /* Predicted covariance */
    memset(ukf->P, 0, sizeof(ukf->P));
    for (int i = 0; i < n_sigma; i++) {
        double diff[KF_STATE_MAX_DIM];
        for (int j = 0; j < n; j++) diff[j] = sigma_x_pred[i][j] - ukf->x[j];
        for (int r = 0; r < n; r++)
            for (int c = 0; c < n; c++)
                ukf->P[r][c] += Wc[i] * diff[r] * diff[c];
    }
    /* Add process noise */
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            ukf->P[i][j] += ukf->Q[i][j];

    ukf->initialized = 1;
}

void ukf_update(ukf_t *ukf, const double *z, ekf_measurement_fn h,
                void *user_data) {
    if (!ukf || !z || !h) return;
    int n = ukf->n_state;
    int m = ukf->n_measure;
    double lambda = ukf->lambda;
    int n_sigma = 2 * n + 1;

    /* Generate sigma points from current state */
    double sqrtP[KF_STATE_MAX_DIM * KF_STATE_MAX_DIM];
    memcpy(sqrtP, ukf->P, n * n * sizeof(double));
    cholesky_decompose(sqrtP, n);
    double scale = sqrt(n + lambda);
    for (int i = 0; i < n * n; i++) sqrtP[i] *= scale;

    double sigma_x[UKF_SIGMA_PTS_MAX][KF_STATE_MAX_DIM];
    for (int j = 0; j < n; j++) sigma_x[0][j] = ukf->x[j];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            sigma_x[i+1][j] = ukf->x[j] + sqrtP[j * n + i];
            sigma_x[i+n+1][j] = ukf->x[j] - sqrtP[j * n + i];
        }
    }

    /* Propagate sigma points through measurement function */
    double sigma_z[UKF_SIGMA_PTS_MAX][KF_STATE_MAX_DIM];
    for (int i = 0; i < n_sigma; i++) {
        h(sigma_x[i], user_data, sigma_z[i]);
    }

    /* Weights */
    double Wm[UKF_SIGMA_PTS_MAX], Wc[UKF_SIGMA_PTS_MAX];
    Wm[0] = lambda/(n+lambda);
    Wc[0] = lambda/(n+lambda) + (1.0 - ukf->alpha*ukf->alpha + ukf->beta);
    for (int i = 1; i < n_sigma; i++) Wm[i] = Wc[i] = 0.5/(n+lambda);

    /* Predicted measurement mean */
    double z_mean[KF_STATE_MAX_DIM] = {0};
    for (int i = 0; i < n_sigma; i++)
        for (int j = 0; j < m; j++)
            z_mean[j] += Wm[i] * sigma_z[i][j];

    /* Innovation covariance P_zz */
    double Pzz[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM] = {{0}};
    for (int i = 0; i < n_sigma; i++) {
        double diff[KF_STATE_MAX_DIM];
        for (int j = 0; j < m; j++) diff[j] = sigma_z[i][j] - z_mean[j];
        for (int r = 0; r < m; r++)
            for (int c = 0; c < m; c++)
                Pzz[r][c] += Wc[i] * diff[r] * diff[c];
    }
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
            Pzz[i][j] += ukf->R[i][j];

    /* Cross covariance P_xz */
    double Pxz[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM] = {{0}};
    for (int i = 0; i < n_sigma; i++) {
        double diff_x[KF_STATE_MAX_DIM];
        double diff_z[KF_STATE_MAX_DIM];
        for (int j = 0; j < n; j++) diff_x[j] = sigma_x[i][j] - ukf->x[j];
        for (int j = 0; j < m; j++) diff_z[j] = sigma_z[i][j] - z_mean[j];
        for (int r = 0; r < n; r++)
            for (int c = 0; c < m; c++)
                Pxz[r][c] += Wc[i] * diff_x[r] * diff_z[c];
    }

    /* Kalman gain K = Pxz * Pzz^{-1} */
    double K[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM] = {{0}};
    for (int col = 0; col < m; col++) {
        double rhs[KF_STATE_MAX_DIM];
        for (int i = 0; i < m; i++) rhs[i] = Pzz[i][col];
        double sol[KF_STATE_MAX_DIM];
        if (solve_sym_posdef((const double *)Pzz, rhs, sol, m) == 0) {
            for (int r = 0; r < n; r++) {
                double sum = 0.0;
                for (int p = 0; p < m; p++) sum += Pxz[r][p] * sol[p];
                K[r][col] = sum;
            }
        }
    }

    /* State update */
    double innovation[KF_STATE_MAX_DIM];
    for (int i = 0; i < m; i++) innovation[i] = z[i] - z_mean[i];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            ukf->x[i] += K[i][j] * innovation[j];

    /* Covariance update: P = P - K * Pzz * K^T */
    double KPzz[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM] = {{0}};
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            for (int p = 0; p < m; p++)
                KPzz[i][j] += K[i][p] * Pzz[p][j];

    double P_new[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM] = {{0}};
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int p = 0; p < m; p++) sum += KPzz[i][p] * K[j][p];
            P_new[i][j] = ukf->P[i][j] - sum;
        }
    memcpy(ukf->P, P_new, sizeof(ukf->P));
}

void ukf_get_state(const ukf_t *ukf, double *x) {
    if (!ukf || !x) return;
    memcpy(x, ukf->x, ukf->n_state * sizeof(double));
}

/* ============================================================================
 * L5 - Complementary Filters
 * ============================================================================ */

double complementary_filter_angle(double current_angle, double gyro_rate,
                                  double accel_angle, double alpha, double dt) {
    return alpha * (current_angle + gyro_rate * dt) + (1.0 - alpha) * accel_angle;
}

double complementary_filter_position(double fused_pos, double vel,
                                     double absolute_pos, double alpha, double dt) {
    double dead_reckoned = fused_pos + vel * dt;
    return alpha * dead_reckoned + (1.0 - alpha) * absolute_pos;
}

/* ============================================================================
 * L8 - Particle Filter
 * ============================================================================ */

void pf_init_uniform(particle_filter_t *pf, int n_particles, int n_state,
                     const double *x_min, const double *x_max) {
    if (!pf || !x_min || !x_max) return;
    if (n_particles > PF_MAX_PARTICLES) n_particles = PF_MAX_PARTICLES;
    if (n_state > KF_STATE_MAX_DIM) n_state = KF_STATE_MAX_DIM;

    memset(pf, 0, sizeof(particle_filter_t));
    pf->n_particles = n_particles;
    pf->n_state = n_state;
    pf->effective_n_thresh = 0.5;  /* Resample when N_eff < N/2 */

    /* Initialize particles uniformly */
    for (int i = 0; i < n_particles; i++) {
        for (int j = 0; j < n_state; j++) {
            /* Uniform random in [x_min[j], x_max[j]] */
            double range = x_max[j] - x_min[j];
            /* Simple LCG random number */
            double r = (double)((i * 1103515245 + 12345 + j * 2531011) & 0x7FFFFFFF)
                       / (double)0x7FFFFFFF;
            pf->particles[i].x[j] = x_min[j] + r * range;
        }
        pf->particles[i].weight = 1.0 / n_particles;
    }
    pf->initialized = 1;
}

void pf_predict(particle_filter_t *pf, ekf_transition_fn f,
                double dt, void *user_data) {
    if (!pf || !f || !pf->initialized) return;
    for (int i = 0; i < pf->n_particles; i++) {
        double x_new[KF_STATE_MAX_DIM];
        f(pf->particles[i].x, dt, user_data, x_new);
        memcpy(pf->particles[i].x, x_new, pf->n_state * sizeof(double));
    }
}

void pf_update_gaussian(particle_filter_t *pf, const double *z,
                        int n_measure, double measurement_std) {
    if (!pf || !z || !pf->initialized) return;

    double var = measurement_std * measurement_std;
    double weight_sum = 0.0;

    for (int i = 0; i < pf->n_particles; i++) {
        /* Gaussian likelihood: exp(-0.5 * ||x_pos - z||^2 / var) */
        double sum_sq = 0.0;
        for (int j = 0; j < n_measure && j < pf->n_state; j++) {
            double diff = pf->particles[i].x[j] - z[j];
            sum_sq += diff * diff;
        }
        double log_lik = -0.5 * sum_sq / var;
        double lik = exp(log_lik);
        pf->particles[i].weight *= lik;
        weight_sum += pf->particles[i].weight;
    }

    /* Normalize weights */
    if (weight_sum > 1e-30) {
        for (int i = 0; i < pf->n_particles; i++) {
            pf->particles[i].weight /= weight_sum;
        }
    } else {
        /* All weights zero: reinitialize uniform */
        for (int i = 0; i < pf->n_particles; i++) {
            pf->particles[i].weight = 1.0 / pf->n_particles;
        }
    }
}

void pf_resample(particle_filter_t *pf) {
    if (!pf || !pf->initialized) return;

    int N = pf->n_particles;

    /* Allocate on heap to avoid stack overflow */
    particle_t *new_particles = (particle_t *)malloc(N * sizeof(particle_t));
    double *cumulative = (double *)malloc(N * sizeof(double));
    if (!new_particles || !cumulative) {
        free(new_particles);
        free(cumulative);
        return;
    }

    /* Systematic resampling */
    cumulative[0] = pf->particles[0].weight;
    for (int i = 1; i < N; i++) {
        cumulative[i] = cumulative[i-1] + pf->particles[i].weight;
    }

    double u0 = (double)(rand()) / RAND_MAX / N;
    int j = 0;

    for (int i = 0; i < N; i++) {
        double u = u0 + (double)i / N;
        while (j < N - 1 && u > cumulative[j]) j++;
        memcpy(&new_particles[i], &pf->particles[j], sizeof(particle_t));
        new_particles[i].weight = 1.0 / N;
    }

    memcpy(pf->particles, new_particles, N * sizeof(particle_t));
    free(new_particles);
    free(cumulative);
}

void pf_get_mean(const particle_filter_t *pf, double *x_mean) {
    if (!pf || !x_mean) return;
    memset(x_mean, 0, pf->n_state * sizeof(double));
    for (int i = 0; i < pf->n_particles; i++) {
        for (int j = 0; j < pf->n_state; j++) {
            x_mean[j] += pf->particles[i].weight * pf->particles[i].x[j];
        }
    }
}

void pf_get_covariance(const particle_filter_t *pf, double *cov) {
    if (!pf || !cov) return;
    int n = pf->n_state;
    double mean[KF_STATE_MAX_DIM];
    pf_get_mean(pf, mean);

    memset(cov, 0, n * n * sizeof(double));
    for (int i = 0; i < pf->n_particles; i++) {
        for (int r = 0; r < n; r++) {
            for (int c = 0; c < n; c++) {
                double dr = pf->particles[i].x[r] - mean[r];
                double dc = pf->particles[i].x[c] - mean[c];
                cov[r * n + c] += pf->particles[i].weight * dr * dc;
            }
        }
    }
}

double pf_effective_particles(const particle_filter_t *pf) {
    if (!pf) return 0.0;
    double sum_sq = 0.0;
    for (int i = 0; i < pf->n_particles; i++) {
        sum_sq += pf->particles[i].weight * pf->particles[i].weight;
    }
    if (sum_sq < 1e-30) return pf->n_particles;
    return 1.0 / sum_sq;
}
