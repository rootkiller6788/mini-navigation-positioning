/**
 * @file nav_kalman.c
 * @brief Kalman Filter Framework Implementation
 *
 * Implements linear Kalman filter, EKF infrastructure, information filter,
 * and covariance intersection for multi-sensor fusion.
 *
 * L3: Matrix operations (multiply, transpose, Cholesky)
 * L4: Kalman optimal filter, information filter
 * L8: Covariance intersection for decentralized fusion
 *
 * Reference: Kalman (1960), Maybeck (1979), Julier & Uhlmann (1997)
 */

#include "nav_kalman.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ---------- Matrix Utilities (L3) ------------------------------------- */

void nav_matrix_multiply(NAV_PRECISION *C, const NAV_PRECISION *A,
                          const NAV_PRECISION *B, int m, int k, int n) {
    if (!C || !A || !B) return;
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++) {
            NAV_PRECISION sum = 0.0;
            for (int p = 0; p < k; p++)
                sum += A[i * k + p] * B[p * n + j];
            C[i * n + j] = sum;
        }
}

void nav_matrix_transpose(NAV_PRECISION *B, const NAV_PRECISION *A, int m, int n) {
    if (!B || !A) return;
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            B[j * m + i] = A[i * n + j];
}

void nav_matrix_add(NAV_PRECISION *C, const NAV_PRECISION *A,
                     const NAV_PRECISION *B, int m, int n) {
    if (!C || !A || !B) return;
    for (int i = 0; i < m * n; i++)
        C[i] = A[i] + B[i];
}

void nav_matrix_subtract(NAV_PRECISION *C, const NAV_PRECISION *A,
                          const NAV_PRECISION *B, int m, int n) {
    if (!C || !A || !B) return;
    for (int i = 0; i < m * n; i++)
        C[i] = A[i] - B[i];
}

int nav_cholesky(NAV_PRECISION *A, int n) {
    if (!A || n <= 0) return -1;
    for (int j = 0; j < n; j++) {
        NAV_PRECISION sum = 0.0;
        for (int k = 0; k < j; k++)
            sum += A[j * n + k] * A[j * n + k];
        NAV_PRECISION diag = A[j * n + j] - sum;
        if (diag <= 0.0) return -1;
        A[j * n + j] = sqrt(diag);
        NAV_PRECISION inv_diag = 1.0 / A[j * n + j];
        for (int i = j + 1; i < n; i++) {
            sum = 0.0;
            for (int k = 0; k < j; k++)
                sum += A[i * n + k] * A[j * n + k];
            A[i * n + j] = (A[i * n + j] - sum) * inv_diag;
        }
    }
    return 0;
}

void nav_cholesky_solve(NAV_PRECISION *x, const NAV_PRECISION *L,
                         const NAV_PRECISION *b, int n) {
    if (!x || !L || !b) return;
    NAV_PRECISION *y = (NAV_PRECISION*)malloc(n * sizeof(NAV_PRECISION));
    if (!y) return;
    /* Forward substitution: L*y = b */
    for (int i = 0; i < n; i++) {
        NAV_PRECISION sum = 0.0;
        for (int j = 0; j < i; j++)
            sum += L[i * n + j] * y[j];
        y[i] = (b[i] - sum) / L[i * n + i];
    }
    /* Back substitution: L^T*x = y */
    for (int i = n - 1; i >= 0; i--) {
        NAV_PRECISION sum = 0.0;
        for (int j = i + 1; j < n; j++)
            sum += L[j * n + i] * x[j];
        x[i] = (y[i] - sum) / L[i * n + i];
    }
    free(y);
}

int nav_matrix_inverse_spd(NAV_PRECISION *A_inv, const NAV_PRECISION *A, int n) {
    if (!A_inv || !A || n <= 0) return -1;
    NAV_PRECISION *L = (NAV_PRECISION*)malloc(n * n * sizeof(NAV_PRECISION));
    if (!L) return -1;
    memcpy(L, A, n * n * sizeof(NAV_PRECISION));
    if (nav_cholesky(L, n) != 0) {
        free(L);
        return -1;
    }
    /* Invert by solving for each column of identity */
    for (int j = 0; j < n; j++) {
        NAV_PRECISION ej[64]; /* max dimension for stack allocation */
        NAV_PRECISION *b = (n <= 64) ? ej : (NAV_PRECISION*)malloc(n * sizeof(NAV_PRECISION));
        if (!b) { free(L); return -1; }
        memset(b, 0, n * sizeof(NAV_PRECISION));
        b[j] = 1.0;
        nav_cholesky_solve(&A_inv[j * n], L, b, n);
        if (n > 64) free(b);
    }
    /* Transpose to get correct inverse (since we solved column-wise into rows) */
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++) {
            NAV_PRECISION tmp = A_inv[i * n + j];
            A_inv[i * n + j] = A_inv[j * n + i];
            A_inv[j * n + i] = tmp;
        }
    free(L);
    return 0;
}

/* ---------- Linear Kalman Filter (L4) --------------------------------- */

nav_kf_t *nav_kf_alloc(int n, int m) {
    nav_kf_t *kf = (nav_kf_t*)calloc(1, sizeof(nav_kf_t));
    if (!kf) return NULL;
    kf->n = n; kf->m = m;
    kf->x      = (NAV_PRECISION*)calloc(n, sizeof(NAV_PRECISION));
    kf->P      = (NAV_PRECISION*)calloc(n*n, sizeof(NAV_PRECISION));
    kf->F      = (NAV_PRECISION*)calloc(n*n, sizeof(NAV_PRECISION));
    kf->H      = (NAV_PRECISION*)calloc(m*n, sizeof(NAV_PRECISION));
    kf->Q      = (NAV_PRECISION*)calloc(n*n, sizeof(NAV_PRECISION));
    kf->R      = (NAV_PRECISION*)calloc(m*m, sizeof(NAV_PRECISION));
    kf->K      = (NAV_PRECISION*)calloc(n*m, sizeof(NAV_PRECISION));
    kf->tmp_nn = (NAV_PRECISION*)calloc(n*n, sizeof(NAV_PRECISION));
    kf->tmp_nm = (NAV_PRECISION*)calloc(n*m, sizeof(NAV_PRECISION));
    kf->tmp_mm = (NAV_PRECISION*)calloc(m*m, sizeof(NAV_PRECISION));
    kf->tmp_mn = (NAV_PRECISION*)calloc(m*n, sizeof(NAV_PRECISION));
    kf->tmp_n1 = (NAV_PRECISION*)calloc(n, sizeof(NAV_PRECISION));
    if (!kf->x || !kf->P || !kf->F || !kf->H || !kf->Q || !kf->R ||
        !kf->K || !kf->tmp_nn || !kf->tmp_mm || !kf->tmp_n1) {
        nav_kf_free(kf);
        return NULL;
    }
    /* Initialize F and P as identity */
    for (int i = 0; i < n; i++) {
        kf->F[i*n+i] = 1.0;
        kf->P[i*n+i] = 1.0;
    }
    /* Initialize R as identity */
    for (int i = 0; i < m; i++)
        kf->R[i*m+i] = 1.0;
    return kf;
}

void nav_kf_free(nav_kf_t *kf) {
    if (!kf) return;
    free(kf->x); free(kf->P); free(kf->F); free(kf->H);
    free(kf->Q); free(kf->R); free(kf->K);
    free(kf->tmp_nn); free(kf->tmp_nm); free(kf->tmp_mm);
    free(kf->tmp_mn); free(kf->tmp_n1);
    free(kf);
}

void nav_kf_predict(nav_kf_t *kf) {
    /* x = F*x */
    /* P = F*P*F^T + Q */
    if (!kf) return;
    int n = kf->n;
    /* tmp_n1 = F * x */
    nav_matrix_multiply(kf->tmp_n1, kf->F, kf->x, n, n, 1);
    memcpy(kf->x, kf->tmp_n1, n * sizeof(NAV_PRECISION));
    /* tmp_nn = F * P */
    nav_matrix_multiply(kf->tmp_nn, kf->F, kf->P, n, n, n);
    /* tmp_mn = (F*P)^T = P^T*F^T */
    nav_matrix_transpose(kf->tmp_mn, kf->tmp_nn, n, n);
    /* P = (F*P) * F^T = tmp_nn * tmp_mn^T... actually need F*P*F^T */
    /* tmp_nm = tmp_nn * F^T */
    nav_matrix_multiply(kf->P, kf->tmp_nn, kf->tmp_mn, n, n, n);
    /* Wait, we need: P_new = F*P*F^T + Q */
    /* tmp_nn = F*P, we need F*P*F^T. Let F_P_Ft = tmp_nn * F^T */
    /* F^T stored in tmp_mn (n x n) */
    /* Use tmp_nm for result: F*P*F^T */
    /* redo cleaner: */
    /* Step 1: tmp_nn = F * P */
    nav_matrix_multiply(kf->tmp_nn, kf->F, kf->P, n, n, n);
    /* Step 2: F^T in tmp_mn */
    nav_matrix_transpose(kf->tmp_mn, kf->F, n, n);
    /* Step 3: P = tmp_nn * tmp_mn */
    nav_matrix_multiply(kf->P, kf->tmp_nn, kf->tmp_mn, n, n, n);
    /* Step 4: P += Q */
    nav_matrix_add(kf->P, kf->P, kf->Q, n, n);
}

int nav_kf_update(nav_kf_t *kf, const NAV_PRECISION *z) {
    /* y = z - H*x                 (innovation)
     * S = H*P*H^T + R            (innovation covariance)
     * K = P*H^T * inv(S)         (Kalman gain)
     * x = x + K*y                (state update)
     * P = (I-K*H)*P*(I-K*H)^T + K*R*K^T  (Joseph form) */
    if (!kf || !z) return -1;
    int n = kf->n, m = kf->m;
    /* H^T in tmp_mn */
    nav_matrix_transpose(kf->tmp_mn, kf->H, m, n);
    /* y = z - H*x */
    nav_matrix_multiply(kf->tmp_n1, kf->H, kf->x, m, n, 1);
    for (int i = 0; i < m; i++)
        kf->tmp_n1[i] = z[i] - kf->tmp_n1[i];
    /* S = H*P*H^T + R */
    /* tmp_nm = H * P */
    nav_matrix_multiply(kf->tmp_nm, kf->H, kf->P, m, n, n);
    /* S = tmp_nm * H^T */
    nav_matrix_multiply(kf->tmp_mm, kf->tmp_nm, kf->tmp_mn, m, n, m);
    /* S += R */
    nav_matrix_add(kf->tmp_mm, kf->tmp_mm, kf->R, m, m);
    /* K = P * H^T * inv(S) */
    /* tmp_nm = P * H^T */
    nav_matrix_multiply(kf->tmp_nm, kf->P, kf->tmp_mn, n, n, m);
    /* inv(S) -> tmp_mm (overwriting S with its inverse) */
    if (nav_matrix_inverse_spd(kf->tmp_mm, kf->tmp_mm, m) != 0)
        return -1;
    /* K = tmp_nm * inv(S) */
    nav_matrix_multiply(kf->K, kf->tmp_nm, kf->tmp_mm, n, m, m);
    /* x = x + K*y */
    nav_matrix_multiply(kf->tmp_n1 + m, kf->K, kf->tmp_n1, n, m, 1);
    for (int i = 0; i < n; i++)
        kf->x[i] += kf->tmp_n1[m + i];
    /* Joseph form: P = (I-K*H)*P*(I-K*H)^T + K*R*K^T */
    /* I_KH = I - K*H */
    /* tmp_nn = K * H */
    nav_matrix_multiply(kf->tmp_nn, kf->K, kf->H, n, m, n);
    /* tmp_nn = I - K*H */
    for (int i = 0; i < n*n; i++)
        kf->tmp_nn[i] = -kf->tmp_nn[i];
    for (int i = 0; i < n; i++)
        kf->tmp_nn[i*n+i] += 1.0;
    /* tmp_mn = (I-K*H) * P */
    nav_matrix_multiply(kf->tmp_mn, kf->tmp_nn, kf->P, n, n, n);
    /* P = tmp_mn * (I-K*H)^T */
    /* (I-K*H)^T in tmp_nm */
    nav_matrix_transpose(kf->tmp_nm, kf->tmp_nn, n, n);
    nav_matrix_multiply(kf->P, kf->tmp_mn, kf->tmp_nm, n, n, n);
    /* K*R*K^T */
    /* tmp_nm = K * R */
    nav_matrix_multiply(kf->tmp_nm, kf->K, kf->R, n, m, m);
    /* K^T in tmp_mn */
    nav_matrix_transpose(kf->tmp_mn, kf->K, n, m);
    /* P += (K*R) * K^T */
    nav_matrix_multiply(kf->tmp_nn, kf->tmp_nm, kf->tmp_mn, n, m, n);
    nav_matrix_add(kf->P, kf->P, kf->tmp_nn, n, n);
    return 0;
}

int nav_kf_update_scalar(nav_kf_t *kf, NAV_PRECISION z, int row) {
    /* Scalar measurement update: more efficient for diagonal R */
    if (!kf || row < 0 || row >= kf->m) return -1;
    int n = kf->n;
    /* H_row is row-th row of H */
    NAV_PRECISION *h = (NAV_PRECISION*)malloc(n * sizeof(NAV_PRECISION));
    NAV_PRECISION *Ph = (NAV_PRECISION*)malloc(n * sizeof(NAV_PRECISION));
    if (!h || !Ph) { free(h); free(Ph); return -1; }
    for (int i = 0; i < n; i++)
        h[i] = kf->H[row * n + i];
    /* y = z - h*x */
    NAV_PRECISION y = z;
    for (int i = 0; i < n; i++)
        y -= h[i] * kf->x[i];
    /* S = h*P*h^T + R(row,row) */
    NAV_PRECISION S = kf->R[row * kf->m + row];
    for (int i = 0; i < n; i++) {
        NAV_PRECISION sum = 0.0;
        for (int j = 0; j < n; j++)
            sum += kf->P[i*n+j] * h[j];
        Ph[i] = sum;
        S += h[i] * sum;
    }
    if (fabs(S) < 1e-15) { free(h); free(Ph); return -1; }
    NAV_PRECISION invS = 1.0 / S;
    /* K = Ph / S */
    /* x = x + K*y */
    for (int i = 0; i < n; i++)
        kf->x[i] += Ph[i] * invS * y;
    /* P = P - K*S*K^T = P - Ph*Ph^T / S */
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            kf->P[i*n+j] -= Ph[i] * Ph[j] * invS;
    free(h); free(Ph);
    return 0;
}

void nav_kf_set_F(nav_kf_t *kf, const NAV_PRECISION *F) {
    if (!kf || !F) return;
    memcpy(kf->F, F, kf->n * kf->n * sizeof(NAV_PRECISION));
}

void nav_kf_set_H(nav_kf_t *kf, const NAV_PRECISION *H) {
    if (!kf || !H) return;
    memcpy(kf->H, H, kf->m * kf->n * sizeof(NAV_PRECISION));
}

void nav_kf_set_Q(nav_kf_t *kf, const NAV_PRECISION *Q) {
    if (!kf || !Q) return;
    memcpy(kf->Q, Q, kf->n * kf->n * sizeof(NAV_PRECISION));
}

void nav_kf_set_R(nav_kf_t *kf, const NAV_PRECISION *R) {
    if (!kf || !R) return;
    memcpy(kf->R, R, kf->m * kf->m * sizeof(NAV_PRECISION));
}

void nav_kf_set_x(nav_kf_t *kf, const NAV_PRECISION *x) {
    if (!kf || !x) return;
    memcpy(kf->x, x, kf->n * sizeof(NAV_PRECISION));
}

void nav_kf_set_P(nav_kf_t *kf, const NAV_PRECISION *P) {
    if (!kf || !P) return;
    memcpy(kf->P, P, kf->n * kf->n * sizeof(NAV_PRECISION));
}

const NAV_PRECISION *nav_kf_get_x(const nav_kf_t *kf) {
    return kf ? kf->x : NULL;
}

const NAV_PRECISION *nav_kf_get_P(const nav_kf_t *kf) {
    return kf ? kf->P : NULL;
}

nav_ekf_t *nav_ekf_alloc(int n, int m, int nu) {
    nav_ekf_t *ekf = (nav_ekf_t*)calloc(1, sizeof(nav_ekf_t));
    if (!ekf) return NULL;
    ekf->nu = nu;
    ekf->x_nominal = (NAV_PRECISION*)calloc(n, sizeof(NAV_PRECISION));
    ekf->u = (nu > 0) ? (NAV_PRECISION*)calloc(nu, sizeof(NAV_PRECISION)) : NULL;
    ekf->jac_F = (NAV_PRECISION*)calloc(n*n, sizeof(NAV_PRECISION));
    ekf->jac_H = (NAV_PRECISION*)calloc(m*n, sizeof(NAV_PRECISION));
    if (!ekf->x_nominal || (nu>0 && !ekf->u) || !ekf->jac_F || !ekf->jac_H) {
        nav_ekf_free(ekf); return NULL;
    }
    for (int i = 0; i < n; i++) {
        ekf->jac_F[i*n+i] = 1.0;
        if (i < m) ekf->jac_H[i*n+i] = 1.0;
    }
    return ekf;
}

void nav_ekf_free(nav_ekf_t *ekf) {
    if (!ekf) return;
    free(ekf->x_nominal); free(ekf->u);
    free(ekf->jac_F); free(ekf->jac_H);
    free(ekf);
}

void nav_ekf_predict(nav_ekf_t *ekf, NAV_PRECISION dt) {
    if (!ekf || !ekf->f || !ekf->compute_F) return;
    NAV_PRECISION *x_pred = (NAV_PRECISION*)malloc(ekf->base.n * sizeof(NAV_PRECISION));
    if (!x_pred) return;
    ekf->f(x_pred, ekf->base.x, ekf->u, dt, ekf->user_data);
    memcpy(ekf->base.x, x_pred, ekf->base.n * sizeof(NAV_PRECISION));
    ekf->compute_F(ekf->jac_F, ekf->base.x, ekf->u, dt, ekf->user_data);
    memcpy(ekf->base.F, ekf->jac_F, ekf->base.n * ekf->base.n * sizeof(NAV_PRECISION));
    nav_kf_predict(&ekf->base);
    free(x_pred);
}

int nav_ekf_update(nav_ekf_t *ekf, const NAV_PRECISION *z) {
    if (!ekf || !z || !ekf->h || !ekf->compute_H) return -1;
    ekf->compute_H(ekf->jac_H, ekf->base.x, ekf->user_data);
    memcpy(ekf->base.H, ekf->jac_H, ekf->base.m * ekf->base.n * sizeof(NAV_PRECISION));
    return nav_kf_update(&ekf->base, z);
}

nav_info_filter_t *nav_info_filter_alloc(int n) {
    nav_info_filter_t *inf = (nav_info_filter_t*)calloc(1, sizeof(nav_info_filter_t));
    if (!inf) return NULL;
    inf->n = n;
    inf->y = (NAV_PRECISION*)calloc(n, sizeof(NAV_PRECISION));
    inf->Y = (NAV_PRECISION*)calloc(n*n, sizeof(NAV_PRECISION));
    inf->F = (NAV_PRECISION*)calloc(n*n, sizeof(NAV_PRECISION));
    inf->Q = (NAV_PRECISION*)calloc(n*n, sizeof(NAV_PRECISION));
    inf->tmp_nn = (NAV_PRECISION*)calloc(n*n, sizeof(NAV_PRECISION));
    if (!inf->y || !inf->Y || !inf->F || !inf->Q || !inf->tmp_nn) {
        nav_info_filter_free(inf); return NULL;
    }
    for (int i = 0; i < n; i++) { inf->F[i*n+i] = 1.0; inf->Y[i*n+i] = 1.0; }
    return inf;
}

void nav_info_filter_free(nav_info_filter_t *inf) {
    if (!inf) return;
    free(inf->y); free(inf->Y); free(inf->F);
    free(inf->Q); free(inf->tmp_nn);
    free(inf);
}

void nav_info_filter_predict(nav_info_filter_t *inf) {
    if (!inf) return;
    int n = inf->n;
    /* P = Y^{-1}, P_pred = P + Q, Y = P_pred^{-1} (for F=I) */
    NAV_PRECISION *P = inf->tmp_nn;
    memcpy(P, inf->Y, n*n*sizeof(NAV_PRECISION));
    if (nav_matrix_inverse_spd(P, P, n) == 0) {
        nav_matrix_add(P, P, inf->Q, n, n);
        nav_matrix_inverse_spd(inf->Y, P, n);
    }
}

void nav_info_filter_update(nav_info_filter_t *inf,
                             const NAV_PRECISION *H, const NAV_PRECISION *R,
                             const NAV_PRECISION *z, int m) {
    if (!inf || !H || !R || !z || m <= 0) return;
    int n = inf->n;
    NAV_PRECISION *R_inv = (NAV_PRECISION*)malloc(m*m*sizeof(NAV_PRECISION));
    NAV_PRECISION *Ht = (NAV_PRECISION*)malloc(n*m*sizeof(NAV_PRECISION));
    NAV_PRECISION *HtR = (NAV_PRECISION*)malloc(n*m*sizeof(NAV_PRECISION));
    NAV_PRECISION *HtRH = (NAV_PRECISION*)malloc(n*n*sizeof(NAV_PRECISION));
    if (!R_inv || !Ht || !HtR || !HtRH) {
        free(R_inv); free(Ht); free(HtR); free(HtRH); return;
    }
    memcpy(R_inv, R, m*m*sizeof(NAV_PRECISION));
    if (nav_matrix_inverse_spd(R_inv, R_inv, m) != 0) {
        free(R_inv); free(Ht); free(HtR); free(HtRH); return;
    }
    nav_matrix_transpose(Ht, H, m, n);
    nav_matrix_multiply(HtR, Ht, R_inv, n, m, m);
    nav_matrix_multiply(HtRH, HtR, H, n, m, n);
    nav_matrix_add(inf->Y, inf->Y, HtRH, n, n);
    NAV_PRECISION *Rz = (NAV_PRECISION*)malloc(m*sizeof(NAV_PRECISION));
    if (Rz) {
        nav_matrix_multiply(Rz, R_inv, z, m, m, 1);
        for (int i = 0; i < n; i++) {
            NAV_PRECISION sum = 0.0;
            for (int j = 0; j < m; j++) sum += Ht[i*m+j] * Rz[j];
            inf->y[i] += sum;
        }
        free(Rz);
    }
    free(R_inv); free(Ht); free(HtR); free(HtRH);
}

void nav_info_filter_get_x(NAV_PRECISION *x, const nav_info_filter_t *inf) {
    if (!x || !inf) return;
    NAV_PRECISION *Yi = (NAV_PRECISION*)malloc(inf->n*inf->n*sizeof(NAV_PRECISION));
    if (!Yi) return;
    memcpy(Yi, inf->Y, inf->n*inf->n*sizeof(NAV_PRECISION));
    if (nav_matrix_inverse_spd(Yi, Yi, inf->n) == 0)
        nav_matrix_multiply(x, Yi, inf->y, inf->n, inf->n, 1);
    free(Yi);
}

void nav_covariance_intersection(NAV_PRECISION *x_fused, NAV_PRECISION *P_fused,
                                  const NAV_PRECISION *x_a, const NAV_PRECISION *P_a,
                                  const NAV_PRECISION *x_b, const NAV_PRECISION *P_b, int n) {
    if (!x_fused || !P_fused || !x_a || !P_a || !x_b || !P_b || n <= 0) return;
    NAV_PRECISION *ia = (NAV_PRECISION*)malloc(n*n*sizeof(NAV_PRECISION));
    NAV_PRECISION *ib = (NAV_PRECISION*)malloc(n*n*sizeof(NAV_PRECISION));
    NAV_PRECISION *Pci = (NAV_PRECISION*)malloc(n*n*sizeof(NAV_PRECISION));
    if (!ia || !ib || !Pci) { free(ia); free(ib); free(Pci); return; }
    memcpy(ia, P_a, n*n*sizeof(NAV_PRECISION));
    memcpy(ib, P_b, n*n*sizeof(NAV_PRECISION));
    if (nav_matrix_inverse_spd(ia, ia, n) || nav_matrix_inverse_spd(ib, ib, n)) {
        free(ia); free(ib); free(Pci); return;
    }
    NAV_PRECISION best = 0.5, best_det = 1e300;
    for (int k=0; k<20; k++) {
        NAV_PRECISION w = 0.01 + k*0.049;
        for (int i=0; i<n*n; i++) Pci[i] = w*ia[i] + (1.0-w)*ib[i];
        NAV_PRECISION *L = (NAV_PRECISION*)malloc(n*n*sizeof(NAV_PRECISION));
        if (L) {
            memcpy(L, Pci, n*n*sizeof(NAV_PRECISION));
            if (nav_cholesky(L, n) == 0) {
                NAV_PRECISION det = 1.0;
                for (int i=0; i<n; i++) det *= L[i*n+i]*L[i*n+i];
                if (det < best_det) { best_det = det; best = w; }
            }
            free(L);
        }
    }
    for (int i=0; i<n*n; i++) Pci[i] = best*ia[i] + (1.0-best)*ib[i];
    NAV_PRECISION *invPci = (NAV_PRECISION*)malloc(n*n*sizeof(NAV_PRECISION));
    if (invPci) {
        memcpy(invPci, Pci, n*n*sizeof(NAV_PRECISION));
        if (nav_matrix_inverse_spd(P_fused, invPci, n) == 0) {
            NAV_PRECISION *t1=calloc(n,sizeof(NAV_PRECISION));
            NAV_PRECISION *t2=calloc(n,sizeof(NAV_PRECISION));
            if (t1&&t2) {
                nav_matrix_multiply(t1, ia, x_a, n, n, 1);
                nav_matrix_multiply(t2, ib, x_b, n, n, 1);
                for(int i=0;i<n;i++) t1[i]=best*t1[i]+(1.0-best)*t2[i];
                nav_matrix_multiply(x_fused, P_fused, t1, n, n, 1);
                free(t1); free(t2);
            }
        }
        free(invPci);
    }
    free(ia); free(ib); free(Pci);
}
