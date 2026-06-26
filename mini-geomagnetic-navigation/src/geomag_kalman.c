/**
 * geomag_kalman.c -- Kalman Filtering for Magnetic-Aided Navigation
 *
 * L5: Linear Kalman filter (predict/update)
 * L5: Extended Kalman filter for nonlinear systems
 * L6: INS/Magnetometer integration (15-state EKF)
 *
 * Reference:
 *   Maybeck, "Stochastic Models, Estimation, and Control" (1979)
 *   Brown & Hwang, "Introduction to Random Signals and Applied
 *     Kalman Filtering" (2012)
 *   Groves, "Principles of GNSS, Inertial, and Multisensor
 *     Integrated Navigation Systems" (2013)
 */

#include "geomag_kalman.h"
#include "geomag_math.h"
#include "geomag_model.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ========================================================================
 * L5: Linear Kalman Filter -- internal utility: matrix multiply
 *
 * C = A * B, dimensions: A=mxk, B=kxn, C=mxn (all row-major)
 * ======================================================================== */
static void matrix_mult(const double *A, const double *B, double *C,
                         int m, int k, int n)
{
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++) {
            C[i * n + j] = 0.0;
            for (int p = 0; p < k; p++)
                C[i * n + j] += A[i * k + p] * B[p * n + j];
        }
}

/* Internal: C = A * B^T, A=mxk, B=nxk, C=mxn */
static void matrix_mult_bt(const double *A, const double *B, double *C,
                            int m, int k, int n)
{
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++) {
            C[i * n + j] = 0.0;
            for (int p = 0; p < k; p++)
                C[i * n + j] += A[i * k + p] * B[j * k + p];
        }
}

/* Internal: C = A + B, both mxn */
static void matrix_add(const double *A, const double *B, double *C, int m, int n)
{
    for (int i = 0; i < m * n; i++) C[i] = A[i] + B[i];
}

/* Internal: C = A - B */
static void matrix_sub(const double *A, const double *B, double *C, int m, int n)
{
    for (int i = 0; i < m * n; i++) C[i] = A[i] - B[i];
}

/* Internal: Solve Ax=b for small n via Gauss-Jordan (in-place on A and b) */
static int small_solve(double *A, double *b, int n)
{
    /* Augmented matrix [A|b] */
    double *aug = (double *)malloc(n * (n + 1) * sizeof(double));
    if (!aug) return -1;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) aug[i*(n+1)+j] = A[i*n+j];
        aug[i*(n+1)+n] = b[i];
    }
    for (int col = 0; col < n; col++) {
        int pivot = col;
        double maxv = fabs(aug[col*(n+1)+col]);
        for (int row = col+1; row < n; row++) {
            double av = fabs(aug[row*(n+1)+col]);
            if (av > maxv) { maxv = av; pivot = row; }
        }
        if (maxv < 1e-15) { free(aug); return -1; }
        if (pivot != col)
            for (int j = 0; j <= n; j++) {
                double tmp = aug[col*(n+1)+j];
                aug[col*(n+1)+j] = aug[pivot*(n+1)+j];
                aug[pivot*(n+1)+j] = tmp;
            }
        double piv = aug[col*(n+1)+col];
        for (int j = 0; j <= n; j++) aug[col*(n+1)+j] /= piv;
        for (int row = 0; row < n; row++) {
            if (row != col) {
                double factor = aug[row*(n+1)+col];
                for (int j = 0; j <= n; j++)
                    aug[row*(n+1)+j] -= factor * aug[col*(n+1)+j];
            }
        }
    }
    for (int i = 0; i < n; i++) b[i] = aug[i*(n+1)+n];
    free(aug);
    return 0;
}


/* ========================================================================
 * L5: Linear Kalman Filter -- initialization
 * ======================================================================== */
int kalman_init(KalmanFilter *kf, int n_state, int n_meas)
{
    if (!kf || n_state < 1 || n_meas < 1) return -1;
    memset(kf, 0, sizeof(KalmanFilter));
    kf->n_states = n_state;
    kf->n_meas = n_meas;
    int nn = n_state * n_state;
    int nm = n_state * n_meas;
    int mm = n_meas * n_meas;

    kf->x = (double *)calloc(n_state, sizeof(double));
    kf->P = (double *)calloc(nn, sizeof(double));
    kf->F = (double *)calloc(nn, sizeof(double));
    kf->H = (double *)calloc(nm, sizeof(double));
    kf->Q = (double *)calloc(nn, sizeof(double));
    kf->R = (double *)calloc(mm, sizeof(double));
    kf->K = (double *)calloc(nm, sizeof(double));
    kf->I = (double *)calloc(nn, sizeof(double));
    kf->scratch_nn = (double *)calloc(nn, sizeof(double));
    kf->scratch_nm = (double *)calloc(nm, sizeof(double));
    kf->scratch_mm = (double *)calloc(mm, sizeof(double));

    if (!kf->x || !kf->P || !kf->F || !kf->H || !kf->Q || !kf->R
        || !kf->K || !kf->I || !kf->scratch_nn || !kf->scratch_nm
        || !kf->scratch_mm) {
        kalman_free(kf);
        return -1;
    }

    for (int i = 0; i < n_state; i++) kf->I[i * n_state + i] = 1.0;
    return 0;
}

void kalman_free(KalmanFilter *kf)
{
    if (!kf) return;
    free(kf->x); free(kf->P); free(kf->F); free(kf->H);
    free(kf->Q); free(kf->R); free(kf->K); free(kf->I);
    free(kf->scratch_nn); free(kf->scratch_nm); free(kf->scratch_mm);
    memset(kf, 0, sizeof(KalmanFilter));
}

void kalman_set_initial(KalmanFilter *kf, const double *x0, const double *P0)
{
    if (!kf || !x0 || !P0) return;
    memcpy(kf->x, x0, kf->n_states * sizeof(double));
    memcpy(kf->P, P0, kf->n_states * kf->n_states * sizeof(double));
}

void kalman_set_F(KalmanFilter *kf, const double *F) {
    if (kf && F) memcpy(kf->F, F, kf->n_states*kf->n_states*sizeof(double));
}
void kalman_set_H(KalmanFilter *kf, const double *H) {
    if (kf && H) memcpy(kf->H, H, kf->n_states*kf->n_meas*sizeof(double));
}
void kalman_set_Q(KalmanFilter *kf, const double *Q) {
    if (kf && Q) memcpy(kf->Q, Q, kf->n_states*kf->n_states*sizeof(double));
}
void kalman_set_R(KalmanFilter *kf, const double *R) {
    if (kf && R) memcpy(kf->R, R, kf->n_meas*kf->n_meas*sizeof(double));
}

/* ========================================================================
 * L5: Kalman predict: x = F*x, P = F*P*F^T + Q
 * ======================================================================== */
void kalman_predict(KalmanFilter *kf)
{
    if (!kf) return;
    int n = kf->n_states;

    /* x_pred = F * x */
    matrix_mult(kf->F, kf->x, kf->scratch_nn, n, n, 1);
    memcpy(kf->x, kf->scratch_nn, n * sizeof(double));

    /* P_pred = F * P * F^T + Q */
    /* Step 1: tmp = F * P */
    matrix_mult(kf->F, kf->P, kf->scratch_nn, n, n, n);
    /* Step 2: P_pred = tmp * F^T */
    matrix_mult_bt(kf->scratch_nn, kf->F, kf->P, n, n, n);
    /* Step 3: P_pred += Q */
    matrix_add(kf->P, kf->Q, kf->P, n, n);
}

/* ========================================================================
 * L5: Kalman update: K = P*H^T/(H*P*H^T+R), x = x+K*(z-H*x), P = (I-K*H)*P
 * Uses Joseph form for numerical stability.
 * ======================================================================== */
void kalman_update(KalmanFilter *kf, const double *z)
{
    if (!kf || !z) return;
    int n = kf->n_states, m = kf->n_meas;

    /* Innovation: y = z - H*x */
    matrix_mult(kf->H, kf->x, kf->scratch_mm, m, n, 1);
    for (int i = 0; i < m; i++) kf->scratch_mm[i] = z[i] - kf->scratch_mm[i];

    /* S = H*P*H^T + R */
    /* tmp1 = P * H^T  [n x m] */
    matrix_mult_bt(kf->P, kf->H, kf->scratch_nm, n, n, m);
    /* S = H * tmp1 + R */
    matrix_mult(kf->H, kf->scratch_nm, kf->scratch_mm, m, n, m);
    matrix_add(kf->scratch_mm, kf->R, kf->scratch_mm, m, m);

    /* K = tmp1 * inv(S)  [n x m] */
    double Scopy[25]; /* max 5x5 measurements */
    double Svec[5];
    if (m <= 5) {
        for (int i = 0; i < m*m; i++) Scopy[i] = kf->scratch_mm[i];
        for (int i = 0; i < m; i++) Svec[i] = kf->scratch_nm[i];
        /* Solve S * Krow_i^T = tmp1row_i^T for each row of K */
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) Svec[j] = kf->scratch_nm[i * m + j];
            /* Small linear solve: Scopy * k_row = Svec */
            double *A_copy = (double *)malloc(m * m * sizeof(double));
            if (!A_copy) continue;
            memcpy(A_copy, Scopy, m * m * sizeof(double));
            if (small_solve(A_copy, Svec, m) == 0)
                for (int j = 0; j < m; j++)
                    kf->K[i * m + j] = Svec[j];
            free(A_copy);
        }
    }

    /* State update: x = x + K * y */
    matrix_mult(kf->K, kf->scratch_mm, kf->scratch_nn, n, m, 1);
    for (int i = 0; i < n; i++) kf->x[i] += kf->scratch_nn[i];

    /* Covariance update (Joseph form): P = (I-K*H)*P*(I-K*H)^T + K*R*K^T */
    /* KH = K * H */
    matrix_mult(kf->K, kf->H, kf->scratch_nn, n, m, n);
    /* I_KH = I - KH */
    matrix_sub(kf->I, kf->scratch_nn, kf->scratch_nn, n, n);
    /* tmp = I_KH * P */
    double *tmp = (double *)malloc(n * n * sizeof(double));
    if (!tmp) {
        /* Fallback: P = (I-KH)*P */
        matrix_mult(kf->scratch_nn, kf->P, kf->P, n, n, n);
    } else {
        matrix_mult(kf->scratch_nn, kf->P, tmp, n, n, n);
        /* P = tmp * I_KH^T */
        matrix_mult_bt(tmp, kf->scratch_nn, kf->P, n, n, n);
        /* KRK = K * R * K^T */
        matrix_mult(kf->K, kf->R, kf->scratch_nm, n, m, m);
        matrix_mult_bt(kf->scratch_nm, kf->K, tmp, n, m, n);
        /* P += KRK */
        matrix_add(kf->P, tmp, kf->P, n, n);
        free(tmp);
    }
}


/* ========================================================================
 * L5: Extended Kalman Filter -- initialization
 * ======================================================================== */
int ekf_init(ExtendedKalmanFilter *ekf, int n_state, int n_meas)
{
    if (!ekf || n_state < 1 || n_meas < 1) return -1;
    memset(ekf, 0, sizeof(ExtendedKalmanFilter));
    ekf->n_states = n_state;
    ekf->n_meas = n_meas;
    int nn = n_state * n_state;
    int nm = n_state * n_meas;
    int mm = n_meas * n_meas;

    ekf->x = (double *)calloc(n_state, sizeof(double));
    ekf->P = (double *)calloc(nn, sizeof(double));
    ekf->Q = (double *)calloc(nn, sizeof(double));
    ekf->R = (double *)calloc(mm, sizeof(double));
    ekf->K = (double *)calloc(nm, sizeof(double));
    ekf->I = (double *)calloc(nn, sizeof(double));
    ekf->scratch_nn = (double *)calloc(nn, sizeof(double));
    ekf->scratch_nm = (double *)calloc(nm, sizeof(double));
    ekf->scratch_mm = (double *)calloc(mm, sizeof(double));
    ekf->F_jac = (double *)calloc(nn, sizeof(double));
    ekf->H_jac = (double *)calloc(nm, sizeof(double));

    if (!ekf->x || !ekf->P || !ekf->Q || !ekf->R || !ekf->K || !ekf->I
        || !ekf->scratch_nn || !ekf->scratch_nm || !ekf->scratch_mm
        || !ekf->F_jac || !ekf->H_jac) {
        ekf_free(ekf);
        return -1;
    }
    for (int i = 0; i < n_state; i++) ekf->I[i * n_state + i] = 1.0;
    return 0;
}

void ekf_free(ExtendedKalmanFilter *ekf)
{
    if (!ekf) return;
    free(ekf->x); free(ekf->P); free(ekf->Q); free(ekf->R);
    free(ekf->K); free(ekf->I);
    free(ekf->scratch_nn); free(ekf->scratch_nm); free(ekf->scratch_mm);
    free(ekf->F_jac); free(ekf->H_jac);
    memset(ekf, 0, sizeof(ExtendedKalmanFilter));
}

void ekf_set_initial(ExtendedKalmanFilter *ekf, const double *x0, const double *P0)
{
    if (!ekf || !x0 || !P0) return;
    memcpy(ekf->x, x0, ekf->n_states * sizeof(double));
    memcpy(ekf->P, P0, ekf->n_states * ekf->n_states * sizeof(double));
}

void ekf_set_Q(ExtendedKalmanFilter *ekf, const double *Q) {
    if (ekf && Q) memcpy(ekf->Q, Q, ekf->n_states*ekf->n_states*sizeof(double));
}
void ekf_set_R(ExtendedKalmanFilter *ekf, const double *R) {
    if (ekf && R) memcpy(ekf->R, R, ekf->n_meas*ekf->n_meas*sizeof(double));
}

/* ========================================================================
 * L5: EKF predict: x = f(x, dt), P = F * P * F^T + Q
 * ======================================================================== */
void ekf_predict(ExtendedKalmanFilter *ekf,
                  void (*f)(const double *x, double dt,
                            double *x_pred, double *F_jac, void *user_data),
                  double dt, void *user_data)
{
    if (!ekf || !f) return;
    int n = ekf->n_states;

    /* Compute x_pred and Jacobian F_jac from user function */
    f(ekf->x, dt, ekf->scratch_nn, ekf->F_jac, user_data);
    memcpy(ekf->x, ekf->scratch_nn, n * sizeof(double));

    /* P = F * P * F^T + Q */
    matrix_mult(ekf->F_jac, ekf->P, ekf->scratch_nn, n, n, n);
    matrix_mult_bt(ekf->scratch_nn, ekf->F_jac, ekf->P, n, n, n);
    matrix_add(ekf->P, ekf->Q, ekf->P, n, n);
}

/* ========================================================================
 * L5: EKF update: standard innovation update with user nonlinear h(x)
 * ======================================================================== */
void ekf_update(ExtendedKalmanFilter *ekf,
                 void (*h)(const double *x, double *z_pred, double *H_jac,
                           void *user_data),
                 const double *z, void *user_data)
{
    if (!ekf || !h || !z) return;
    int n = ekf->n_states, m = ekf->n_meas;

    /* Compute predicted measurement and Jacobian */
    h(ekf->x, ekf->scratch_mm, ekf->H_jac, user_data);

    /* Innovation y = z - h(x) */
    for (int i = 0; i < m; i++) ekf->scratch_mm[i] = z[i] - ekf->scratch_mm[i];

    /* S = H*P*H^T + R */
    matrix_mult_bt(ekf->P, ekf->H_jac, ekf->scratch_nm, n, n, m);
    matrix_mult(ekf->H_jac, ekf->scratch_nm, ekf->scratch_mm, m, n, m);
    for (int i = 0; i < m*m; i++) ekf->scratch_mm[i] += ekf->R[i];

    /* K = P*H^T * inv(S) */
    double Scopy[25], Svec[5];
    if (m <= 5) {
        for (int i = 0; i < m*m; i++) Scopy[i] = ekf->scratch_mm[i];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) Svec[j] = ekf->scratch_nm[i*m+j];
            double *Ac = (double *)malloc(m*m*sizeof(double));
            if (!Ac) continue;
            memcpy(Ac, Scopy, m*m*sizeof(double));
            if (small_solve(Ac, Svec, m) == 0)
                for (int j = 0; j < m; j++) ekf->K[i*m+j] = Svec[j];
            free(Ac);
        }
    }

    /* x = x + K * y */
    matrix_mult(ekf->K, ekf->scratch_mm, ekf->scratch_nn, n, m, 1);
    for (int i = 0; i < n; i++) ekf->x[i] += ekf->scratch_nn[i];

    /* P = (I - K*H) * P (Joseph form) */
    matrix_mult(ekf->K, ekf->H_jac, ekf->scratch_nn, n, m, n);
    double *IKH = (double *)malloc(n*n*sizeof(double));
    if (IKH) {
        for (int i = 0; i < n*n; i++) IKH[i] = ekf->I[i] - ekf->scratch_nn[i];
        matrix_mult(IKH, ekf->P, ekf->scratch_nn, n, n, n);
        memcpy(ekf->P, ekf->scratch_nn, n*n*sizeof(double));
        free(IKH);
    }
}

/* ========================================================================
 * L6: Magnetic measurement Jacobian for EKF
 *
 * H = [dF/dlat, dF/dlon, dF/dalt, 0...0]
 *
 * The magnetic field depends only on position, so derivatives w.r.t.
 * velocity, attitude, and bias states are zero.
 *
 * Jacobian is n_meas x n_states, row-major.
 * ======================================================================== */
void mag_measurement_jacobian(const IGRFModel *model,
                               double lat, double lon, double alt,
                               double *H_jac, int n_states)
{
    if (!model || !H_jac || n_states < 3) return;
    int nmeas = 3; /* F, D, I */

    /* Zero out Jacobian */
    memset(H_jac, 0, nmeas * n_states * sizeof(double));

    GeodeticCoord loc0 = { lat, lon, alt };
    double eps_deg = 0.05;
    double dh = 10.0;

    /* For each measurement type (F, D, I), compute derivatives w.r.t position */
    for (int meas = 0; meas < 3; meas++) {
        MagneticElements elem0, elem_plus;
        MagVector B;

        /* Base value */
        if (igrf_compute_field(model, &loc0, &B) != 0) continue;
        compute_magnetic_elements(&B, &elem0);

        double base_val;
        if (meas == 0) base_val = elem0.total_intensity;
        else if (meas == 1) base_val = elem0.declination;
        else base_val = elem0.inclination;

        /* d/dlat */
        GeodeticCoord loc_p = loc0; loc_p.lat += eps_deg;
        if (igrf_compute_field(model, &loc_p, &B) == 0) {
            compute_magnetic_elements(&B, &elem_plus);
            double val_plus = (meas==0) ? elem_plus.total_intensity
                           : (meas==1) ? elem_plus.declination
                           : elem_plus.inclination;
            H_jac[meas * n_states + EKF_IDX_POS_N] = (val_plus - base_val)
                / (eps_deg * DEG2RAD * GEOMAG_EARTH_RADIUS);
        }

        /* d/dlon */
        loc_p = loc0; loc_p.lon += eps_deg;
        if (igrf_compute_field(model, &loc_p, &B) == 0) {
            compute_magnetic_elements(&B, &elem_plus);
            double val_plus = (meas==0) ? elem_plus.total_intensity
                           : (meas==1) ? elem_plus.declination
                           : elem_plus.inclination;
            double cos_lat = cos(lat * DEG2RAD);
            if (fabs(cos_lat) < 1e-10) cos_lat = 1e-10;
            H_jac[meas * n_states + EKF_IDX_POS_E] = (val_plus - base_val)
                / (eps_deg * DEG2RAD * GEOMAG_EARTH_RADIUS * cos_lat);
        }

        /* d/dalt */
        loc_p = loc0; loc_p.alt += dh;
        if (igrf_compute_field(model, &loc_p, &B) == 0) {
            compute_magnetic_elements(&B, &elem_plus);
            double val_plus = (meas==0) ? elem_plus.total_intensity
                           : (meas==1) ? elem_plus.declination
                           : elem_plus.inclination;
            H_jac[meas * n_states + EKF_IDX_POS_D] = (val_plus - base_val) / dh;
        }
    }
}

/* L6: Predict magnetic measurement from state */
void predict_magnetic_measurement(const IGRFModel *model,
                                   double lat, double lon, double alt,
                                   double z_pred[3])
{
    if (!model || !z_pred) return;
    GeodeticCoord loc = { lat, lon, alt };
    MagVector B;
    MagneticElements elem;
    if (igrf_compute_field(model, &loc, &B) == 0) {
        compute_magnetic_elements(&B, &elem);
        z_pred[0] = elem.total_intensity;
        z_pred[1] = elem.declination;
        z_pred[2] = elem.inclination;
    } else {
        z_pred[0] = z_pred[1] = z_pred[2] = 0.0;
    }
}

/* L6: Full INS/MAG EKF update (convenience wrapper) */
void ins_mag_ekf_update(ExtendedKalmanFilter *ekf, const IGRFModel *model,
                         const double z_mag[3],
                         double lat, double lon, double alt)
{
    if (!ekf || !model || !z_mag) return;

    /* Measurement function wrapper */
    /* We use ekf_update which needs a callback. For simplicity, we do a
     * manual update since we have the specific magnetic measurement model. */

    double z_pred[3], H_jac[EKF_N_MEAS_MAG * EKF_N_STATES];
    predict_magnetic_measurement(model, lat, lon, alt, z_pred);
    mag_measurement_jacobian(model, lat, lon, alt, H_jac, ekf->n_states);

    int n = ekf->n_states, m = ekf->n_meas;
    if (m != 3) return;

    /* Innovation */
    double y[3];
    for (int i = 0; i < 3; i++) y[i] = z_mag[i] - z_pred[i];

    /* S = H*P*H^T + R */
    double *tmp_nm = ekf->scratch_nm;
    double *tmp_mm = ekf->scratch_mm;
    matrix_mult_bt(ekf->P, H_jac, tmp_nm, n, n, m);
    matrix_mult(H_jac, tmp_nm, tmp_mm, m, n, m);
    for (int i = 0; i < m*m; i++) tmp_mm[i] += ekf->R[i];

    /* K = P*H^T * inv(S) */
    double Scopy[9];
    for (int i = 0; i < 9; i++) Scopy[i] = tmp_mm[i];
    for (int i = 0; i < n; i++) {
        double Svec[3];
        for (int j = 0; j < 3; j++) Svec[j] = tmp_nm[i*3+j];
        double Ac[9]; memcpy(Ac, Scopy, 9*sizeof(double));
        if (small_solve(Ac, Svec, 3) == 0)
            for (int j = 0; j < 3; j++) ekf->K[i*3+j] = Svec[j];
    }

    /* State update */
    double *dx = ekf->scratch_nn;
    matrix_mult(ekf->K, y, dx, n, m, 1);
    for (int i = 0; i < n; i++) ekf->x[i] += dx[i];

    /* Covariance update: P = (I-KH)*P */
    double *KH = (double *)malloc(n*n*sizeof(double));
    if (KH) {
        matrix_mult(ekf->K, H_jac, KH, n, m, n);
        for (int i = 0; i < n*n; i++) KH[i] = ekf->I[i] - KH[i];
        matrix_mult(KH, ekf->P, ekf->scratch_nn, n, n, n);
        memcpy(ekf->P, ekf->scratch_nn, n*n*sizeof(double));
        free(KH);
    }
}

/* L6: Extract position accuracy from EKF covariance */
void ekf_position_accuracy(const double *P, int n_states,
                            double *sigma_major, double *sigma_minor,
                            double *orientation)
{
    if (!P || !sigma_major || !sigma_minor || !orientation || n_states < 3)
        return;

    /* Extract 2x2 position sub-covariance (North, East) */
    double Pnn = P[EKF_IDX_POS_N * n_states + EKF_IDX_POS_N];
    double Pee = P[EKF_IDX_POS_E * n_states + EKF_IDX_POS_E];
    double Pne = P[EKF_IDX_POS_N * n_states + EKF_IDX_POS_E];

    /* Eigenvalues of 2x2 covariance [Pnn, Pne; Pne, Pee] */
    double tr = Pnn + Pee;
    double det = Pnn * Pee - Pne * Pne;
    double disc = sqrt(tr*tr - 4.0*det);
    double lambda1 = (tr + disc) / 2.0;
    double lambda2 = (tr - disc) / 2.0;

    *sigma_major = sqrt(lambda1);
    *sigma_minor = sqrt(lambda2);

    /* Orientation of major axis from North */
    *orientation = 0.5 * atan2(2.0 * Pne, Pee - Pnn) * RAD2DEG;
}
