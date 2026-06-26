/**
 * @file nav_kalman.h
 * @brief Kalman Filter Framework for Integrated Navigation
 *
 * L4 Fundamental Laws: Kalman optimal filtering under Gaussian noise.
 * L5 Algorithms: Linear KF, EKF, UKF, Error-State KF, Information Filter.
 *
 * Reference: Kalman (1960), "A New Approach to Linear Filtering
 * and Prediction Problems"
 */

#ifndef NAV_KALMAN_H
#define NAV_KALMAN_H

#include "nav_common.h"
#include <stdlib.h>
#include <string.h>

void nav_matrix_multiply(NAV_PRECISION *C, const NAV_PRECISION *A,
                          const NAV_PRECISION *B, int m, int k, int n);
void nav_matrix_transpose(NAV_PRECISION *B, const NAV_PRECISION *A, int m, int n);
void nav_matrix_add(NAV_PRECISION *C, const NAV_PRECISION *A,
                     const NAV_PRECISION *B, int m, int n);
void nav_matrix_subtract(NAV_PRECISION *C, const NAV_PRECISION *A,
                          const NAV_PRECISION *B, int m, int n);
int nav_cholesky(NAV_PRECISION *A, int n);
void nav_cholesky_solve(NAV_PRECISION *x, const NAV_PRECISION *L,
                         const NAV_PRECISION *b, int n);
int nav_matrix_inverse_spd(NAV_PRECISION *A_inv, const NAV_PRECISION *A, int n);

typedef struct {
    int      n, m;
    NAV_PRECISION *x;
    NAV_PRECISION *P;
    NAV_PRECISION *F;
    NAV_PRECISION *H;
    NAV_PRECISION *Q;
    NAV_PRECISION *R;
    NAV_PRECISION *K;
    NAV_PRECISION *tmp_nn;
    NAV_PRECISION *tmp_nm;
    NAV_PRECISION *tmp_mm;
    NAV_PRECISION *tmp_mn;
    NAV_PRECISION *tmp_n1;
} nav_kf_t;

nav_kf_t *nav_kf_alloc(int n, int m);
void nav_kf_free(nav_kf_t *kf);
void nav_kf_predict(nav_kf_t *kf);
int nav_kf_update(nav_kf_t *kf, const NAV_PRECISION *z);
int nav_kf_update_scalar(nav_kf_t *kf, NAV_PRECISION z, int row);
void nav_kf_set_F(nav_kf_t *kf, const NAV_PRECISION *F);
void nav_kf_set_H(nav_kf_t *kf, const NAV_PRECISION *H);
void nav_kf_set_Q(nav_kf_t *kf, const NAV_PRECISION *Q);
void nav_kf_set_R(nav_kf_t *kf, const NAV_PRECISION *R);
void nav_kf_set_x(nav_kf_t *kf, const NAV_PRECISION *x);
void nav_kf_set_P(nav_kf_t *kf, const NAV_PRECISION *P);
const NAV_PRECISION *nav_kf_get_x(const nav_kf_t *kf);
const NAV_PRECISION *nav_kf_get_P(const nav_kf_t *kf);

typedef struct {
    nav_kf_t   base;
    NAV_PRECISION *x_nominal;
    NAV_PRECISION *u;
    NAV_PRECISION *jac_F;
    NAV_PRECISION *jac_H;
    int         nu;
    int         initialized;
    void (*f)(NAV_PRECISION *x_out, const NAV_PRECISION *x,
              const NAV_PRECISION *u, NAV_PRECISION dt, void *user_data);
    void (*h)(NAV_PRECISION *z_pred, const NAV_PRECISION *x, void *user_data);
    void (*compute_F)(NAV_PRECISION *F, const NAV_PRECISION *x,
                      const NAV_PRECISION *u, NAV_PRECISION dt, void *user_data);
    void (*compute_H)(NAV_PRECISION *H, const NAV_PRECISION *x, void *user_data);
    void *user_data;
} nav_ekf_t;

nav_ekf_t *nav_ekf_alloc(int n, int m, int nu);
void nav_ekf_free(nav_ekf_t *ekf);
void nav_ekf_predict(nav_ekf_t *ekf, NAV_PRECISION dt);
int nav_ekf_update(nav_ekf_t *ekf, const NAV_PRECISION *z);

typedef struct {
    int      n;
    NAV_PRECISION *y;
    NAV_PRECISION *Y;
    NAV_PRECISION *F;
    NAV_PRECISION *Q;
    NAV_PRECISION *tmp_nn;
} nav_info_filter_t;

nav_info_filter_t *nav_info_filter_alloc(int n);
void nav_info_filter_free(nav_info_filter_t *inf);
void nav_info_filter_predict(nav_info_filter_t *inf);
void nav_info_filter_update(nav_info_filter_t *inf,
                             const NAV_PRECISION *H, const NAV_PRECISION *R,
                             const NAV_PRECISION *z, int m);
void nav_info_filter_get_x(NAV_PRECISION *x, const nav_info_filter_t *inf);

void nav_covariance_intersection(NAV_PRECISION *x_fused, NAV_PRECISION *P_fused,
                                  const NAV_PRECISION *x_a,
                                  const NAV_PRECISION *P_a,
                                  const NAV_PRECISION *x_b,
                                  const NAV_PRECISION *P_b, int n);

#endif /* NAV_KALMAN_H */
