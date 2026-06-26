/**
 * @file nav_federated.c
 * @brief Federated Kalman Filter for Multi-Sensor Integrated Navigation
 *
 * Implements Carlson's federated filter with information sharing.
 * L8: Federated (decentralized) Kalman filtering
 *
 * Reference: Carlson (1990), "Federated Square Root Filter for
 * Decentralized Parallel Processes"
 */

#include "nav_kalman.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

/**
 * @brief Federated filter structure: master filter plus N local filters.
 */
typedef struct {
    int   n;              /* common state dimension */
    int   num_local;      /* number of local filters */
    nav_kf_t  *master;    /* master filter */
    nav_kf_t **local;     /* array of local filters */
    NAV_PRECISION *beta;  /* information sharing factors, sum = 1 */
    NAV_PRECISION *x_m;   /* master state */
    NAV_PRECISION *P_m;   /* master covariance */
} nav_federated_t;

/**
 * @brief Initialize federated filter.
 *
 * Master filter maintains the global solution. Each local filter
 * processes a subset of sensors independently.
 *
 * @param ff [out] federated filter
 * @param n common state dimension
 * @param num_local number of local filters
 * @param m_dim array of measurement dimensions for each local filter
 * @return 0 on success
 */
int nav_federated_init(nav_federated_t *ff, int n, int num_local,
                        const int *m_dim) {
    if (!ff || !m_dim || n <= 0 || num_local <= 0) return -1;
    memset(ff, 0, sizeof(nav_federated_t));
    ff->n = n;
    ff->num_local = num_local;
    ff->x_m = (NAV_PRECISION*)calloc(n, sizeof(NAV_PRECISION));
    ff->P_m = (NAV_PRECISION*)calloc(n*n, sizeof(NAV_PRECISION));
    ff->beta = (NAV_PRECISION*)calloc(num_local, sizeof(NAV_PRECISION));
    ff->local = (nav_kf_t**)calloc(num_local, sizeof(nav_kf_t*));
    if (!ff->x_m || !ff->P_m || !ff->beta || !ff->local) return -1;
    for (int i = 0; i < n; i++) ff->P_m[i*n+i] = 1.0;
    /* Default: equal sharing */
    NAV_PRECISION b = 1.0 / num_local;
    for (int i = 0; i < num_local; i++) ff->beta[i] = b;
    /* Allocate local filters */
    for (int i = 0; i < num_local; i++) {
        ff->local[i] = nav_kf_alloc(n, m_dim[i]);
        if (!ff->local[i]) return -1;
    }
    /* Master filter: 0 measurements (prediction only) */
    ff->master = nav_kf_alloc(n, 1);
    if (!ff->master) return -1;
    return 0;
}

/**
 * @brief Perform federated filter time update.
 *
 * Each local filter and the master filter predict independently.
 */
void nav_federated_predict(nav_federated_t *ff) {
    if (!ff) return;
    for (int i = 0; i < ff->num_local; i++)
        if (ff->local[i]) nav_kf_predict(ff->local[i]);
    if (ff->master) nav_kf_predict(ff->master);
}

/**
 * @brief Local filter measurement update.
 *
 * @param ff [in/out] federated filter
 * @param local_idx [in] which local filter to update
 * @param z [in] measurement vector
 * @return 0 on success
 */
int nav_federated_update_local(nav_federated_t *ff, int local_idx,
                                const NAV_PRECISION *z) {
    if (!ff || local_idx < 0 || local_idx >= ff->num_local ||
        !ff->local[local_idx] || !z) return -1;
    return nav_kf_update(ff->local[local_idx], z);
}

/**
 * @brief Federated filter fusion (master update).
 *
 * Fuses all local filter solutions using information-based fusion:
 * P_m^{-1} = sum(beta_i * P_i^{-1})
 * x_m = P_m * sum(beta_i * P_i^{-1} * x_i)
 *
 * Then resets local filters with the fused solution (feedback mode).
 */
int nav_federated_fuse(nav_federated_t *ff) {
    if (!ff) return -1;
    int n = ff->n;
    /* Reset master information */
    NAV_PRECISION *Pm_inv = (NAV_PRECISION*)calloc(n*n, sizeof(NAV_PRECISION));
    NAV_PRECISION *Pmx = (NAV_PRECISION*)calloc(n, sizeof(NAV_PRECISION));
    NAV_PRECISION *Pi_inv = (NAV_PRECISION*)malloc(n*n*sizeof(NAV_PRECISION));
    NAV_PRECISION *tmp = (NAV_PRECISION*)malloc(n*sizeof(NAV_PRECISION));
    if (!Pm_inv || !Pmx || !Pi_inv || !tmp) {
        free(Pm_inv); free(Pmx); free(Pi_inv); free(tmp);
        return -1;
    }
    for (int k = 0; k < ff->num_local; k++) {
        nav_kf_t *kf = ff->local[k];
        if (!kf) continue;
        memcpy(Pi_inv, kf->P, n*n*sizeof(NAV_PRECISION));
        if (nav_matrix_inverse_spd(Pi_inv, Pi_inv, n) != 0) continue;
        /* Pm_inv += beta_k * P_k_inv */
        for (int i = 0; i < n*n; i++)
            Pm_inv[i] += ff->beta[k] * Pi_inv[i];
        /* Pmx += beta_k * P_k_inv * x_k */
        nav_matrix_multiply(tmp, Pi_inv, kf->x, n, n, 1);
        for (int i = 0; i < n; i++)
            Pmx[i] += ff->beta[k] * tmp[i];
    }
    /* P_m = Pm_inv^{-1} */
    NAV_PRECISION *Pm = (NAV_PRECISION*)malloc(n*n*sizeof(NAV_PRECISION));
    if (!Pm) { free(Pm_inv); free(Pmx); free(Pi_inv); free(tmp); return -1; }
    memcpy(Pm, Pm_inv, n*n*sizeof(NAV_PRECISION));
    if (nav_matrix_inverse_spd(ff->P_m, Pm, n) != 0) {
        free(Pm); free(Pm_inv); free(Pmx); free(Pi_inv); free(tmp);
        return -1;
    }
    free(Pm);
    /* x_m = P_m * Pmx */
    nav_matrix_multiply(ff->x_m, ff->P_m, Pmx, n, n, 1);
    /* Feedback: reset all local filters with fused solution */
    for (int k = 0; k < ff->num_local; k++) {
        if (ff->local[k]) {
            memcpy(ff->local[k]->x, ff->x_m, n*sizeof(NAV_PRECISION));
            /* Covariance reset with information sharing factor */
            memcpy(Pi_inv, ff->P_m, n*n*sizeof(NAV_PRECISION));
            if (nav_matrix_inverse_spd(Pi_inv, Pi_inv, n) == 0) {
                for (int i = 0; i < n*n; i++)
                    Pi_inv[i] /= ff->beta[k];
                if (nav_matrix_inverse_spd(ff->local[k]->P, Pi_inv, n) != 0) {
                    memcpy(ff->local[k]->P, ff->P_m, n*n*sizeof(NAV_PRECISION));
                }
            }
        }
    }
    if (ff->master) {
        memcpy(ff->master->x, ff->x_m, n*sizeof(NAV_PRECISION));
        memcpy(ff->master->P, ff->P_m, n*n*sizeof(NAV_PRECISION));
    }
    free(Pm_inv); free(Pmx); free(Pi_inv); free(tmp);
    return 0;
}

/**
 * @brief Free federated filter memory.
 */
void nav_federated_free(nav_federated_t *ff) {
    if (!ff) return;
    if (ff->master) nav_kf_free(ff->master);
    if (ff->local) {
        for (int i = 0; i < ff->num_local; i++)
            if (ff->local[i]) nav_kf_free(ff->local[i]);
        free(ff->local);
    }
    free(ff->x_m); free(ff->P_m); free(ff->beta);
}
