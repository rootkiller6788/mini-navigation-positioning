/**
 * @file nav_integrity.c
 * @brief GNSS Integrity Monitoring (RAIM) and Fault Detection
 *
 * L5: RAIM algorithms (residual-based, solution separation)
 * L8: Advanced integrity monitoring
 *
 * Reference: Parkinson & Axelrad (1988), "Autonomous GPS Integrity
 * Monitoring Using the Pseudorange Residual"
 *            RTCA DO-229D, "MOPS for GPS/WAAS Airborne Equipment"
 */

#include "nav_common.h"
#include "nav_gnss.h"
#include "nav_kalman.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

/**
 * @brief RAIM: Residual-based fault detection.
 *
 * Computes test statistic from pseudorange residuals and compares
 * against threshold based on chi-square distribution.
 *
 * @param fault_detected [out] 1 if fault detected
 * @param test_stat [out] test statistic value
 * @param threshold [out] detection threshold
 * @param G [in] geometry matrix (n x 4)
 * @param residuals [in] pseudorange residuals, m
 * @param n_svs [in] number of satellites
 * @param pfa [in] probability of false alarm (typical: 1e-5)
 * @return 0 on success
 */
int nav_raim_residual(int *fault_detected, NAV_PRECISION *test_stat,
                       NAV_PRECISION *threshold,
                       const NAV_PRECISION *G, const NAV_PRECISION *residuals,
                       int n_svs, NAV_PRECISION pfa) {
    (void)pfa;
    if (!fault_detected || !test_stat || !threshold || !G || !residuals || n_svs < 5) {
        if (fault_detected) *fault_detected = 0;
        return -1;
    }
    /* S = I - G*(G^T*G)^{-1}*G^T */
    /* Compute G^T*G */
    NAV_PRECISION GtG[16] = {0};
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            for (int k = 0; k < n_svs; k++)
                GtG[i*4+j] += G[k*4+i] * G[k*4+j];
    /* Invert GtG */
    NAV_PRECISION GtGi[16];
    if (nav_matrix_inverse_spd(GtGi, GtG, 4) != 0) {
        *fault_detected = 0;
        return -1;
    }
    /* SSE = sum of squared weighted residuals */
    NAV_PRECISION sse = 0.0;
    for (int i = 0; i < n_svs; i++)
        sse += residuals[i] * residuals[i];
    /* Test statistic = sqrt(SSE / (n-4)) */
    int dof = n_svs - 4;
    *test_stat = sqrt(sse / dof);
    /* Threshold: approximate from chi-square. For pfa ~ 1e-5:
     * threshold ~ sigma * sqrt(chi2_inv(1-pfa, dof) / dof).
     * Typical: ~3.0-5.0 for aviation-grade RAIM */
    /* Simplified: use constant threshold based on DOF */
    NAV_PRECISION k_chi2[] = {0, 0, 0, 0, 3.5, 3.2, 3.0, 2.9, 2.8, 2.7, 2.6};
    if (dof <= 10 && dof >= 1)
        *threshold = k_chi2[dof];
    else
        *threshold = 2.5;
    *fault_detected = (*test_stat > *threshold) ? 1 : 0;
    return 0;
}

/**
 * @brief RAIM: Fault Detection and Exclusion (FDE).
 *
 * Attempts to identify and exclude a faulty satellite by performing
 * consistency checks on each subset with one SV removed.
 *
 * @param excluded_sv [out] index of excluded satellite (-1 if none)
 * @param G [in] geometry matrix (n x 4)
 * @param residuals [in] residuals, m
 * @param n_svs [in] number of satellites
 * @param pfa [in] false alarm probability
 * @return 0 on success, -1 if insufficient data
 */
int nav_raim_fde(int *excluded_sv,
                  const NAV_PRECISION *G, const NAV_PRECISION *residuals,
                  int n_svs, NAV_PRECISION pfa) {
    (void)pfa;
    if (!excluded_sv || !G || !residuals || n_svs < 6) {
        if (excluded_sv) *excluded_sv = -1;
        return -1;
    }
    *excluded_sv = -1;
    /* For each SV, compute test statistic with that SV removed */
    NAV_PRECISION best_stat = 1e10;
    for (int skip = 0; skip < n_svs; skip++) {
        /* Build reduced G and residuals without SV 'skip' */
        NAV_PRECISION Gr[64], rr[16];
        int nr = 0;
        for (int i = 0; i < n_svs; i++) {
            if (i == skip) continue;
            for (int j = 0; j < 4; j++)
                Gr[nr*4+j] = G[i*4+j];
            rr[nr] = residuals[i];
            nr++;
        }
        int fault;
        NAV_PRECISION stat, thresh;
        if (nav_raim_residual(&fault, &stat, &thresh, Gr, rr, nr, pfa) == 0) {
            if (stat < best_stat) {
                best_stat = stat;
                if (!fault) *excluded_sv = skip;
            }
        }
    }
    /* If best subset has no fault, that skip SV is the faulty one */
    return 0;
}

/**
 * @brief Innovation sequence monitoring for Kalman filter integrity.
 *
 * Monitors the normalized innovation squared (NIS) for filter consistency.
 * NIS = y^T * S^{-1} * y should follow chi-square(m) distribution.
 *
 * @param nis [out] Normalized Innovation Squared
 * @param innovation [in] measurement innovation, length m
 * @param S [in] innovation covariance, m x m
 * @param m [in] measurement dimension
 * @return 0 on success
 */
int nav_nis_test(NAV_PRECISION *nis,
                  const NAV_PRECISION *innovation,
                  const NAV_PRECISION *S, int m) {
    if (!nis || !innovation || !S || m <= 0) return -1;
    /* nis = y^T * S^{-1} * y */
    NAV_PRECISION *Sinv = (NAV_PRECISION*)malloc(m*m*sizeof(NAV_PRECISION));
    if (!Sinv) return -1;
    memcpy(Sinv, S, m*m*sizeof(NAV_PRECISION));
    if (nav_matrix_inverse_spd(Sinv, Sinv, m) != 0) {
        free(Sinv);
        return -1;
    }
    *nis = 0.0;
    for (int i = 0; i < m; i++) {
        NAV_PRECISION sum = 0.0;
        for (int j = 0; j < m; j++)
            sum += Sinv[i*m+j] * innovation[j];
        *nis += innovation[i] * sum;
    }
    free(Sinv);
    return 0;
}

/**
 * @brief Solution Separation RAIM for multi-constellation.
 *
 * Compares position solutions from subsets to detect inconsistencies.
 *
 * @param sep_m [out] maximum horizontal separation, m
 * @param sols [in] array of subset solutions (n_subsets x 4)
 * @param n_subsets [in] number of subsets
 * @return 0 on success
 */
int nav_solution_separation(NAV_PRECISION *sep_m,
                             const NAV_PRECISION *sols, int n_subsets) {
    if (!sep_m || !sols || n_subsets < 2) return -1;
    NAV_PRECISION max_sep = 0.0;
    for (int i = 0; i < n_subsets; i++) {
        for (int j = i+1; j < n_subsets; j++) {
            NAV_PRECISION dx = sols[i*4]-sols[j*4];
            NAV_PRECISION dy = sols[i*4+1]-sols[j*4+1];
            NAV_PRECISION sep = sqrt(dx*dx+dy*dy);
            if (sep > max_sep) max_sep = sep;
        }
    }
    *sep_m = max_sep;
    return 0;
}
