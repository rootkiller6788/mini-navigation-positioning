#ifndef GNSS_DOP_H
#define GNSS_DOP_H
#include "gnss_common.h"

#ifdef __cplusplus
extern "C" {
#endif

/* =========================================================================
 * DOP Theory (L1-L4)
 *
 * The covariance of the position error is:
 *   Cov(Δx) = σ_UERE² · (HᵀH)⁻¹  = σ_UERE² · Q
 *
 * DOP values are derived from the diagonal elements of Q:
 *
 *   GDOP = √(qₓₓ + qᵧᵧ + q_zz + q_tt)    = √(q₁₁+q₂₂+q₃₃+q₄₄)
 *   PDOP = √(qₓₓ + qᵧᵧ + q_zz)            = √(q₁₁+q₂₂+q₃₃)
 *   HDOP = √(qₓₓ + qᵧᵧ)                  (in local ENU frame)
 *   VDOP = √(q_zz)                         (in local ENU frame)
 *   TDOP = √(q_tt)
 *
 * Position accuracy: σ_position = PDOP × σ_UERE
 *
 * DOP Rating (common interpretation):
 *   1-2: Excellent, 2-4: Good, 4-6: Fair, 6-10: Marginal, >10: Poor
 * ========================================================================= */

/** DOP decomposition (ECEF and local-level frames) */
typedef struct {
    double gdop;        /* Geometric DOP */
    double pdop;        /* Position DOP (ECEF) */
    double tdop;        /* Time DOP */
    double hdop;        /* Horizontal DOP (ENU/Local) */
    double vdop;        /* Vertical DOP (ENU/Local) */
} gnss_dop_t;

/** DOP rating category */
typedef enum {
    DOP_RATING_EXCELLENT = 0,   /* DOP < 2 */
    DOP_RATING_GOOD      = 1,   /* 2 ≤ DOP < 4 */
    DOP_RATING_FAIR      = 2,   /* 4 ≤ DOP < 6 */
    DOP_RATING_MARGINAL  = 3,   /* 6 ≤ DOP < 10 */
    DOP_RATING_POOR      = 4    /* DOP ≥ 10 */
} gnss_dop_rating_t;

/* -------------------------------------------------------------------------
 * L1: Q-Matrix (4×4 covariance factor in ECEF)
 * ------------------------------------------------------------------------- */

typedef struct {
    double q[4][4];     /* Q = (HᵀH)⁻¹ in ECEF frame */
} gnss_qmatrix_t;

/* -------------------------------------------------------------------------
 * L2: Satellite geometry analysis
 * ------------------------------------------------------------------------- */

/** Satellite sky-plot entry (used for constellation geometry) */
typedef struct {
    int    prn;
    double azimuth_deg;
    double elevation_deg;
    double weight;
} gnss_skyplot_entry_t;

/** Constellation geometry statistics */
typedef struct {
    int     n_sats_above_mask;  /* Number of sats above elevation mask */
    double  max_dilution;       /* Largest DOP source */
    gnss_dop_t dop;             /* Full DOP decomposition */
    gnss_dop_rating_t pdop_rating;
    gnss_dop_rating_t hdop_rating;
    double  max_gap_azimuth_deg;/* Largest azimuth gap [deg] */
    int     has_4plus;          /* 1 = ≥4 sats (can compute PVT) */
} gnss_geometry_stats_t;

/* -------------------------------------------------------------------------
 * API: DOP computation
 * ------------------------------------------------------------------------- */

/**
 * @brief Compute DOP from Q = (HᵀH)⁻¹ in ECEF frame
 *
 * GDOP = √(q₁₁ + q₂₂ + q₃₃ + q₄₄)
 * PDOP = √(q₁₁ + q₂₂ + q₃₃)
 * TDOP = √(q₄₄)
 *
 * For HDOP/VDOP, the Q matrix is rotated from ECEF to ENU
 * using the receiver LLA position.
 */
gnss_dop_t gnss_dop_compute(const gnss_qmatrix_t *q_ECEF, gnss_lla_t rx_lla);

/**
 * @brief Compute DOP directly from satellite ECEF positions
 *
 * Convenience function: builds H matrix internally and computes DOP.
 *
 * @param sat_pos Array of satellite ECEF positions
 * @param n_sats  Number of satellites (≥4)
 * @param rx_pos  Approximate receiver position (for LOS vectors)
 * @param dop     Output DOP values
 * @return 0 = success, <0 = error (singular geometry)
 */
int gnss_dop_from_sat_positions(const gnss_ecef_t sat_pos[], int n_sats,
                                 gnss_ecef_t rx_pos, gnss_dop_t *dop);

/** @brief Rate DOP value */
gnss_dop_rating_t gnss_dop_rate_pdop(double pdop);
gnss_dop_rating_t gnss_dop_rate_hdop(double hdop);

/** @brief Human-readable DOP rating string */
const char *gnss_dop_rating_string(gnss_dop_rating_t rating);

/**
 * @brief Compute constellation sky-plot and geometry statistics
 *
 * @param sat_pos Array of satellite ECEF positions
 * @param n_sats  Number of satellites
 * @param rx_pos  Receiver ECEF position
 * @param elev_mask Minimum elevation angle [rad] (typ. 5° = 0.087)
 * @param entries  Output sky-plot entries (caller-allocated [n_sats])
 * @param stats    Output geometry statistics
 */
int gnss_geometry_analyze(const gnss_ecef_t sat_pos[], int n_sats,
                           gnss_ecef_t rx_pos, double elev_mask,
                           gnss_skyplot_entry_t entries[],
                           gnss_geometry_stats_t *stats);

/**
 * @brief Select best satellite subset (for RAIM / optimal geometry)
 *
 * Given n_sats > 4, greedily remove the satellite that contributes
 * most to PDOP until only keep_n remain. Used when n_sats >> 4 and
 * computation resources are limited.
 *
 * Complexity: O(n_sats² · matrix-inversion)
 */
int gnss_select_best_subset(const gnss_ecef_t sat_pos[], int n_sats,
                             gnss_ecef_t rx_pos, int keep_n,
                             int selected_indices[]);

#ifdef __cplusplus
}
#endif
#endif /* GNSS_DOP_H */
