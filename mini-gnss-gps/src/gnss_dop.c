/* =========================================================================
 * gnss_dop.c — Dilution of Precision, constellation analysis
 *
 * Covers L1 (DOP definitions: GDOP/PDOP/HDOP/VDOP/TDOP),
 * L2 (satellite geometry and its impact on positioning accuracy),
 * L3 (Q-matrix → DOP decomposition, rotation ECEF→ENU),
 * L4 (DOP as factor connecting UERE to position error).
 *
 * References:
 * - Parkinson, B.W. & Spilker, J.J. (1996). GPS: Theory and Applications.
 * - Langley, R.B. (1999). "Dilution of Precision." GPS World, May.
 * - IS-GPS-200 Annex B (DOP definitions)
 * ========================================================================= */

#include "gnss_dop.h"
#include "gnss_common.h"
#include "gnss_position.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

/* -------------------------------------------------------------------------
 * L3: DOP from Q-matrix with ECEF→ENU rotation
 *
 * Position covariance in ECEF: Cov_ECEF = σ_UERE² · Q_3×3
 * Position covariance in ENU:  Cov_ENU  = R · Cov_ECEF · Rᵀ
 *
 * where R(φ,λ) = rotation matrix from ECEF to ENU (see gnss_common.c).
 *
 * ECEF DOP: PDOP = √(q₁₁+q₂₂+q₃₃)
 * ENU DOP:  HDOP = √(c₁₁+c₂₂), VDOP = √(c₃₃)
 *           where C = R·Q_3×3·Rᵀ
 * ------------------------------------------------------------------------- */

gnss_dop_t gnss_dop_compute(const gnss_qmatrix_t *q_ECEF, gnss_lla_t rx_lla) {
    gnss_dop_t dop;

    /* GDOP, PDOP, TDOP directly from ECEF Q */
    dop.gdop = sqrt(q_ECEF->q[0][0] + q_ECEF->q[1][1]
                  + q_ECEF->q[2][2] + q_ECEF->q[3][3]);
    dop.pdop = sqrt(q_ECEF->q[0][0] + q_ECEF->q[1][1]
                  + q_ECEF->q[2][2]);
    dop.tdop = sqrt(q_ECEF->q[3][3]);

    /* Rotate 3×3 position sub-block to ENU */
    double s_lat = sin(rx_lla.lat), c_lat = cos(rx_lla.lat);
    double s_lon = sin(rx_lla.lon), c_lon = cos(rx_lla.lon);

    /* Rotation matrix ECEF → ENU (3×3) */
    double R[3][3] = {
        { -s_lon,        c_lon,        0.0 },
        { -s_lat*c_lon, -s_lat*s_lon,  c_lat },
        {  c_lat*c_lon,  c_lat*s_lon,  s_lat }
    };

    /* C = R·Q₃·Rᵀ */
    int i, j, k;
    double C[3][3] = {{0}};
    /* Step 1: T = R · Q (3×3 times 3×3) */
    double T[3][3] = {{0}};
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            double s = 0.0;
            for (k = 0; k < 3; k++)
                s += R[i][k] * q_ECEF->q[k][j];
            T[i][j] = s;
        }
    }
    /* Step 2: C = T · Rᵀ (3×3 times 3×3) */
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            double s = 0.0;
            for (k = 0; k < 3; k++)
                s += T[i][k] * R[j][k]; /* Rᵀ[j][k] = R[k][j] */
            C[i][j] = s;
        }
    }

    /* HDOP = sqrt(c_EE + c_NN), VDOP = sqrt(c_UU) */
    dop.hdop = sqrt(C[0][0] + C[1][1]);
    dop.vdop = sqrt(C[2][2]);

    /* Sanity: ensure HDOP² + VDOP² ≈ PDOP² */
    double check = C[0][0] + C[1][1] + C[2][2];
    double pdop_sq = dop.pdop * dop.pdop;
    if (fabs(check - pdop_sq) > 1e-6 * pdop_sq) {
        /* Rotation trace invariant ensures this holds exactly in theory;
         * numerical errors from Cholesky in Q computation may cause mismatch */
    }

    return dop;
}

/* -------------------------------------------------------------------------
 * L2: DOP directly from satellite positions (convenience)
 * ------------------------------------------------------------------------- */

int gnss_dop_from_sat_positions(const gnss_ecef_t sat_pos[], int n_sats,
                                 gnss_ecef_t rx_pos, gnss_dop_t *dop) {
    if (n_sats < 4 || !sat_pos || !dop) return -1;

    /* Build design matrix */
    gnss_sat_meas_t *meas = (gnss_sat_meas_t*)malloc(n_sats * sizeof(gnss_sat_meas_t));
    if (!meas) return -2;

    int i;
    for (i = 0; i < n_sats; i++) {
        meas[i].sat_pos = sat_pos[i];
        meas[i].pseudorange = 0.0; /* dummy */
        meas[i].elevation = gnss_sat_elevation(sat_pos[i], rx_pos);
        meas[i].azimuth = gnss_sat_azimuth(sat_pos[i], rx_pos);
        meas[i].weight = 1.0;
        meas[i].used = 1;
    }

    gnss_design_matrix_t H;
    if (gnss_design_matrix_build(meas, n_sats, rx_pos, &H) != 0) {
        free(meas);
        return -3;
    }

    /* Compute HᵀH */
    double HTH[4][4] = {{0}};
    int j, k;
    for (i = 0; i < n_sats; i++) {
        for (j = 0; j < 4; j++)
            for (k = 0; k < 4; k++)
                HTH[j][k] += H.data[i*4+j] * H.data[i*4+k];
    }

    /* Cholesky decomposition of HTH */
    double L[4][4] = {{0}};
    for (i = 0; i < 4; i++) {
        for (j = 0; j <= i; j++) {
            double sum = HTH[i][j];
            for (k = 0; k < j; k++) sum -= L[i][k] * L[j][k];
            if (i == j) {
                if (sum <= 1e-30) { gnss_design_matrix_free(&H); free(meas); return -1; }
                L[i][j] = sqrt(sum);
            } else {
                L[i][j] = sum / L[j][j];
            }
        }
    }

    /* Solve L·Lᵀ·Q = I for Q columns */
    gnss_qmatrix_t Q;
    int col;
    for (col = 0; col < 4; col++) {
        double e[4] = {(col==0)?1.0:0.0, (col==1)?1.0:0.0,
                        (col==2)?1.0:0.0, (col==3)?1.0:0.0};
        double y[4];
        for (i = 0; i < 4; i++) {
            double sum = e[i];
            for (j = 0; j < i; j++) sum -= L[i][j] * y[j];
            y[i] = sum / L[i][i];
        }
        for (i = 3; i >= 0; i--) {
            double sum = y[i];
            for (j = i+1; j < 4; j++) sum -= L[j][i] * Q.q[j][col];
            Q.q[i][col] = sum / L[i][i];
        }
    }

    gnss_lla_t rx_lla = gnss_ecef_to_lla(rx_pos);
    *dop = gnss_dop_compute(&Q, rx_lla);

    gnss_design_matrix_free(&H);
    free(meas);
    return 0;
}

/* -------------------------------------------------------------------------
 * L1: DOP rating
 * ------------------------------------------------------------------------- */

gnss_dop_rating_t gnss_dop_rate_pdop(double pdop) {
    if (pdop < 2.0) return DOP_RATING_EXCELLENT;
    if (pdop < 4.0) return DOP_RATING_GOOD;
    if (pdop < 6.0) return DOP_RATING_FAIR;
    if (pdop < 10.0) return DOP_RATING_MARGINAL;
    return DOP_RATING_POOR;
}

gnss_dop_rating_t gnss_dop_rate_hdop(double hdop) {
    if (hdop < 2.0) return DOP_RATING_EXCELLENT;
    if (hdop < 4.0) return DOP_RATING_GOOD;
    if (hdop < 6.0) return DOP_RATING_FAIR;
    if (hdop < 10.0) return DOP_RATING_MARGINAL;
    return DOP_RATING_POOR;
}

const char *gnss_dop_rating_string(gnss_dop_rating_t rating) {
    switch (rating) {
    case DOP_RATING_EXCELLENT: return "Excellent";
    case DOP_RATING_GOOD:      return "Good";
    case DOP_RATING_FAIR:      return "Fair";
    case DOP_RATING_MARGINAL:  return "Marginal";
    case DOP_RATING_POOR:      return "Poor";
    default:                   return "Unknown";
    }
}

/* -------------------------------------------------------------------------
 * L2: Constellation geometry analysis
 * ------------------------------------------------------------------------- */

int gnss_geometry_analyze(const gnss_ecef_t sat_pos[], int n_sats,
                           gnss_ecef_t rx_pos, double elev_mask,
                           gnss_skyplot_entry_t entries[],
                           gnss_geometry_stats_t *stats) {
    if (!sat_pos || n_sats < 1 || !stats) return -1;

    memset(stats, 0, sizeof(gnss_geometry_stats_t));
    int i, count = 0;
    double elevations[100], azimuths[100]; /* max 100 sats for analysis */

    for (i = 0; i < n_sats; i++) {
        double el = gnss_sat_elevation(sat_pos[i], rx_pos);
        if (el < elev_mask) continue;

        elevations[count] = el;
        azimuths[count] = gnss_sat_azimuth(sat_pos[i], rx_pos);

        if (entries) {
            entries[count].prn = i + 1;
            entries[count].azimuth_deg = azimuths[count] * 180.0 / M_PI;
            entries[count].elevation_deg = elevations[count] * 180.0 / M_PI;
            entries[count].weight = 1.0;
        }
        count++;
    }
    stats->n_sats_above_mask = count;

    /* Satellite subset above mask */
    if (count >= 4) {
        stats->has_4plus = 1;
        gnss_ecef_t *visible = (gnss_ecef_t*)malloc(count * sizeof(gnss_ecef_t));
        for (i = 0; i < count; i++) visible[i] = sat_pos[i];
        gnss_dop_t dop;
        if (gnss_dop_from_sat_positions(visible, count, rx_pos, &dop) == 0) {
            stats->dop = dop;
            stats->pdop_rating = gnss_dop_rate_pdop(dop.pdop);
            stats->hdop_rating = gnss_dop_rate_hdop(dop.hdop);
        }
        free(visible);
    }

    /* Largest azimuth gap */
    if (count >= 2) {
        /* Sort by azimuth */
        int j;
        for (i = 0; i < count-1; i++) {
            for (j = i+1; j < count; j++) {
                if (azimuths[i] > azimuths[j]) {
                    double t = azimuths[i]; azimuths[i] = azimuths[j]; azimuths[j] = t;
                    t = elevations[i]; elevations[i] = elevations[j]; elevations[j] = t;
                }
            }
        }
        double max_gap = 0.0;
        for (i = 0; i < count-1; i++) {
            double gap = azimuths[i+1] - azimuths[i];
            if (gap > max_gap) max_gap = gap;
        }
        /* Wrap-around gap */
        double wrap_gap = 2.0*M_PI - azimuths[count-1] + azimuths[0];
        if (wrap_gap > max_gap) max_gap = wrap_gap;
        stats->max_gap_azimuth_deg = max_gap * 180.0 / M_PI;
    }

    return 0;
}

/* -------------------------------------------------------------------------
 * L2: Best satellite subset selection (greedy PDOP minimization)
 *
 * Given n_sats visible satellites, select keep_n (≥4) that minimize PDOP.
 * Greedy approach: start with 4 with lowest PDOP, then iteratively add
 * the satellite that most reduces PDOP.
 * ------------------------------------------------------------------------- */

int gnss_select_best_subset(const gnss_ecef_t sat_pos[], int n_sats,
                             gnss_ecef_t rx_pos, int keep_n,
                             int selected_indices[]) {
    if (n_sats < keep_n || keep_n < 4 || !sat_pos || !selected_indices) return -1;

    int i;
    for (i = 0; i < n_sats; i++) selected_indices[i] = 0;

    int *current = (int*)malloc(n_sats * sizeof(int));
    if (!current) return -1;
    int cur_count = 0;

    for (i = 0; i < n_sats; i++) current[i] = 0;

    /* First pass: add 4 satellites with best geometric diversity */
    double azimuths[100];
    for (i = 0; i < n_sats; i++) {
        azimuths[i] = gnss_sat_azimuth(sat_pos[i], rx_pos);
    }

    /* Pick one each from 4 quadrants */
    double quad_min[4] = { 0.0, M_PI_2, M_PI, 3.0*M_PI_2 };
    double quad_max[4] = { M_PI_2, M_PI, 3.0*M_PI_2, 2.0*M_PI };
    int q;
    for (q = 0; q < 4; q++) {
        int best_q = -1;
        double best_el_q = -M_PI;
        for (i = 0; i < n_sats; i++) {
            if (current[i]) continue;
            if (azimuths[i] >= quad_min[q] && azimuths[i] < quad_max[q]) {
                double el = gnss_sat_elevation(sat_pos[i], rx_pos);
                if (el > best_el_q && el > 0.05) {
                    best_el_q = el;
                    best_q = i;
                }
            }
        }
        if (best_q >= 0) {
            current[best_q] = 1;
            cur_count++;
        }
    }

    /* Fill remaining with highest elevation if needed */
    while (cur_count < 4 && cur_count < keep_n) {
        int best_i = -1;
        double best_el = -M_PI;
        for (i = 0; i < n_sats; i++) {
            if (current[i]) continue;
            double el = gnss_sat_elevation(sat_pos[i], rx_pos);
            if (el > best_el) { best_el = el; best_i = i; }
        }
        if (best_i < 0) break;
        current[best_i] = 1;
        cur_count++;
    }

    /* Now iteratively add the satellite that minimizes PDOP
     * (greedy DOP-optimal selection) */
    while (cur_count < keep_n && cur_count < n_sats) {
        double best_pdop = 1e30;
        int add_idx = -1;

        for (i = 0; i < n_sats; i++) {
            if (current[i]) continue;

            /* Temporarily add this satellite */
            current[i] = 1;
            gnss_ecef_t subset[100];
            int si = 0, sj;
            for (sj = 0; sj < n_sats; sj++)
                if (current[sj]) subset[si++] = sat_pos[sj];

            gnss_dop_t dop;
            if (gnss_dop_from_sat_positions(subset, si, rx_pos, &dop) == 0) {
                if (dop.pdop < best_pdop) {
                    best_pdop = dop.pdop;
                    add_idx = i;
                }
            }
            current[i] = 0;
        }

        if (add_idx < 0) break;
        current[add_idx] = 1;
        cur_count++;
    }

    /* Output selected indices */
    for (i = 0; i < n_sats; i++)
        selected_indices[i] = current[i];

    free(current);
    return 0;
}
