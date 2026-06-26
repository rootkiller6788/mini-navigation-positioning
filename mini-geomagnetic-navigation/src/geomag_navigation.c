/**
 * geomag_navigation.c -- Magnetic Navigation Algorithms
 *
 * L5: Magnetic contour matching (MAGCOM), correlation-based
 * L6: IGRF single-point positioning, compass dead-reckoning
 * L6: Magnetic gradient-based positioning updates
 * L7: GPS-denied navigation applications
 *
 * Reference:
 *   Goldenberg, "Geomagnetic Navigation beyond the Magnetic Compass", IEEE PLANS 2006.
 *   Canciani & Raquet, "Absolute Positioning Using the Earth Magnetic Field",
 *     NAVIGATION, Vol.63, No.2, 2016, pp. 111-126.
 *   Shockley & Raquet, "Navigation Using Magnetic Field Gradients",
 *     IEEE Trans. AES, Vol.50, No.1, 2014.
 */

#include "geomag_navigation.h"
#include "geomag_math.h"
#include "geomag_model.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ========================================================================
 * L6: MAGCOM via Normalized Cross-Correlation
 *
 * Slides measurement profile across magnetic map and computes NCC.
 * Best position = position of maximum correlation.
 *
 * NCC = sum[(z_i-z_bar)*(m_i-m_bar)] / sqrt(sum[(z_i-z_bar)^2]*sum[(m_i-m_bar)^2])
 *
 * Complexity: O(N_search * N_points) with bilinear interpolation per lookup.
 * ======================================================================== */
int magcom_correlation_match(const MagneticMap *map,
                              const double *measurements,
                              const double *trajectory_lat,
                              const double *trajectory_lon,
                              int N,
                              double search_radius,
                              double search_step,
                              double *est_lat, double *est_lon,
                              double *confidence)
{
    if (!map || !measurements || !est_lat || !est_lon || !confidence || N < 3)
        return -1;

    double z_bar = 0.0;
    for (int i = 0; i < N; i++) z_bar += measurements[i];
    z_bar /= N;

    double z_ssq = 0.0;
    for (int i = 0; i < N; i++) {
        double dz = measurements[i] - z_bar;
        z_ssq += dz * dz;
    }
    if (z_ssq < 1e-10) z_ssq = 1.0;

    double best_ncc = -2.0;
    double best_lat = 0.0, best_lon = 0.0;
    int nsteps = (int)(search_radius / fabs(search_step)) + 1;

    for (int i = -nsteps; i <= nsteps; i++) {
        double off_lat = i * search_step;
        for (int j = -nsteps; j <= nsteps; j++) {
            double off_lon = j * search_step;
            if (off_lat*off_lat + off_lon*off_lon > search_radius*search_radius)
                continue;

            double m_bar = 0.0;
            double *m_vals = (double *)malloc(N * sizeof(double));
            if (!m_vals) continue;
            int valid = 1;

            for (int k = 0; k < N && valid; k++) {
                GeodeticCoord pt = { off_lat+trajectory_lat[k],
                                     off_lon+trajectory_lon[k], 0.0 };
                double map_val;
                if (magmap_bilinear_interpolate(map, &pt, &map_val) != 0)
                    valid = 0;
                else { m_vals[k] = map_val; m_bar += map_val; }
            }
            if (!valid) { free(m_vals); continue; }
            m_bar /= N;

            double m_ssq = 0.0, cross = 0.0;
            for (int k = 0; k < N; k++) {
                double dz = measurements[k] - z_bar;
                double dm = m_vals[k] - m_bar;
                cross += dz * dm;
                m_ssq += dm * dm;
            }
            if (m_ssq < 1e-10) { free(m_vals); continue; }

            double ncc = cross / sqrt(z_ssq * m_ssq);
            if (ncc > best_ncc) {
                best_ncc = ncc;
                best_lat = off_lat;
                best_lon = off_lon;
            }
            free(m_vals);
        }
    }
    *est_lat = best_lat;
    *est_lon = best_lon;
    *confidence = (best_ncc > 0.0) ? best_ncc : 0.0;
    return 0;
}

/* ========================================================================
 * L5: MAGCOM via Mean Absolute Difference (robust to outliers)
 * ======================================================================== */
int magcom_mad_match(const MagneticMap *map,
                      const double *measurements,
                      const double *trajectory_lat,
                      const double *trajectory_lon,
                      int N,
                      double search_radius, double search_step,
                      double *est_lat, double *est_lon,
                      double *min_mad)
{
    if (!map || !measurements || !est_lat || !est_lon || !min_mad || N < 3)
        return -1;

    double best_mad = 1e100;
    double best_lat = 0.0, best_lon = 0.0;
    int nsteps = (int)(search_radius / fabs(search_step)) + 1;

    for (int i = -nsteps; i <= nsteps; i++) {
        double off_lat = i * search_step;
        for (int j = -nsteps; j <= nsteps; j++) {
            double off_lon = j * search_step;
            if (off_lat*off_lat + off_lon*off_lon > search_radius*search_radius)
                continue;

            double sum_abs = 0.0;
            int valid = 1;
            for (int k = 0; k < N && valid; k++) {
                GeodeticCoord pt = { off_lat+trajectory_lat[k],
                                     off_lon+trajectory_lon[k], 0.0 };
                double map_val;
                if (magmap_bilinear_interpolate(map, &pt, &map_val) != 0)
                    valid = 0;
                else sum_abs += fabs(measurements[k] - map_val);
            }
            if (!valid) continue;
            double mad = sum_abs / N;
            if (mad < best_mad) { best_mad = mad; best_lat = off_lat; best_lon = off_lon; }
        }
    }
    *est_lat = best_lat; *est_lon = best_lon; *min_mad = best_mad;
    return 0;
}

/* ========================================================================
 * L5: Generate expected magnetic profile along trajectory
 * Complexity: O(N) bilinear lookups.
 * ======================================================================== */
int generate_magnetic_profile(const MagneticMap *map,
                               const GeodeticCoord *start,
                               double heading, double step_size,
                               int N, double *profile,
                               double *profile_lats, double *profile_lons)
{
    if (!map || !start || !profile || N < 1) return -1;
    GeodeticCoord cur = *start;
    for (int i = 0; i < N; i++) {
        if (i > 0) great_circle_destination(start, heading, i*step_size, &cur);
        if (profile_lats) profile_lats[i] = cur.lat;
        if (profile_lons) profile_lons[i] = cur.lon;
        if (magmap_bilinear_interpolate(map, &cur, &profile[i]) != 0)
            profile[i] = 0.0;
    }
    return 0;
}


/* ========================================================================
 * L6: IGRF-based single-point positioning
 *
 * Uses total field F and inclination I to estimate position via grid search
 * + gradient descent on IGRF model. Accuracy ~100-500 km.
 * ======================================================================== */

typedef struct {
    const IGRFModel *model;
    double F_meas, I_meas, w_I, alt;
} IGRFPosCtx;

static double igrf_pos_cost_func(double lat, double lon, void *ctx_ptr)
{
    IGRFPosCtx *ctx = (IGRFPosCtx *)ctx_ptr;
    GeodeticCoord loc = { lat, lon, ctx->alt };
    MagVector B; MagneticElements elem;
    if (igrf_compute_field(ctx->model, &loc, &B) != 0) return 1e100;
    compute_magnetic_elements(&B, &elem);
    double cost = (ctx->F_meas - elem.total_intensity)
                * (ctx->F_meas - elem.total_intensity);
    if (ctx->w_I > 0.0)
        cost += ctx->w_I * (ctx->I_meas - elem.inclination)
                         * (ctx->I_meas - elem.inclination);
    return cost;
}

int igrf_single_point_position(const IGRFModel *model,
                                double F_measured, double I_measured,
                                double weight_I, double alt,
                                double init_lat, double init_lon,
                                double search_radius,
                                double *est_lat, double *est_lon,
                                double *residual)
{
    if (!model || !est_lat || !est_lon || !residual) return -1;
    IGRFPosCtx ctx = { model, F_measured, I_measured, weight_I, alt };
    double best_cost = 1e100;
    double best_lat = init_lat, best_lon = init_lon;
    int nsteps = 10;
    double step = search_radius / nsteps;

    for (int i = -nsteps; i <= nsteps; i++) {
        double lat = init_lat + i * step;
        for (int j = -nsteps; j <= nsteps; j++) {
            double lon = init_lon + j * step;
            double cost = igrf_pos_cost_func(lat, lon, &ctx);
            if (cost < best_cost) {
                best_cost = cost; best_lat = lat; best_lon = lon;
            }
        }
    }

    double gd_step = search_radius / 20.0;
    for (int iter = 0; iter < 20; iter++) {
        double eps = 0.05;
        double dfdlat = (igrf_pos_cost_func(best_lat+eps, best_lon, &ctx)
                       - igrf_pos_cost_func(best_lat-eps, best_lon, &ctx)) / (2.0*eps);
        double dfdlon = (igrf_pos_cost_func(best_lat, best_lon+eps, &ctx)
                       - igrf_pos_cost_func(best_lat, best_lon-eps, &ctx)) / (2.0*eps);
        double gn = sqrt(dfdlat*dfdlat + dfdlon*dfdlon);
        if (gn < 1e-12) break;
        double nl = best_lat - gd_step*dfdlat/gn;
        double no = best_lon - gd_step*dfdlon/gn;
        double nc = igrf_pos_cost_func(nl, no, &ctx);
        if (nc < best_cost) { best_cost = nc; best_lat = nl; best_lon = no; }
        else gd_step *= 0.5;
    }
    *est_lat = best_lat; *est_lon = best_lon; *residual = best_cost;
    return 0;
}

/* ========================================================================
 * L6: Dead reckoning with magnetic compass heading + speed
 * ======================================================================== */
void mag_compass_dr_update(NavSolution *state, double heading,
                            double speed, double dt)
{
    if (!state || dt <= 0.0) return;
    double lat_rad = state->position.lat * DEG2RAD;
    double hdg_rad = heading * DEG2RAD;
    double dist = speed * dt;
    double sin_lat = sin(lat_rad);
    double Rm = WGS84_A * (1.0 - WGS84_E2)
              / pow(1.0 - WGS84_E2 * sin_lat * sin_lat, 1.5);
    double Rn = WGS84_A / sqrt(1.0 - WGS84_E2 * sin_lat * sin_lat);
    state->position.lat += dist * cos(hdg_rad) / Rm * RAD2DEG;
    state->position.lon += dist * sin(hdg_rad) / (Rn * cos(lat_rad)) * RAD2DEG;
    state->timestamp += dt;
}

/* L5: Level magnetic heading */
double mag_heading_from_triaxial(const MagVector *B_body, double declination)
{
    double mag_rad = atan2(-B_body->by, B_body->bx);
    return mag_to_true_heading(mag_rad * RAD2DEG, declination);
}

/* L5: Tilt-compensated heading */
double tilt_compensated_heading(const MagVector *B_body,
                                 double roll, double pitch,
                                 double declination)
{
    double sr = sin(roll), cr = cos(roll);
    double sp = sin(pitch), cp = cos(pitch);
    double B_hx = B_body->bx * cp + B_body->by * sr * sp + B_body->bz * cr * sp;
    double B_hy = B_body->by * cr - B_body->bz * sr;
    double mag_rad = atan2(-B_hy, B_hx);
    return mag_to_true_heading(mag_rad * RAD2DEG, declination);
}

/* ========================================================================
 * L6: Magnetic gradient-based position update
 *
 * Kalman innovation: dz = F_meas - F_pred
 * K = P*H^T / (H*P*H^T + R)
 * dx = K*dz, P_post = (I-K*H)*P_prior
 * ======================================================================== */
void magnetic_gradient_update(double F_measured, double F_predicted,
                               const double gradient[3],
                               const double P_prior[9], double R_meas,
                               double dx[3], double P_post[9])
{
    double dz = F_measured - F_predicted;
    double HP[3];
    for (int i = 0; i < 3; i++)
        HP[i] = P_prior[i*3+0]*gradient[0] + P_prior[i*3+1]*gradient[1]
              + P_prior[i*3+2]*gradient[2];

    double S = gradient[0]*HP[0] + gradient[1]*HP[1] + gradient[2]*HP[2] + R_meas;
    if (S < 1e-15) {
        dx[0] = dx[1] = dx[2] = 0.0;
        memcpy(P_post, P_prior, 9 * sizeof(double));
        return;
    }

    double K[3] = { HP[0]/S, HP[1]/S, HP[2]/S };
    dx[0] = K[0]*dz; dx[1] = K[1]*dz; dx[2] = K[2]*dz;

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            double KH = K[i]*gradient[0]*P_prior[0*3+j]
                      + K[i]*gradient[1]*P_prior[1*3+j]
                      + K[i]*gradient[2]*P_prior[2*3+j];
            P_post[i*3+j] = P_prior[i*3+j] - KH;
        }
}

/* L5: Spatial gradient of F from IGRF via central differences */
int compute_field_spatial_gradient(const IGRFModel *model,
                                    const GeodeticCoord *loc,
                                    double eps_deg, double dh,
                                    double gradient[3])
{
    if (!model || !loc || !gradient) return -1;
    double lat_rad = loc->lat * DEG2RAD, R = GEOMAG_EARTH_RADIUS;
    GeodeticCoord p1, p2; MagVector B; MagneticElements elem;

    p1 = *loc; p2 = *loc; p1.lat += eps_deg; p2.lat -= eps_deg;
    igrf_compute_field(model, &p1, &B); compute_magnetic_elements(&B, &elem);
    double F1 = elem.total_intensity;
    igrf_compute_field(model, &p2, &B); compute_magnetic_elements(&B, &elem);
    double F2 = elem.total_intensity;
    gradient[0] = (F1-F2) / (2.0*eps_deg*DEG2RAD*R);

    p1 = *loc; p2 = *loc; p1.lon += eps_deg; p2.lon -= eps_deg;
    igrf_compute_field(model, &p1, &B); compute_magnetic_elements(&B, &elem);
    F1 = elem.total_intensity;
    igrf_compute_field(model, &p2, &B); compute_magnetic_elements(&B, &elem);
    F2 = elem.total_intensity;
    gradient[1] = (F1-F2) / (2.0*eps_deg*DEG2RAD*R*cos(lat_rad));

    p1 = *loc; p2 = *loc; p1.alt += dh; p2.alt -= dh;
    igrf_compute_field(model, &p1, &B); compute_magnetic_elements(&B, &elem);
    F1 = elem.total_intensity;
    igrf_compute_field(model, &p2, &B); compute_magnetic_elements(&B, &elem);
    F2 = elem.total_intensity;
    gradient[2] = (F1-F2) / (2.0*dh);
    return 0;
}

/* L6: Map roughness metric (spatial std dev for navigation quality) */
int magmap_roughness_metric(const MagneticMap *map,
                             const GeodeticCoord *center,
                             double radius, double *roughness)
{
    if (!map || !center || !roughness) return -1;
    double sum = 0.0, sum_sq = 0.0;
    int count = 0, n = 20;
    for (int i = -n; i <= n; i++) {
        double dlat = i * radius / n;
        for (int j = -n; j <= n; j++) {
            double dlon = j * radius / n;
            if (dlat*dlat + dlon*dlon > radius*radius) continue;
            GeodeticCoord pt = { center->lat+dlat, center->lon+dlon, 0.0 };
            double val;
            if (magmap_bilinear_interpolate(map, &pt, &val) == 0) {
                sum += val; sum_sq += val*val; count++;
            }
        }
    }
    if (count < 4) return -1;
    double mean = sum / count;
    *roughness = sqrt(sum_sq/count - mean*mean);
    return 0;
}
