/**
 * @file nav_utils.c
 * @brief Navigation Utility Functions
 *
 * L2: Navigation data processing utilities.
 * L5: Outlier detection, moving average, data validation.
 */

#include "nav_common.h"
#include "nav_rotation.h"
#include <math.h>
#include <string.h>

/**
 * @brief 3-point moving average filter for IMU data smoothing.
 */
void nav_moving_average_3(NAV_PRECISION out[3],
                           const NAV_PRECISION prev[3],
                           const NAV_PRECISION curr[3],
                           const NAV_PRECISION next[3]) {
    if (!out || !prev || !curr || !next) return;
    for (int i = 0; i < 3; i++)
        out[i] = (prev[i] + curr[i] + next[i]) / 3.0;
}

/**
 * @brief Simple outlier detection using median absolute deviation.
 * Returns 1 if any element exceeds median + k * MAD, 0 otherwise.
 */
int nav_outlier_detect(const NAV_PRECISION data[3],
                        const NAV_PRECISION history[][3], int n_hist,
                        NAV_PRECISION k) {
    if (!data || !history || n_hist < 3) return 0;
    /* Compute component-wise median of recent history */
    NAV_PRECISION med[3];
    for (int c = 0; c < 3; c++) {
        NAV_PRECISION vals[16];
        int nv = (n_hist < 16) ? n_hist : 16;
        for (int i = 0; i < nv; i++)
            vals[i] = history[i][c];
        /* Simple bubble sort for median */
        for (int i = 0; i < nv-1; i++)
            for (int j = i+1; j < nv; j++)
                if (vals[i] > vals[j]) {
                    NAV_PRECISION t = vals[i];
                    vals[i] = vals[j];
                    vals[j] = t;
                }
        med[c] = vals[nv/2];
    }
    /* MAD */
    NAV_PRECISION mad[3] = {0,0,0};
    for (int c = 0; c < 3; c++) {
        NAV_PRECISION ad[16];
        int nv = (n_hist < 16) ? n_hist : 16;
        for (int i = 0; i < nv; i++)
            ad[i] = fabs(history[i][c] - med[c]);
        for (int i = 0; i < nv-1; i++)
            for (int j = i+1; j < nv; j++)
                if (ad[i] > ad[j]) {
                    NAV_PRECISION t = ad[i];
                    ad[i] = ad[j];
                    ad[j] = t;
                }
        mad[c] = ad[nv/2] * 1.4826; /* scale factor for normal distribution */
    }
    for (int c = 0; c < 3; c++)
        if (fabs(data[c] - med[c]) > k * (mad[c] + 1e-10))
            return 1;
    return 0;
}

/**
 * @brief Interpolate between two navigation solutions.
 * Linear interpolation for position and velocity; SLERP for attitude.
 */
void nav_solution_interpolate(nav_solution_t *out,
                               const nav_solution_t *a,
                               const nav_solution_t *b,
                               NAV_PRECISION t) {
    if (!out || !a || !b) return;
    NAV_PRECISION t1 = 1.0 - t;
    /* Position: linear interpolation on sphere */
    out->pos.latitude  = t1*a->pos.latitude  + t*b->pos.latitude;
    out->pos.longitude = t1*a->pos.longitude + t*b->pos.longitude;
    out->pos.altitude  = t1*a->pos.altitude  + t*b->pos.altitude;
    /* Velocity: linear */
    out->vel_ned.x = t1*a->vel_ned.x + t*b->vel_ned.x;
    out->vel_ned.y = t1*a->vel_ned.y + t*b->vel_ned.y;
    out->vel_ned.z = t1*a->vel_ned.z + t*b->vel_ned.z;
    /* Attitude: SLERP */
    nav_quat_t qa, qb, qr;
    nav_quat_set(&qa, a->quat[0], a->quat[1], a->quat[2], a->quat[3]);
    nav_quat_set(&qb, b->quat[0], b->quat[1], b->quat[2], b->quat[3]);
    nav_quat_slerp(&qr, &qa, &qb, t);
    out->quat[0] = qr.w; out->quat[1] = qr.x;
    out->quat[2] = qr.y; out->quat[3] = qr.z;
}

/**
 * @brief Compute horizontal position difference in meters.
 * Uses the Haversine formula for great-circle distance.
 */
NAV_PRECISION nav_horizontal_distance(NAV_PRECISION lat1, NAV_PRECISION lon1,
                                       NAV_PRECISION lat2, NAV_PRECISION lon2) {
    NAV_PRECISION dlat = lat2 - lat1;
    NAV_PRECISION dlon = lon2 - lon1;
    NAV_PRECISION a = sin(dlat/2.0)*sin(dlat/2.0) +
                      cos(lat1)*cos(lat2)*sin(dlon/2.0)*sin(dlon/2.0);
    NAV_PRECISION c = 2.0 * atan2(sqrt(a), sqrt(1.0 - a));
    return NAV_EARTH_RADIUS_M * c;
}

/**
 * @brief Compute course over ground (heading) from NED velocity.
 * @return heading in radians [0, 2*PI], 0 = North, PI/2 = East
 */
NAV_PRECISION nav_course_over_ground(NAV_PRECISION vn, NAV_PRECISION ve) {
    NAV_PRECISION cog = atan2(ve, vn);
    if (cog < 0.0) cog += 2.0 * NAV_PI;
    return cog;
}

/**
 * @brief Compute speed over ground from NED velocity, m/s.
 */
NAV_PRECISION nav_speed_over_ground(NAV_PRECISION vn, NAV_PRECISION ve,
                                      NAV_PRECISION vd) {
    return sqrt(vn*vn + ve*ve + vd*vd);
}
