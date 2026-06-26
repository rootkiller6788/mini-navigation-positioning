/**
 * geomag_navigation.h -- Magnetic Navigation Algorithms
 *
 * L5: Magnetic map matching (MAGCOM), gradient navigation
 * L6: Underwater/underground positioning, compass navigation
 * L7: GPS-denied navigation applications
 *
 * Magnetic navigation uses spatial variations in Earth magnetic field
 * as a positioning signal. Key methods:
 *
 * 1. MAGCOM (Magnetic Contour Matching):
 *    Correlates measured magnetic profile along trajectory against
 *    a pre-surveyed magnetic map (analogous to TERCOM for terrain).
 *    Accuracy: 100-1000 m with good maps, limited by map resolution.
 *
 * 2. Magnetic Gradient Navigation:
 *    Uses spatial gradient of the magnetic field for dead-reckoning
 *    correction. Gradient measurements are insensitive to diurnal
 *    variations (common-mode rejection).
 *    Accuracy: ~10-100 m drift per hour when integrated with INS.
 *
 * 3. Magnetic Compass + Map:
 *    Compass provides heading; combined with magnetic inclination
 *    and total field measurements for position estimation via
 *    matching against global models (IGRF/WMM) or local maps.
 *
 * Reference:
 *   Goldenberg, "Geomagnetic Navigation beyond the Magnetic Compass",
 *     IEEE PLANS (2006)
 *   Canciani & Raquet, "Absolute Positioning Using the Earth Magnetic
 *     Field", NAVIGATION (2016)
 *   Shockley & Raquet, "Navigation Using Magnetic Field Gradients",
 *     IEEE Trans. Aerospace (2014)
 */

#ifndef GEOMAG_NAVIGATION_H
#define GEOMAG_NAVIGATION_H

#include "geomag_core.h"

/* ========================================================================
 * L6: Magnetic Contour Matching (MAGCOM)
 * ======================================================================== */

/**
 * L6: MAGCOM positioning via correlation-based map matching.
 *
 * Algorithm:
 *   1. Collect a sequence of N magnetic measurements along trajectory
 *      (total field, or total field + inclination + declination)
 *   2. Slide the measurement profile across the magnetic map
 *   3. Compute correlation at each candidate position
 *   4. Position estimate = location of maximum correlation
 *
 * Correlation metric: Normalized Cross-Correlation (NCC)
 *   NCC = sum[(z_i - z_bar)*(m_i(p) - m_bar(p))] /
 *         sqrt(sum[(z_i-z_bar)^2] * sum[(m_i(p)-m_bar(p))^2])
 *
 * @param map           Pre-loaded magnetic map
 * @param measurements  Array of measured total field values [nT]
 * @param trajectory_lat Array of latitude offsets along trajectory [deg]
 * @param trajectory_lon Array of longitude offsets along trajectory [deg]
 * @param N             Number of measurements in profile
 * @param search_radius Search radius around initial position [deg]
 * @param search_step   Search step size [deg]
 * @param est_lat       Output estimated latitude [deg]
 * @param est_lon       Output estimated longitude [deg]
 * @param confidence    Output match confidence [0-1]
 * @return 0 on success, -1 on error
 */
int magcom_correlation_match(const MagneticMap *map,
                              const double *measurements,
                              const double *trajectory_lat,
                              const double *trajectory_lon,
                              int N,
                              double search_radius,
                              double search_step,
                              double *est_lat, double *est_lon,
                              double *confidence);

/**
 * L5: Mean Absolute Difference (MAD) map matching metric.
 *
 * MAD(x,y) = (1/N) * sum_{i=1}^{N} |z_i - m_i(x,y)|
 *
 * Position estimate = argmin MAD(x,y)
 *
 * More robust to outliers than correlation-based matching.
 *
 * @param map           Magnetic map
 * @param measurements  Measured total field [nT]
 * @param trajectory_lat Trajectory latitude offsets [deg]
 * @param trajectory_lon Trajectory longitude offsets [deg]
 * @param N             Profile length
 * @param search_radius Search radius [deg]
 * @param search_step   Search step [deg]
 * @param est_lat       Output latitude estimate [deg]
 * @param est_lon       Output longitude estimate [deg]
 * @param min_mad       Output minimum MAD value
 * @return 0 on success
 */
int magcom_mad_match(const MagneticMap *map,
                      const double *measurements,
                      const double *trajectory_lat,
                      const double *trajectory_lon,
                      int N,
                      double search_radius, double search_step,
                      double *est_lat, double *est_lon,
                      double *min_mad);

/**
 * L5: Magnetic field profile generation from map along a trajectory.
 *
 * Extracts the expected magnetic field measurements that would be observed
 * along a given trajectory using bilinear interpolation on the map.
 *
 * @param map           Magnetic map
 * @param start         Start position
 * @param heading       Trajectory heading [deg]
 * @param step_size     Distance between samples [m]
 * @param N             Number of samples
 * @param profile       Output: expected total field values [nT]
 * @param profile_lats  Output: latitudes of sample points [deg]
 * @param profile_lons  Output: longitudes of sample points [deg]
 * @return 0 on success
 */
int generate_magnetic_profile(const MagneticMap *map,
                               const GeodeticCoord *start,
                               double heading, double step_size,
                               int N,
                               double *profile,
                               double *profile_lats,
                               double *profile_lons);

/* ========================================================================
 * L6: IGRF-based single-point positioning
 * ======================================================================== */

/**
 * L6: Single-point magnetic positioning using IGRF model.
 *
 * Given measurements of total field (F) and inclination (I),
 * find the position that best matches these scalar values using
 * the IGRF global model as reference map.
 *
 * This is a 2D search problem: find (lat, lon) minimizing
 *   cost(lat,lon) = (F_measured - F_igrf(lat,lon))^2
 *                  + w * (I_measured - I_igrf(lat,lon))^2
 *
 * Uses coarse grid search followed by gradient descent refinement.
 *
 * Accuracy: ~100-500 km (limited by low spatial gradient of main field).
 *
 * @param model       IGRF model
 * @param F_measured  Measured total field [nT]
 * @param I_measured  Measured inclination [deg] (NaN if unavailable)
 * @param weight_I    Weight for inclination term (0 if not using)
 * @param alt         Assumed altitude [m]
 * @param init_lat,init_lon Initial position guess (search center) [deg]
 * @param search_radius Search radius [deg]
 * @param est_lat,est_lon Output estimated position [deg]
 * @param residual    Output final cost value
 * @return 0 on success
 */
int igrf_single_point_position(const IGRFModel *model,
                                double F_measured, double I_measured,
                                double weight_I, double alt,
                                double init_lat, double init_lon,
                                double search_radius,
                                double *est_lat, double *est_lon,
                                double *residual);

/* ========================================================================
 * L6: Magnetic compass navigation
 * ======================================================================== */

/**
 * L6: Dead reckoning using magnetic compass heading and speed.
 *
 * Updates position given heading (from compass, corrected for declination),
 * speed over ground, and time step.
 *
 * lat_new = lat + (speed * dt * cos(heading)) / R_earth
 * lon_new = lon + (speed * dt * sin(heading)) / (R_earth * cos(lat))
 *
 * @param state     Input/output navigation state
 * @param heading   True heading [deg]
 * @param speed     Speed over ground [m/s]
 * @param dt        Time step [s]
 */
void mag_compass_dr_update(NavSolution *state, double heading,
                            double speed, double dt);

/**
 * L5: Estimate heading from triaxial magnetometer in level attitude.
 *
 * heading_mag = atan2(-By, Bx)   [rad]
 * heading_true = heading_mag + declination
 *
 * where Bx, By are horizontal components in body frame.
 * This formula assumes the sensor is level (roll=pitch=0).
 *
 * @param B_body      Magnetometer reading in body frame [nT]
 * @param declination Local declination [deg]
 * @return True heading [deg]
 */
double mag_heading_from_triaxial(const MagVector *B_body, double declination);

/**
 * L5: Tilt-compensated magnetic heading.
 *
 * For non-level sensor, project magnetometer readings to horizontal
 * plane using attitude (roll, pitch) from accelerometer or INS.
 *
 * B_horiz_x = Bx*cos(pitch) + By*sin(roll)*sin(pitch) + Bz*cos(roll)*sin(pitch)
 * B_horiz_y = By*cos(roll) - Bz*sin(roll)
 * heading = atan2(-B_horiz_y, B_horiz_x)
 *
 * @param B_body    Magnetometer reading [nT]
 * @param roll      Roll angle [rad]
 * @param pitch     Pitch angle [rad]
 * @param declination Local declination [deg]
 * @return Tilt-compensated true heading [deg]
 */
double tilt_compensated_heading(const MagVector *B_body,
                                 double roll, double pitch,
                                 double declination);

/* ========================================================================
 * L6: Magnetic gradient-based positioning
 * ======================================================================== */

/**
 * L6: Gradient-based magnetic navigation update.
 *
 * Uses spatial gradient of the magnetic field to correct dead-reckoning
 * drift. The innovation is the difference between measured and predicted
 * field values.
 *
 * Concept:
 *   z - h(x) = H * dx + v
 *   where H = [dF/dx, dF/dy, dF/dz] is the gradient vector
 *         dx = position correction
 *
 * Used within an EKF framework for integrated INS/MAG navigation.
 *
 * @param F_measured    Measured total field [nT]
 * @param F_predicted   Predicted total field [nT]
 * @param gradient      Spatial gradient of F [nT/m] (3 components in NED)
 * @param P_prior       Prior position error covariance (3x3)
 * @param R_meas        Measurement noise variance
 * @param dx            Output position correction [m] (NED)
 * @param P_post        Output posterior covariance (3x3)
 */
void magnetic_gradient_update(double F_measured, double F_predicted,
                               const double gradient[3],
                               const double P_prior[9], double R_meas,
                               double dx[3], double P_post[9]);

/**
 * L5: Compute spatial gradient of total field from IGRF model.
 *
 * Numerical gradient via central differences on the WGS84 ellipsoid.
 *
 * dF/dx ≈ (F(lat+eps, lon) - F(lat-eps, lon)) / (2*eps*R_Earth*cos(lat))
 * dF/dy ≈ (F(lat, lon+eps) - F(lat, lon-eps)) / (2*eps*R_Earth*cos(lat))
 * dF/dz ≈ (F(lat, lon, alt+dh) - F(lat, lon, alt-dh)) / (2*dh)
 *
 * @param model    IGRF model
 * @param loc      Location for gradient computation
 * @param eps_deg  Perturbation [deg]
 * @param dh       Altitude perturbation [m]
 * @param gradient Output gradient [nT/m] in NED frame (3 components)
 * @return 0 on success
 */
int compute_field_spatial_gradient(const IGRFModel *model,
                                    const GeodeticCoord *loc,
                                    double eps_deg, double dh,
                                    double gradient[3]);

/**
 * L6: Evaluate map-matching quality (information content of a map region).
 *
 * Computes the spatial roughness (standard deviation of field) in a
 * region, which determines navigation accuracy.
 *
 * High roughness (complex magnetic terrain) --> better positioning
 * Low roughness (flat magnetic field) --> poor positioning
 *
 * @param map      Magnetic map
 * @param center   Center of region
 * @param radius   Radius of region [deg]
 * @param roughness Output: std dev of total field [nT]
 * @return 0 on success
 */
int magmap_roughness_metric(const MagneticMap *map,
                             const GeodeticCoord *center,
                             double radius,
                             double *roughness);

#endif /* GEOMAG_NAVIGATION_H */
