/**
 * geomag_model.h -- Geomagnetic Field Models (IGRF/WMM)
 *
 * L2: Core concept -- Earth magnetic field as spherical harmonic expansion
 * L4: Laplace equation solution in source-free region, dipole approximation
 * L5: IGRF synthesis algorithm via Schmidt semi-normalized Legendre polynomials
 *
 * The Earth main magnetic field is the negative gradient of a scalar
 * potential V satisfying Laplace equation above Earth surface:
 *
 *   V(r,theta,phi) = a * sum_{n=1}^{Nmax} (a/r)^{n+1} *
 *                    sum_{m=0}^{n} [g(n,m) cos(m phi) + h(n,m) sin(m phi)] *
 *                    P(n,m)(cos theta)
 *
 *   B = -grad(V), giving components:
 *   B_r = sum (a/r)^{n+2} (n+1) * term
 *   B_theta = -(1/r) dV/dtheta
 *   B_phi = -(1/(r sin theta)) dV/dphi
 *
 * Reference: Alken et al., "IGRF-13", Earth Planets Space (2021)
 *            Langel, "The Main Field", Geomagnetism Vol.1 (1987)
 *            NOAA/NCEI, "The World Magnetic Model 2020 Technical Report"
 */

#ifndef GEOMAG_MODEL_H
#define GEOMAG_MODEL_H

#include "geomag_core.h"

/**
 * Compute flat array index for Gauss coefficients at degree n, order m.
 * formula: n*(n+1)/2 + m  (natural index; use -1 for Gauss coeff array)
 */
int schmidt_index(int n, int m);

/**
 * L5: Compute Schmidt semi-normalized Legendre functions P_n^m(cos theta)
 * and their theta-derivatives up to degree nmax.
 *
 * Schmidt normalization: P_n^m_schmidt = S(n,m) * P_n^m_unnormalized
 *   S(n,0) = 1.0
 *   S(n,m) = sqrt(2 * (n-m)! / (n+m)!)  for m>0
 *
 * Recurrence (forward column method):
 *   P(0,0) = 1
 *   P(m,m) = sqrt((2m+1)/(2m)) * sin(theta) * P(m-1,m-1)  for m>=1
 *   P(n,m) = [(2n-1)*cos(theta)*P(n-1,m) - (n+m-1)*P(n-2,m)] / (n-m)  for n>m
 *
 * dP/dtheta: via relation using P(n,m) and P(n-1,m)
 *
 * Complexity: O(nmax^2) time, O(nmax^2) space.
 *
 * @param nmax       Maximum degree
 * @param colatitude Colatitude theta [rad], 0 = geographic North pole
 * @param state      Output Legendre state (caller-allocated)
 */
void compute_schmidt_legendre(int nmax, double colatitude, LegendreState *state);

LegendreState *alloc_legendre_state(int nmax);
void free_legendre_state(LegendreState *state);

/**
 * L4+L5: Compute geomagnetic field at a point from IGRF model.
 *
 * Steps:
 *   1. Convert geodetic (lat,lon,h) to geocentric (r,theta,phi)
 *      theta_colat = PI/2 - geocentric_latitude
 *      r = sqrt(x^2+y^2+z^2) from ECEF
 *   2. Compute P(n,m)(cos theta) and dP/dtheta via Schmidt recurrence
 *   3. Compute partial sums for B_r, B_theta, B_phi
 *      B_r = -sum_{n=1}^{nmax} (a/r)^{n+2} * (n+1) *
 *            sum_{m=0}^{n} [g(n,m)*cos(m phi) + h(n,m)*sin(m phi)] * P(n,m)
 *      B_theta = sum (a/r)^{n+2} * sum [g cos(m phi)+h sin(m phi)] * dP/dtheta
 *      B_phi = -(1/sin theta) * sum (a/r)^{n+2} *
 *              sum m*[-g sin(m phi)+h cos(m phi)] * P(n,m)
 *   4. Transform geocentric spherical -> geodetic NED via rotation by
 *      (geodetic_lat - geocentric_lat)
 *
 * @param model  IGRF model with loaded coefficients
 * @param loc    Geodetic coordinates [deg, deg, m]
 * @param B      Output field vector in NED frame [nT]
 * @return 0 on success, -1 on error
 */
int igrf_compute_field(const IGRFModel *model, const GeodeticCoord *loc, MagVector *B);

/**
 * L5: Compute secular variation dB/dt at a point.
 *
 * Same algorithm as igrf_compute_field but substituting dg/dh for g/h.
 *
 * @param model  IGRF model with SV coefficients
 * @param loc    Geodetic coordinates
 * @param dB_dt  Output secular variation vector [nT/yr]
 * @return 0 on success
 */
int igrf_compute_secular_variation(const IGRFModel *model,
                                    const GeodeticCoord *loc, MagVector *dB_dt);

/**
 * L5: Predict field at future date using linear secular variation.
 * B(t) = B(epoch) + (dB/dt) * (t - epoch)
 *
 * Valid within +/-5 years of model epoch for IGRF.
 *
 * @param model  IGRF model
 * @param loc    Geodetic coordinates
 * @param year   Decimal year (e.g., 2023.5)
 * @param B      Output predicted field [nT]
 * @return 0 on success
 */
int igrf_predict_field(const IGRFModel *model, const GeodeticCoord *loc,
                        double year, MagVector *B);

/**
 * L4: Geocentric tilted dipole approximation.
 *
 * Uses only n=1 coefficients. Accounts for ~90% of total field.
 * Fast O(1) computation for coarse navigation/compass applications.
 *
 * Dipole moment vector from (g10, g11, h11).
 *
 * @param g10, g11, h11  Dipole Gauss coefficients [nT]
 * @param loc             Geodetic coordinates
 * @param B               Output dipole field in NED [nT]
 * @return 0 on success
 */
int compute_dipole_field(double g10, double g11, double h11,
                          const GeodeticCoord *loc, MagVector *B);

/**
 * L5: Find magnetic poles (where inclination = +/-90 degrees).
 *
 * Iterative search on a coarse grid refined by gradient descent.
 *
 * @param model       IGRF model
 * @param north_pole  Output magnetic north dip pole [deg]
 * @param south_pole  Output magnetic south dip pole [deg]
 * @return 0 on success
 */
int compute_magnetic_poles(const IGRFModel *model,
                            GeodeticCoord *north_pole, GeodeticCoord *south_pole);

/**
 * L5: Compute Earth dipole magnetic moment.
 * m = (4*pi*a^3/mu_0) * sqrt(g10^2 + g11^2 + h11^2)
 *
 * @param model  IGRF model
 * @return Dipole moment [A*m^2]
 */
double compute_dipole_moment(const IGRFModel *model);

/**
 * L5: Compute McIlwain L-shell parameter using dipole approximation.
 * L = r / (a * cos^2(lambda_magnetic))
 *
 * Used for radiation belt and auroral zone mapping.
 *
 * @param g10, g11, h11  Dipole coefficients [nT]
 * @param loc             Geodetic coordinates
 * @return L-value [dimensionless, Earth radii]
 */
double compute_l_shell(double g10, double g11, double h11,
                        const GeodeticCoord *loc);

/**
 * L5: Compute total magnetic field energy outside core.
 * W = (4*pi*a^3/mu_0) * sum_{n=1}^{Nmax} (n+1) * sum_{m=0}^{n} [(g_n^m)^2 + (h_n^m)^2]
 *
 * @param model  IGRF model
 * @return Field energy [J]
 */
double compute_field_energy(const IGRFModel *model);

/**
 * L6: South Atlantic Anomaly proximity indicator.
 *
 * Returns normalized "closeness" metric (0 = far, 1 = deep in SAA).
 * Based on total field depression relative to dipole.
 *
 * @param model  IGRF model
 * @param loc    Location to evaluate
 * @return SAA indicator [0-1]
 */
double compute_saa_indicator(const IGRFModel *model, const GeodeticCoord *loc);

/**
 * Initialize IGRF model with built-in coefficients for a given version.
 *
 * Coefficients are compiled from IAGA IGRF-13 dataset for supported epochs.
 * Currently supports IGRF-13 (epoch 2020.0) with full n=1..13 coefficients.
 *
 * @param version  IGRF version number (13 supported)
 * @param model    Output model (caller allocates memory for struct, we fill coeffs)
 * @return 0 on success, -1 on unsupported version
 */
int igrf_init_model(int version, IGRFModel *model);

/**
 * Free IGRF model coefficient memory.
 */
void igrf_free_model(IGRFModel *model);

#endif /* GEOMAG_MODEL_H */
