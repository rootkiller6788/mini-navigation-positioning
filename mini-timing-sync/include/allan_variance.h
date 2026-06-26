/**
 * @file allan_variance.h
 * @brief Allan variance and frequency stability analysis
 *
 * L1: Allan variance, Allan deviation, Modified Allan, Hadamard variance
 * L2: Noise type identification (white PM, flicker PM, white FM, flicker FM, random walk FM)
 * L3: Power-law spectral density model, sigma-tau plots
 * L4: Allan variance formula, relation to power spectral density
 * L5: Overlapping Allan variance computation, noise coefficient estimation
 * L6: Oscillator stability characterization
 * L7: GPS-disciplined oscillator (GPSDO) stability analysis
 * L8: Three-cornered hat method for separating oscillator noises
 *
 * Reference: IEEE 1139-2008 (Standard Definitions of Frequency Stability)
 *            Riley, W.J. "Handbook of Frequency Stability Analysis" (NIST SP 1065)
 *            Allan, D.W. "Statistics of Atomic Frequency Standards" (1966)
 *
 * Course: MIT 6.450, Stanford EE359, TU Munich High-Frequency Eng.
 */

#ifndef ALLAN_VARIANCE_H
#define ALLAN_VARIANCE_H
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif

/* L1: Allan variance result structure */
typedef struct {
    double tau_s;              /* Averaging time [seconds] */
    double adev;               /* Allan deviation (sqrt of AVAR) */
    double avar;               /* Allan variance */
    double mdev;               /* Modified Allan deviation */
    double hdev;               /* Hadamard deviation */
    double adev_uncertainty;   /* 1-sigma uncertainty of ADEV */
    int    num_pairs;          /* Number of frequency pairs used */
} AllanResult;

/* L1: Noise type classification */
typedef enum {
    NOISE_WHITE_PM        = 0,  /* White Phase Modulation (tau^-1) */
    NOISE_FLICKER_PM      = 1,  /* Flicker Phase Modulation (tau^-1) */
    NOISE_WHITE_FM        = 2,  /* White Frequency Modulation (tau^-1/2) */
    NOISE_FLICKER_FM      = 3,  /* Flicker Frequency Modulation (tau^0) */
    NOISE_RANDOM_WALK_FM  = 4,  /* Random Walk Frequency Modulation (tau^+1/2) */
    NOISE_FREQ_DRIFT      = 5   /* Linear frequency drift (tau^+1) */
} NoiseType;

/* L1: Power-law noise coefficients */
typedef struct {
    double h_2;   /* White PM coefficient (alpha = +2) */
    double h_1;   /* Flicker PM coefficient (alpha = +1) */
    double h_0;   /* White FM coefficient (alpha = 0) */
    double h_m1;  /* Flicker FM coefficient (alpha = -1) */
    double h_m2;  /* Random walk FM coefficient (alpha = -2) */
} NoiseCoefficients;

/**
 * L4: Compute overlapping Allan variance.
 *
 * For N fractional frequency measurements y[i] spaced by tau0:
 *   sigma_y^2(tau) = 1/(2*m^2*(M-2*m+1)) * sum_{j=1}^{M-2*m+1} [sum_{i=j}^{j+m-1} (y[i+m] - y[i])]^2
 *
 * where m = tau/tau0 (averaging factor) and M = N (total samples).
 *
 * @param yi          Array of fractional frequency measurements
 * @param N           Number of measurements
 * @param tau0_s      Basic measurement interval [seconds]
 * @param tau_s       Desired averaging time [seconds]
 * @param result      [out] Allan variance result
 * @return 0 on success, -1 if insufficient data
 */
int allan_variance_compute(const double *yi, int N, double tau0_s,
                           double tau_s, AllanResult *result);

/**
 * L5: Compute overlapping Allan deviation for multiple tau values.
 *
 * Generates the sigma-tau plot data points.
 *
 * @param yi           Fractional frequency measurements
 * @param N            Number of measurements
 * @param tau0_s       Basic measurement interval
 * @param taus_s       Array of desired averaging times
 * @param num_taus     Number of tau values
 * @param results      [out] Array of results (length = num_taus)
 * @return 0 on success, -1 on error
 */
int allan_deviation_multi_tau(const double *yi, int N, double tau0_s,
                              const double *taus_s, int num_taus,
                              AllanResult *results);

/**
 * L4: Compute modified Allan variance (MVAR).
 *
 * MVAR provides better discrimination between white PM and flicker PM:
 *   mod_sigma_y^2(tau) = 1/(2*m^4*(N-3*m+2)) *
 *     sum_{j=1}^{N-3*m+2} { sum_{i=j}^{j+m-1} [sum_{k=i}^{i+m-1} (y[k+m] - y[k])] }^2
 *
 * @param yi       Fractional frequency measurements
 * @param N        Number of measurements
 * @param tau0_s   Basic measurement interval
 * @param m        Averaging factor
 * @param mvar     [out] Modified Allan variance
 * @return 0 on success, -1 on error
 */
int modified_allan_variance(const double *yi, int N, double tau0_s,
                            int m, double *mvar);

/**
 * L5: Fit noise coefficients from Allan deviation data.
 *
 * Uses least-squares fit in log-log domain:
 *   log(sigma_y^2(tau)) = log(h_alpha) + alpha * log(tau)
 *
 * For each noise type, sigma_y^2(tau) ~ tau^{mu/2} where
 * mu = -alpha - 1 for FM noises or mu = -alpha - 2 for PM noises.
 *
 * @param taus_s     Array of averaging times
 * @param adevs      Array of Allan deviations
 * @param num_points Number of data points
 * @param coeffs     [out] Fitted noise coefficients
 * @return 0 on success
 */
int allan_fit_noise_coefficients(const double *taus_s,
                                 const double *adevs,
                                 int num_points,
                                 NoiseCoefficients *coeffs);

/**
 * L5: Identify dominant noise type at a given tau.
 *
 * Uses the slope of log(ADEV) vs log(tau):
 *   -1.0  slope -> White PM or Flicker PM
 *   -0.5  slope -> White FM
 *    0.0  slope -> Flicker FM
 *   +0.5  slope -> Random Walk FM
 *   +1.0  slope -> Frequency Drift
 *
 * @param adev1   ADEV at tau1
 * @param tau1_s  First averaging time
 * @param adev2   ADEV at tau2
 * @param tau2_s  Second averaging time
 * @return Identified noise type
 */
NoiseType allan_identify_noise_type(double adev1, double tau1_s,
                                    double adev2, double tau2_s);

/**
 * L6: Compute time deviation (TDEV) from phase measurements.
 *
 * TDEV is the time-domain stability measure preferred for
 * time transfer and synchronization applications.
 * TDEV(tau) = (tau/sqrt(3)) * mod_sigma_y(tau)
 *
 * @param xi           Phase (time) measurements [seconds]
 * @param N            Number of measurements
 * @param tau0_s       Basic measurement interval
 * @param tau_s        Desired averaging time
 * @param tdev         [out] Time deviation in seconds
 * @return 0 on success
 */
int time_deviation_compute(const double *xi, int N, double tau0_s,
                           double tau_s, double *tdev);

/**
 * L6: Compute Hadamard variance for frequency data.
 *
 * Hadamard variance is insensitive to linear frequency drift,
 * making it useful for characterizing cesium and rubidium standards.
 *
 * @param yi       Fractional frequency measurements
 * @param N        Number of measurements
 * @param tau0_s   Basic measurement interval
 * @param tau_s    Averaging time
 * @param hadamard [out] Hadamard deviation
 * @return 0 on success
 */
int hadamard_variance(const double *yi, int N, double tau0_s,
                      double tau_s, double *hadamard);

/**
 * L7: GPSDO stability specification compliance check.
 *
 * Typical GPSDO performance:
 * - ADEV(1s) < 2e-11
 * - ADEV(100s) < 1e-12
 * - ADEV(10000s) < 5e-13
 *
 * @param adev_1s     Allan deviation at 1 s
 * @param adev_100s   Allan deviation at 100 s
 * @param adev_10000s Allan deviation at 10000 s
 * @param grade       'S' (standard), 'P' (precision), 'U' (ultra-stable)
 * @return 1 if compliant with grade specification
 */
int gpsdo_stability_check(double adev_1s, double adev_100s,
                          double adev_10000s, char grade);

/**
 * L8: Three-cornered hat method for separating oscillator noises.
 *
 * Given three oscillators A, B, C with pairwise Allan variances:
 *   sigma_AB^2, sigma_BC^2, sigma_CA^2
 *
 * Individual variances estimated as:
 *   sigma_A^2 = (sigma_AB^2 + sigma_CA^2 - sigma_BC^2) / 2
 *   sigma_B^2 = (sigma_AB^2 + sigma_BC^2 - sigma_CA^2) / 2
 *   sigma_C^2 = (sigma_BC^2 + sigma_CA^2 - sigma_AB^2) / 2
 *
 * @param sigma_AB2  Allan variance of A vs B
 * @param sigma_BC2  Allan variance of B vs C
 * @param sigma_CA2  Allan variance of C vs A
 * @param sigma_A2   [out] Individual variance of oscillator A
 * @param sigma_B2   [out] Individual variance of oscillator B
 * @param sigma_C2   [out] Individual variance of oscillator C
 * @return 0 on success, -1 if triangle inequality violated
 */
int three_cornered_hat(double sigma_AB2, double sigma_BC2, double sigma_CA2,
                       double *sigma_A2, double *sigma_B2, double *sigma_C2);

#ifdef __cplusplus
}
#endif
#endif
