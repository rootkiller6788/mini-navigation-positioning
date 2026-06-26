/**
 * @file allan_variance.c
 * @brief Allan variance and frequency stability analysis implementation
 *
 * Implements: Overlapping Allan variance, Modified Allan variance,
 * Hadamard variance, Time deviation, noise type identification,
 * noise coefficient fitting, GPSDO stability check,
 * three-cornered hat method.
 */

#include "allan_variance.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

/* ===================================================================
 * L4: Overlapping Allan Variance
 * =================================================================== */

int allan_variance_compute(const double *yi, int N, double tau0_s,
                           double tau_s, AllanResult *result)
{
    if (!yi || !result || N < 4) return -1;
    if (tau0_s <= 0.0 || tau_s <= 0.0) return -1;

    /* Averaging factor m = tau / tau0 */
    int m = (int)(tau_s / tau0_s + 0.5);
    if (m < 1) m = 1;
    if (m > N / 2) m = N / 2;

    /* Number of frequency pairs: M = N - 2*m + 1 */
    int M = N - 2 * m + 1;
    if (M < 2) return -1;

    /* Compute overlapping Allan variance:
     *   sigma_y^2(tau) = 1/(2*m^2*(M)) * sum_{j=1}^{M}
     *     [ (1/m)*sum_{i=j}^{j+m-1} y[i] - (1/m)*sum_{i=j+m}^{j+2m-1} y[i] ]^2
     *
     * Using blocks of size m:
     *   block1[j] = average(y[j]..y[j+m-1])
     *   block2[j] = average(y[j+m]..y[j+2*m-1])
     *   diff[j] = block2[j] - block1[j]
     */

    double sum_sq = 0.0;

    for (int j = 0; j < M; j++) {
        /* Average first block */
        double block1 = 0.0;
        for (int i = j; i < j + m; i++) {
            block1 += yi[i];
        }
        block1 /= (double)m;

        /* Average second block */
        double block2 = 0.0;
        for (int i = j + m; i < j + 2 * m; i++) {
            block2 += yi[i];
        }
        block2 /= (double)m;

        double diff = block2 - block1;
        sum_sq += diff * diff;
    }

    double avar = sum_sq / (2.0 * (double)M);

    result->tau_s = tau_s;
    result->avar = avar;
    result->adev = sqrt(avar);
    result->num_pairs = M;

    /* Uncertainty of ADEV: sigma_adev ~= adev / sqrt(M) */
    result->adev_uncertainty = result->adev / sqrt((double)M);
    result->mdev = 0.0;
    result->hdev = 0.0;

    return 0;
}

/* ===================================================================
 * L5: Multi-Tau Allan Deviation
 * =================================================================== */

int allan_deviation_multi_tau(const double *yi, int N, double tau0_s,
                              const double *taus_s, int num_taus,
                              AllanResult *results)
{
    if (!yi || !taus_s || !results) return -1;
    if (N < 4 || num_taus < 1) return -1;

    for (int t = 0; t < num_taus; t++) {
        int ret = allan_variance_compute(yi, N, tau0_s, taus_s[t], &results[t]);
        if (ret != 0) {
            /* If cannot compute for this tau, set to 0 */
            memset(&results[t], 0, sizeof(AllanResult));
            results[t].tau_s = taus_s[t];
        }
    }

    return 0;
}

/* ===================================================================
 * L4: Modified Allan Variance
 * =================================================================== */

int modified_allan_variance(const double *yi, int N, double tau0_s,
                            int m, double *mvar)
{
    if (!yi || !mvar || N < 4 || m < 1) return -1;
    (void)tau0_s; /* tau0 available for time-domain scaling */

    /* Modified Allan variance:
     *   mod_sigma_y^2(tau) = 1/(2*m^4*(N-3*m+2)) *
     *     sum_{j=1}^{N-3*m+2}
     *       { sum_{i=j}^{j+m-1} [ sum_{k=i}^{i+m-1} (y[k+m] - y[k]) ] }^2
     */

    int max_j = N - 3 * m + 2;
    if (max_j < 1) return -1;

    double sum_sq = 0.0;

    for (int j = 0; j < max_j; j++) {
        double inner_sum = 0.0;

        for (int i = j; i < j + m; i++) {
            double block_sum = 0.0;

            for (int k = i; k < i + m; k++) {
                block_sum += (yi[k + m] - yi[k]);
            }
            inner_sum += block_sum;
        }

        sum_sq += inner_sum * inner_sum;
    }

    *mvar = sum_sq / (2.0 * (double)m * m * m * m * (double)max_j);
    return 0;
}

/* ===================================================================
 * L5: Noise Coefficient Fitting
 * =================================================================== */

int allan_fit_noise_coefficients(const double *taus_s,
                                 const double *adevs,
                                 int num_points,
                                 NoiseCoefficients *coeffs)
{
    if (!taus_s || !adevs || !coeffs || num_points < 2) return -1;

    /* Initialize coefficients to zero */
    memset(coeffs, 0, sizeof(NoiseCoefficients));

    /* Fit using least squares in log-log domain.
     *
     * For each noise type alpha:
     *   sigma^2(tau) ~ h_alpha * tau^(mu/2)
     *   log(sigma^2) = log(h_alpha) + (mu/2) * log(tau)
     *
     * We identify the dominant noise by slope and fit each coefficient
     * at the tau where it dominates.
     *
     * Approach: piecewise fit using the known slope signatures:
     *   White PM (alpha=2):  sigma ~ tau^-1
     *   Flicker PM (alpha=1): sigma ~ tau^-1
     *   White FM (alpha=0):  sigma ~ tau^-1/2
     *   Flicker FM (alpha=-1): sigma ~ tau^0 (flat)
     *   Random Walk FM (alpha=-2): sigma ~ tau^+1/2
     */

    /* Simple approach: fit each noise type at its dominant tau region */

    /* White PM dominant at smallest tau (ADC/timestamp noise floor) */
    if (num_points >= 2) {
        double slope_pm = log10(adevs[1] / adevs[0])
                        / log10(taus_s[1] / taus_s[0]);
        /* White PM: slope = -1 in log-log ADEV */
        if (slope_pm < -0.75) {
            /* h_2: White PM, sigma^2 = 3*h_2*f_h^2/(2*pi^2*tau^2)
             * h_2 ~= (2/3)*pi^2*tau^2*sigma^2/f_h^2
             * Simplified estimate: h_2 ~= adev^2 * tau^2
             */
            double h_2_est = adevs[0] * adevs[0] * taus_s[0] * taus_s[0];
            coeffs->h_2 = h_2_est * 3.0;
        }
    }

    /* White FM: tau^-1/2 slope region */
    for (int i = 0; i < num_points - 1; i++) {
        if (taus_s[i] < 1e-12) continue;
        double slope = log10(adevs[i+1] / adevs[i])
                     / log10(taus_s[i+1] / taus_s[i]);
        if (fabs(slope + 0.5) < 0.2) {
            /* h_0: White FM, sigma^2 = h_0/(2*tau)
             * h_0 = 2*tau*sigma^2
             */
            coeffs->h_0 = 2.0 * taus_s[i] * adevs[i] * adevs[i];
            break;
        }
    }

    /* Flicker FM: sigma ~ tau^0 (flat) */
    for (int i = 0; i < num_points - 1; i++) {
        double slope = log10(adevs[i+1] / adevs[i])
                     / log10(taus_s[i+1] / taus_s[i]);
        if (fabs(slope) < 0.15) {
            /* h_-1: Flicker FM, sigma^2 = h_-1 * 2*ln(2)
             * h_-1 = sigma^2 / (2*ln(2))
             */
            coeffs->h_m1 = adevs[i] * adevs[i] / (2.0 * M_LN2);
            break;
        }
    }

    /* Random Walk FM: tau^+1/2 slope region */
    for (int i = 0; i < num_points - 1; i++) {
        double slope = log10(adevs[i+1] / adevs[i])
                     / log10(taus_s[i+1] / taus_s[i]);
        if (fabs(slope - 0.5) < 0.2) {
            /* h_-2: Random Walk FM, sigma^2 = (2*pi^2/3)*h_-2*tau
             * h_-2 = 3*sigma^2 / (2*pi^2*tau)
             */
            coeffs->h_m2 = 3.0 * adevs[i] * adevs[i]
                         / (2.0 * M_PI * M_PI * taus_s[i]);
            break;
        }
    }

    /* Flicker PM (difficult to separate from White PM with ADEV alone,
     * use modified Allan for better discrimination) */
    coeffs->h_1 = 0.0; /* Requires MVAR for good estimate */

    return 0;
}

/* ===================================================================
 * L5: Noise Type Identification
 * =================================================================== */

NoiseType allan_identify_noise_type(double adev1, double tau1_s,
                                    double adev2, double tau2_s)
{
    if (tau1_s <= 0.0 || tau2_s <= 0.0 || adev1 <= 0.0 || adev2 <= 0.0) {
        return NOISE_WHITE_FM;
    }
    if (fabs(tau2_s - tau1_s) < 1e-30) {
        return NOISE_WHITE_FM;
    }

    /* Slope in log-log:
     *   mu = log10(adev2/adev1) / log10(tau2/tau1)
     *
     * Noise type mapping:
     *   mu ~ -1.0  -> White PM or Flicker PM
     *   mu ~ -0.5  -> White FM
     *   mu ~  0.0  -> Flicker FM
     *   mu ~ +0.5  -> Random Walk FM
     *   mu ~ +1.0  -> Frequency Drift
     */
    double mu = log10(adev2 / adev1) / log10(tau2_s / tau1_s);

    if (mu < -0.75) return NOISE_WHITE_PM;
    if (mu < -0.25) return NOISE_WHITE_FM;
    if (mu <  0.25) return NOISE_FLICKER_FM;
    if (mu <  0.75) return NOISE_RANDOM_WALK_FM;
    return NOISE_FREQ_DRIFT;
}

/* ===================================================================
 * L6: Time Deviation (TDEV)
 * =================================================================== */

int time_deviation_compute(const double *xi, int N, double tau0_s,
                           double tau_s, double *tdev)
{
    if (!xi || !tdev || N < 4) return -1;

    /* TDEV from phase data:
     *   TDEV(tau) = sqrt( 1/(6*(N-3m+1)) *
     *     sum_{i=1}^{N-3m+1} [ sum_{j=i}^{i+m-1} (delta^2 x_j)^2 ] )
     *
     * where delta^2 is the second difference operator:
     *   delta^2 x_j = x_{j+2m} - 2*x_{j+m} + x_j
     */

    int m = (int)(tau_s / tau0_s + 0.5);
    if (m < 1) m = 1;
    if (m > N / 3) m = N / 3;

    int max_i = N - 3 * m + 1;
    if (max_i < 1) return -1;

    double sum_sq = 0.0;
    for (int i = 0; i < max_i; i++) {
        double inner_sum = 0.0;
        for (int j = i; j < i + m; j++) {
            /* Second difference: delta^2 x_j */
            double d2 = xi[j + 2*m] - 2.0 * xi[j + m] + xi[j];
            inner_sum += d2 * d2;
        }
        sum_sq += inner_sum;
    }

    *tdev = sqrt(sum_sq / (6.0 * (double)m * (double)m * (double)max_i));
    return 0;
}

/* ===================================================================
 * L6: Hadamard Variance
 * =================================================================== */

int hadamard_variance(const double *yi, int N, double tau0_s,
                      double tau_s, double *hadamard)
{
    if (!yi || !hadamard || N < 5) return -1;

    /* Hadamard variance:
     *   H_sigma_y^2(tau) = 1/(6*(M-2)) *
     *     sum_{i=1}^{M-2} (y_{i+2} - 2*y_{i+1} + y_i)^2
     *
     * where each y_i is an average over tau (m samples).
     * Insensitive to linear frequency drift.
     */

    int m = (int)(tau_s / tau0_s + 0.5);
    if (m < 1) m = 1;
    if (m > N / 3) m = N / 3;

    int M = N / m; /* Number of averaged samples */
    if (M < 3) return -1;

    /* Compute averaged frequency values */
    double *y_avg = (double *)malloc((size_t)M * sizeof(double));
    if (!y_avg) return -1;

    for (int i = 0; i < M; i++) {
        double sum = 0.0;
        for (int k = 0; k < m; k++) {
            sum += yi[i * m + k];
        }
        y_avg[i] = sum / (double)m;
    }

    double sum_sq = 0.0;
    for (int i = 0; i < M - 2; i++) {
        /* Second difference of averaged frequencies */
        double d2 = y_avg[i+2] - 2.0 * y_avg[i+1] + y_avg[i];
        sum_sq += d2 * d2;
    }

    *hadamard = sqrt(sum_sq / (6.0 * (double)(M - 2)));

    free(y_avg);
    return 0;
}

/* ===================================================================
 * L7: GPSDO Stability Check
 * =================================================================== */

int gpsdo_stability_check(double adev_1s, double adev_100s,
                          double adev_10000s, char grade)
{
    /* GPSDO (GPS-Disciplined Oscillator) performance grades:
     *
     * Standard Grade (S):
     *   ADEV(1s)   < 5e-11
     *   ADEV(100s)  < 5e-12
     *   ADEV(10000s) < 1e-12
     *
     * Precision Grade (P):
     *   ADEV(1s)   < 1e-11
     *   ADEV(100s)  < 2e-12
     *   ADEV(10000s) < 5e-13
     *
     * Ultra-Stable Grade (U):
     *   ADEV(1s)   < 5e-12
     *   ADEV(100s)  < 5e-13
     *   ADEV(10000s) < 1e-13
     *
     * Reference: typical GPSDO datasheets (Trimble, Symmetricom/Microsemi)
     */

    switch (grade) {
    case 'S':
    case 's':
        return (adev_1s <= 5.0e-11 && adev_100s <= 5.0e-12
                && adev_10000s <= 1.0e-12) ? 1 : 0;

    case 'P':
    case 'p':
        return (adev_1s <= 1.0e-11 && adev_100s <= 2.0e-12
                && adev_10000s <= 5.0e-13) ? 1 : 0;

    case 'U':
    case 'u':
        return (adev_1s <= 5.0e-12 && adev_100s <= 5.0e-13
                && adev_10000s <= 1.0e-13) ? 1 : 0;

    default:
        return 0;
    }
}

/* ===================================================================
 * L8: Three-Cornered Hat Method
 * =================================================================== */

int three_cornered_hat(double sigma_AB2, double sigma_BC2, double sigma_CA2,
                       double *sigma_A2, double *sigma_B2, double *sigma_C2)
{
    if (!sigma_A2 || !sigma_B2 || !sigma_C2) return -1;

    /* Three-cornered hat:
     *
     * For three oscillators A, B, C with pairwise measurements:
     *   sigma_AB^2 = sigma_A^2 + sigma_B^2
     *   sigma_BC^2 = sigma_B^2 + sigma_C^2
     *   sigma_CA^2 = sigma_C^2 + sigma_A^2
     *
     * Solution:
     *   sigma_A^2 = (sigma_AB^2 + sigma_CA^2 - sigma_BC^2) / 2
     *   sigma_B^2 = (sigma_AB^2 + sigma_BC^2 - sigma_CA^2) / 2
     *   sigma_C^2 = (sigma_BC^2 + sigma_CA^2 - sigma_AB^2) / 2
     *
     * Validity check: each result must be >= 0.
     * Negative values indicate measurement inconsistency
     * (triangle inequality violated).
     */

    double A2 = (sigma_AB2 + sigma_CA2 - sigma_BC2) / 2.0;
    double B2 = (sigma_AB2 + sigma_BC2 - sigma_CA2) / 2.0;
    double C2 = (sigma_BC2 + sigma_CA2 - sigma_AB2) / 2.0;

    /* Triangle inequality check */
    if (A2 < -1.0e-30 || B2 < -1.0e-30 || C2 < -1.0e-30) {
        /* Clamp small negatives to zero (measurement noise) */
        if (A2 < -1.0e-6 || B2 < -1.0e-6 || C2 < -1.0e-6) {
            return -1; /* Significant violation */
        }
    }

    *sigma_A2 = (A2 < 0.0) ? 0.0 : A2;
    *sigma_B2 = (B2 < 0.0) ? 0.0 : B2;
    *sigma_C2 = (C2 < 0.0) ? 0.0 : C2;

    return 0;
}
