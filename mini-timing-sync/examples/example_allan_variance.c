/**
 * @file example_allan_variance.c
 * @brief Allan variance analysis example
 *
 * Demonstrates: Computing Allan deviation for simulated oscillator data,
 * noise type identification, multi-tau sigma-tau plot generation,
 * GPSDO compliance check, three-cornered hat decomposition.
 *
 * L6 - Canonical Problem: Oscillator stability characterization
 * L7 - Application: GPSDO (GPS-Disciplined Oscillator) analysis
 * L8 - Advanced: Three-cornered hat method
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "allan_variance.h"

/* Generate simulated frequency data with mixed noise types */
static void generate_test_data(double *yi, int N, double tau0)
{
    /* White FM noise dominates at short tau */
    double white_fm_level = 1e-10;
    /* Flicker FM adds at mid tau */
    double flicker_fm_level = 5e-12;
    /* Random walk FM dominates at long tau */
    double rw_fm_level = 1e-13;

    double rw_state = 0.0;

    for (int i = 0; i < N; i++) {
        /* White FM */
        double wfm = white_fm_level * ((double)rand() / RAND_MAX - 0.5);

        /* Flicker FM (simplified as filtered white noise) */
        double ffm = flicker_fm_level * ((double)rand() / RAND_MAX - 0.5);

        /* Random walk FM */
        rw_state += rw_fm_level * ((double)rand() / RAND_MAX - 0.5);

        yi[i] = wfm + ffm + rw_state;
    }
}

int main(void)
{
    printf("=== Allan Variance Analysis Example ===\n\n");

    /* Generate 100000 samples at 1 Hz (tau0 = 1s) */
    int N = 100000;
    double tau0_s = 1.0;
    double *yi = (double*)malloc((size_t)N * sizeof(double));
    if (!yi) { printf("Memory allocation failed\n"); return 1; }

    printf("Generating %d simulated frequency measurements...\n", N);
    generate_test_data(yi, N, tau0_s);

    /* Compute Allan deviation at multiple tau values */
    double taus_s[] = {1.0, 10.0, 100.0, 1000.0, 10000.0};
    int num_taus = 5;
    AllanResult results[5];

    printf("\nComputing Allan deviation at multiple tau values:\n");
    printf("  tau(s)    | ADEV       | Noise Type\n");
    printf("  ----------|------------|------------\n");

    allan_deviation_multi_tau(yi, N, tau0_s, taus_s, num_taus, results);

    for (int i = 0; i < num_taus; i++) {
        printf("  %9.1f | %10.3e",
               results[i].tau_s, results[i].adev);

        /* Identify noise type between consecutive tau values */
        if (i > 0) {
            NoiseType nt = allan_identify_noise_type(
                results[i-1].adev, results[i-1].tau_s,
                results[i].adev, results[i].tau_s);
            const char *nt_str;
            switch (nt) {
            case NOISE_WHITE_PM:       nt_str = "White PM"; break;
            case NOISE_FLICKER_PM:     nt_str = "Flicker PM"; break;
            case NOISE_WHITE_FM:       nt_str = "White FM"; break;
            case NOISE_FLICKER_FM:     nt_str = "Flicker FM"; break;
            case NOISE_RANDOM_WALK_FM: nt_str = "RW FM"; break;
            case NOISE_FREQ_DRIFT:     nt_str = "Freq Drift"; break;
            default:                   nt_str = "Unknown"; break;
            }
            printf(" | %s", nt_str);
        }
        printf("\n");
    }

    /* Fit noise coefficients */
    printf("\nFitting noise coefficients...\n");
    NoiseCoefficients coeffs;
    double adevs[5];
    for (int i = 0; i < num_taus; i++) adevs[i] = results[i].adev;
    allan_fit_noise_coefficients(taus_s, adevs, num_taus, &coeffs);

    printf("  h_2  (White PM):     %.3e\n", coeffs.h_2);
    printf("  h_1  (Flicker PM):   %.3e\n", coeffs.h_1);
    printf("  h_0  (White FM):     %.3e\n", coeffs.h_0);
    printf("  h_-1 (Flicker FM):   %.3e\n", coeffs.h_m1);
    printf("  h_-2 (RW FM):        %.3e\n", coeffs.h_m2);

    /* Hadamard variance at tau=10000s (drift-insensitive) */
    double hadamard;
    if (hadamard_variance(yi, N, tau0_s, 10000.0, &hadamard) == 0) {
        printf("\nHadamard deviation at 10000s: %.3e\n", hadamard);
        printf("  (Insensitive to linear frequency drift)\n");
    }

    /* GPSDO grade check */
    printf("\nGPSDO Stability Grade Check:\n");
    double adev_1s = results[0].adev;
    double adev_100s = results[2].adev;
    double adev_10000s = results[4].adev;

    printf("  ADEV(1s)     = %.3e\n", adev_1s);
    printf("  ADEV(100s)   = %.3e\n", adev_100s);
    printf("  ADEV(10000s) = %.3e\n", adev_10000s);

    char grades[] = {'S', 'P', 'U'};
    const char *grade_names[] = {"Standard", "Precision", "Ultra-Stable"};
    for (int g = 0; g < 3; g++) {
        int ok = gpsdo_stability_check(adev_1s, adev_100s, adev_10000s,
                                        grades[g]);
        printf("  %s Grade: %s\n", grade_names[g], ok ? "PASS" : "FAIL");
    }

    /* Three-cornered hat example */
    printf("\nThree-Cornered Hat Example:\n");
    double sA2, sB2, sC2;
    if (three_cornered_hat(4.0e-20, 9.0e-20, 11.0e-20, &sA2, &sB2, &sC2) == 0) {
        printf("  Oscillator A ADEV: %.3e\n", sqrt(sA2));
        printf("  Oscillator B ADEV: %.3e\n", sqrt(sB2));
        printf("  Oscillator C ADEV: %.3e\n", sqrt(sC2));
    }

    free(yi);
    printf("\n=== Example Complete ===\n");
    return 0;
}
