/* =========================================================================
 * example_carrier_smoothing.c — Carrier-smoothed pseudorange and cycle slips
 *
 * Demonstrates: L2 (carrier phase concept), L5 (Hatch smoothing),
 * L6 (cycle slip detection and recovery), L7 (precise navigation).
 * ========================================================================= */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "gnss_common.h"
#include "gnss_carrier.h"

int main(void) {
    printf("=== Carrier-Smoothed Pseudorange Example ===\n\n");

    int n_epochs = 100;
    double *true_ranges = (double*)malloc(n_epochs * sizeof(double));
    double *code_meas   = (double*)malloc(n_epochs * sizeof(double));
    double *phase_meas  = (double*)malloc(n_epochs * sizeof(double));
    double *hatch_smoothed = (double*)malloc(n_epochs * sizeof(double));

    unsigned int seed = 12345;
    double iono_rate = 0.01;

    int i;
    for (i = 0; i < n_epochs; i++) {
        true_ranges[i] = 20000000.0 + i * 500.0;

        double u1 = (double)((seed = seed * 1103515245 + 12345) % 100000) / 100000.0;
        double u2 = (double)((seed = seed * 1103515245 + 12345) % 100000) / 100000.0;
        double n1 = sqrt(-2.0 * log(u1 + 1e-10)) * cos(2.0 * M_PI * u2);
        double n2 = sqrt(-2.0 * log(u1 + 1e-10)) * sin(2.0 * M_PI * u2);

        double iono_delay = iono_rate * i;
        code_meas[i] = true_ranges[i] + iono_delay + 3.0 * n1;
        phase_meas[i] = (true_ranges[i] - iono_delay) / GNSS_L1_WAVELENGTH + 0.005 * n2;

        if (i >= 60) phase_meas[i] += 5.0;
    }

    gnss_hatch_filter_t hatch;
    gnss_hatch_init(&hatch, 50, 1);

    printf("Epoch | Code(m) | Phase(m) | Smoothed(m) | Error(m) | Slip?\n");
    printf("------|---------|----------|-------------|----------|-------\n");

    double rms_raw = 0.0, rms_smooth = 0.0;
    int slip_count = 0;

    for (i = 0; i < n_epochs; i++) {
        double phase_m = phase_meas[i] * GNSS_L1_WAVELENGTH;
        double iono_delay = iono_rate * i;

        int slip = 0;
        if (i > 0) {
            double delta_phi = phase_meas[i] - phase_meas[i-1];
            double expected = -2626.0;
            slip = (fabs(delta_phi - expected) > 3.0) ? 1 : 0;
            if (slip) slip_count++;
        }

        hatch_smoothed[i] = gnss_hatch_smooth(&hatch, code_meas[i], phase_m, slip);

        double raw_err = code_meas[i] - (true_ranges[i] + iono_delay);
        double smooth_err = hatch_smoothed[i] - (true_ranges[i] + iono_delay);

        rms_raw += raw_err * raw_err;
        rms_smooth += smooth_err * smooth_err;

        if (i < 10 || i == 59 || i == 60 || i == 61 || i >= 90) {
            printf("%5d | %7.2f | %8.2f | %11.2f | %8.3f | %5s\n",
                   i, code_meas[i], phase_m, hatch_smoothed[i],
                   smooth_err, slip ? "YES" : "no");
        }
    }

    rms_raw = sqrt(rms_raw / n_epochs);
    rms_smooth = sqrt(rms_smooth / n_epochs);

    printf("\n=== Results ===\n");
    printf("Raw code RMS error:          %.2f m\n", rms_raw);
    printf("Carrier-smoothed RMS error:  %.2f m\n", rms_smooth);
    printf("Noise reduction factor:      %.1fx\n", rms_raw / (rms_smooth + 0.001));
    printf("Cycle slips detected:        %d\n", slip_count);

    printf("\n=== Ionosphere-Free Combination ===\n");
    gnss_ionofree_combo_t ic = gnss_ionofree_combination(GNSS_L1_FREQ, GNSS_L2_FREQ);
    printf("GPS L1/L2 IF combination:\n");
    printf("  P1 coefficient: %.6f\n", ic.coeff_p1);
    printf("  P2 coefficient: %.6f\n", ic.coeff_p2);
    printf("  Noise factor:   %.2f\n", ic.noise_factor);

    double P1 = 20000000.0 + 10.0;
    double P2 = 20000000.0 + 10.0 * pow(GNSS_L1_FREQ/GNSS_L2_FREQ, 2.0);
    double P_IF = gnss_ionofree_code(P1, P2, &ic);
    printf("  P1=%.3f m (with 10m iono), P2=%.3f m, P_IF=%.3f m\n", P1, P2, P_IF);
    printf("  Ionospheric residual: %.4f m\n", P_IF - 20000000.0);

    free(true_ranges);
    free(code_meas);
    free(phase_meas);
    free(hatch_smoothed);

    printf("\n=== Example Complete ===\n");
    return 0;
}
