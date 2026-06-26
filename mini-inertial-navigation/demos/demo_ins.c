/**
 * @file    demo_ins.c
 * @brief   Demo: Visualize INS error accumulation and Schuler oscillation
 *
 * Demonstrates the characteristic 84.4-minute Schuler oscillation
 * of INS errors when un-aided by external measurements.
 *
 * This demo shows why INS cannot operate standalone for extended periods:
 * even navigation-grade INS drifts at ~1 nautical mile per hour.
 */

#include "ins_core.h"
#include "ins_errors.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== INS Error Accumulation Demo ===\n\n");

    /* Compare different IMU grades over time */
    ins_grade_t grades[] = {
        INS_GRADE_CONSUMER,
        INS_GRADE_INDUSTRIAL,
        INS_GRADE_TACTICAL,
        INS_GRADE_NAVIGATION
    };
    const char *names[] = {"Consumer", "Industrial", "Tactical", "Navigation-grade"};
    double times[] = {10.0, 60.0, 300.0, 1800.0, 3600.0, 7200.0};
    const char *time_labels[] = {"10 s", "1 min", "5 min", "30 min", "1 hr", "2 hr"};

    printf("Horizontal Position Drift [m] vs Time:\n\n");
    printf("%-18s", "IMU Grade");
    for (int t = 0; t < 6; t++) printf("%10s", time_labels[t]);
    printf("\n");

    for (int g = 0; g < 4; g++) {
        const ins_grade_spec_t *spec = ins_grade_spec_get(grades[g]);
        printf("%-18s", names[g]);

        for (int t = 0; t < 6; t++) {
            /* Build error model from grade specs */
            ins_imu_error_model_t model;
            memset(&model, 0, sizeof(model));

            double gyro_bias_rads = spec->gyro_bias * (M_PI / 180.0) / 3600.0;
            model.gyro_x.bias_offset = gyro_bias_rads;
            model.gyro_y.bias_offset = gyro_bias_rads;
            model.gyro_z.bias_offset = gyro_bias_rads;
            model.gyro_x.white_noise_std = spec->gyro_arw * (M_PI / 180.0) / 3600.0;
            model.gyro_y.white_noise_std = spec->gyro_arw * (M_PI / 180.0) / 3600.0;
            model.gyro_z.white_noise_std = spec->gyro_arw * (M_PI / 180.0) / 3600.0;

            double pos_drift, vel_drift, att_drift;
            ins_error_predict_drift(&model, times[t], &pos_drift, &vel_drift, &att_drift);

            if (pos_drift < 1000.0) {
                printf("%10.1f", pos_drift);
            } else if (pos_drift < 1e6) {
                printf("%9.1fk", pos_drift / 1000.0);
            } else {
                printf("%9.1fM", pos_drift / 1e6);
            }
        }
        printf("\n");
    }

    printf("\nKey observations:\n");
    printf("  1. Consumer IMU: unusable after seconds (phone orientation only)\n");
    printf("  2. Industrial IMU: ~1 km drift after 1 hour (needs frequent GPS)\n");
    printf("  3. Tactical IMU: ~10 m drift after 1 hour (missile-grade)\n");
    printf("  4. Navigation-grade: ~1 nm/hr drift (aviation standard)\n");
    printf("\nSchuler oscillation period: %.0f seconds (~84.4 minutes)\n", INS_SCHULER_PERIOD);
    printf("This is the natural period of INS error dynamics in the horizontal channel.\n");
    return 0;
}
