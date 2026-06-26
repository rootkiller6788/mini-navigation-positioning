/**
 * @file example_ptp_sync.c
 * @brief End-to-end PTP slave synchronization example
 *
 * Demonstrates: PTP two-way timestamp exchange, slave servo loop,
 * Kalman clock tracking, holdover entry/exit, and
 * sync accuracy monitoring.
 *
 * L6 - Canonical Problem: PTP Master-Slave Synchronization
 * L7 - Application: 5G Fronthaul timing compliance check
 *
 * Run: make examples && ./example_ptp_sync
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "timing_sync.h"
#include "ptp_engine.h"
#include "clock_model.h"

int main(void)
{
    printf("=== PTP Slave Synchronization Example ===\n\n");

    /* Initialize PTP slave */
    PtpSlaveState slave;
    ptp_slave_init(&slave, 500.0, 1.0);  /* Initial offset: 500 ns */

    printf("Initial state:\n");
    printf("  Status: FREE_RUNNING\n");
    printf("  Initial offset estimate: 500 ns\n\n");

    /* Simulate 20 PTP sync exchanges */
    printf("Simulating 20 PTP sync exchanges...\n");
    printf("  Step |  Offset(ns) | Delay(ns) |  Status      | Correction(ns)\n");
    printf("  -----|-------------|-----------|--------------|---------------\n");

    double true_offset = 500.0;  /* True offset starts at 500 ns */
    double true_delay = 1000.0;  /* True path delay */

    for (int step = 0; step < 20; step++) {
        /* Simulate timestamps with true offset and delay plus noise */
        PtpTimestamps ts;

        /* Master sends sync at master time */
        ts.t1.seconds = 1000 + step;
        ts.t1.nanoseconds = 0;

        /* Slave receives: t2 = t1 + delay + offset + noise */
        double noise_t2 = (double)(rand() % 11 - 5); /* +/-5 ns noise */
        ts.t2 = ts.t1;
        timing_timestamp_add_ns(&ts.t2,
            (int64_t)(true_delay + true_offset + noise_t2));

        /* Slave sends delay_req after processing */
        ts.t3 = ts.t2;
        timing_timestamp_add_ns(&ts.t3, 500000); /* 0.5 ms later */

        /* Master receives: t4 = t3 + delay - offset + noise */
        double noise_t4 = (double)(rand() % 11 - 5);
        ts.t4 = ts.t3;
        timing_timestamp_add_ns(&ts.t4,
            (int64_t)(true_delay - true_offset + noise_t4));

        /* Process and get servo correction */
        double correction = ptp_slave_update(&slave, &ts);

        /* Apply correction to true offset (simulate clock adjustment) */
        true_offset -= correction;

        const char *status_str;
        switch (slave.sync_status) {
        case SYNC_FREE_RUNNING: status_str = "FREE_RUNNING"; break;
        case SYNC_ACQUIRING:    status_str = "ACQUIRING   "; break;
        case SYNC_LOCKED:       status_str = "LOCKED      "; break;
        case SYNC_HOLDOVER:     status_str = "HOLDOVER    "; break;
        default:                status_str = "LOSS        "; break;
        }

        printf("  %4d | %10.1f | %9.1f | %-12s | %13.1f\n",
               step + 1,
               slave.offset_from_master_ns,
               slave.mean_path_delay_ns,
               status_str,
               correction);

        /* Check 5G Class B compliance */
        if (slave.sync_status == SYNC_LOCKED) {
            int class_b_ok = ptp_5g_fronthaul_check(
                slave.offset_from_master_ns, 'B');
            printf("         -> 5G Class B compliant: %s\n",
                   class_b_ok ? "YES" : "NO");
        }
    }

    /* Final state */
    printf("\nFinal synchronization state:\n");
    printf("  Status:           ");
    switch (slave.sync_status) {
    case SYNC_LOCKED: printf("LOCKED\n"); break;
    default:          printf("NOT LOCKED\n"); break;
    }
    printf("  Remaining offset: %.1f ns\n", slave.offset_from_master_ns);
    printf("  Path delay:       %.1f ns\n", slave.mean_path_delay_ns);

    double accuracy = ptp_sync_accuracy(&slave);
    printf("  Sync accuracy:    %.1f ns (1-sigma)\n", accuracy);
    printf("  95%% confidence:    %.1f ns\n", accuracy * 2.0);

    /* Power grid compliance */
    int grid_ok = ptp_power_grid_check(slave.offset_from_master_ns);
    printf("  IEC 61850 power grid compliant: %s\n",
           grid_ok ? "YES" : "NO");

    printf("\n=== Example Complete ===\n");
    return 0;
}
