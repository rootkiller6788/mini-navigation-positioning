/**
 * @file    bench_ins.c
 * @brief   Performance benchmark for INS core computations
 *
 * Benchmarks key INS operations:
 *   - Quaternion update throughput
 *   - Mechanization step rate
 *   - Allan variance computation time
 */

#include "ins_core.h"
#include "ins_attitude.h"
#include "ins_mechanization.h"
#include "ins_errors.h"
#include <stdio.h>
#include <time.h>

#define BENCH_ITERATIONS 1000000

static double get_time_ms(void) {
    return (double)clock() / CLOCKS_PER_SEC * 1000.0;
}

int main(void) {
    printf("=== INS Performance Benchmarks ===\n\n");

    /* Benchmark 1: Quaternion update */
    {
        ins_quat_t q;
        ins_quat_identity(&q);
        ins_vec3_t omega = {0.01, 0.02, 0.03};
        double t0 = get_time_ms();
        for (int i = 0; i < BENCH_ITERATIONS; i++) {
            ins_quat_update_exact(&q, &omega, 0.01);
        }
        double t1 = get_time_ms();
        printf("Quaternion update:       %10.3f ns/op (%d iterations in %.1f ms)\n",
               (t1 - t0) / BENCH_ITERATIONS * 1e6, BENCH_ITERATIONS, t1 - t0);
    }

    /* Benchmark 2: Vector cross product */
    {
        ins_vec3_t a = {1, 2, 3}, b = {4, 5, 6}, c;
        double t0 = get_time_ms();
        for (int i = 0; i < BENCH_ITERATIONS; i++) {
            ins_vec3_cross(&a, &b, &c);
        }
        double t1 = get_time_ms();
        printf("Vector cross product:    %10.3f ns/op\n",
               (t1 - t0) / BENCH_ITERATIONS * 1e6);
    }

    /* Benchmark 3: Mechanization step */
    {
        ins_mech_state_t ms;
        ins_mech_init(&ms, INS_MECH_EXACT, 0.5, -2.0, 100.0);
        ins_imu_sample_t imu;
        imu.accel.x = 0; imu.accel.y = 0; imu.accel.z = 9.81;
        imu.gyro.x = 0; imu.gyro.y = 0; imu.gyro.z = 0.001;
        imu.dt = 0.01;

        int iter = BENCH_ITERATIONS / 10;
        double t0 = get_time_ms();
        for (int i = 0; i < iter; i++) {
            ins_mech_step(&ms, &imu);
        }
        double t1 = get_time_ms();
        printf("Mechanization step:      %10.3f us/op (%d iterations in %.1f ms)\n",
               (t1 - t0) / iter * 1e3, iter, t1 - t0);
    }

    /* Benchmark 4: Allan Variance */
    {
        double data[2048];
        for (int i = 0; i < 2048; i++) {
            data[i] = (double)(rand() % 1000) / 100000.0;
        }
        double taus[128], adevs[128];
        double t0 = get_time_ms();
        ins_allan_variance(data, 2048, 0.01, taus, adevs, 128);
        double t1 = get_time_ms();
        printf("Allan variance (N=2048): %10.3f ms total\n", t1 - t0);
    }

    printf("\nAll benchmarks complete.\n");
    return 0;
}
