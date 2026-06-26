/* Mini bench: measure key GNSS operations throughput */
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "gnss_common.h"
#include "gnss_signal.h"
#include "gnss_position.h"

int main(void) {
    clock_t start, end;
    int i, n_iter = 10000;

    /* 1. C/A code generation */
    gnss_ca_code_t code;
    start = clock();
    for (i = 0; i < n_iter; i++) gnss_ca_code_generate((i%32)+1, &code);
    end = clock();
    printf("C/A code gen x%d: %.2f ms\n", n_iter,
           1000.0*(double)(end-start)/CLOCKS_PER_SEC);

    /* 2. Kepler solver */
    start = clock();
    for (i = 0; i < n_iter; i++) gnss_kepler_solve(0.5, 0.01, 10, 1e-12);
    end = clock();
    printf("Kepler solve x%d: %.2f ms\n", n_iter,
           1000.0*(double)(end-start)/CLOCKS_PER_SEC);

    /* 3. ECEF→LLA conversion */
    gnss_ecef_t pos = {4200000, -170000, 4780000};
    start = clock();
    for (i = 0; i < n_iter; i++) gnss_ecef_to_lla(pos);
    end = clock();
    printf("ECEF→LLA x%d: %.2f ms\n", n_iter,
           1000.0*(double)(end-start)/CLOCKS_PER_SEC);

    printf("Bench complete.\n");
    return 0;
}
