#include <stdio.h>
#include <math.h>
#include "../include/uwb_types.h"
#include "../include/uwb_positioning.h"

int main(void) {
    /* 5 anchors in an office room (10m x 10m), tag walks through */
    uwb_pos3d_t anchors[5] = {
        {0.0, 0.0, 2.5},     /* Corner ceiling */
        {10.0, 0.0, 2.5},    /* Corner ceiling */
        {10.0, 10.0, 2.5},   /* Corner ceiling */
        {0.0, 10.0, 2.5},    /* Corner ceiling */
        {5.0, 5.0, 0.5}      /* Center floor */
    };

    /* Tag positions along a path */
    uwb_pos3d_t path[] = {
        {1.0, 1.0, 1.0},
        {3.0, 2.0, 1.0},
        {5.0, 5.0, 1.0},
        {7.0, 8.0, 1.0},
        {9.0, 9.0, 1.0}
    };
    int n_path = 5;
    int i, j;

    printf("Multilateration (Linear LS) tracking demo:\n");
    printf("Room: 10m x 10m, 5 anchors at corners + center\n\n");

    for (i = 0; i < n_path; i++) {
        double ranges[5];
        for (j = 0; j < 5; j++) {
            double dx = path[i].x - anchors[j].x;
            double dy = path[i].y - anchors[j].y;
            double dz = path[i].z - anchors[j].z;
            ranges[j] = sqrt(dx*dx + dy*dy + dz*dz);
            /* Add small Gaussian noise */
            ranges[j] += ((double)((j*7+13*i)%100)/100.0 - 0.5) * 0.1;
        }

        uwb_positioning_result_t result;
        result.position.x = result.position.y = result.position.z = 0.0;

        if (multilateration_linear_ls(anchors, ranges, 5, POS_DIM_3D, &result)) {
            double err = sqrt(
                (result.position.x-path[i].x)*(result.position.x-path[i].x) +
                (result.position.y-path[i].y)*(result.position.y-path[i].y) +
                (result.position.z-path[i].z)*(result.position.z-path[i].z));
            printf("Step %d: est=(%.3f, %.3f, %.3f) true=(%.1f, %.1f, %.1f) "
                   "err=%.3fm GDOP=%.2f\n",
                   i, result.position.x, result.position.y, result.position.z,
                   path[i].x, path[i].y, path[i].z, err, result.gdop);
        }
    }

    return 0;
}
