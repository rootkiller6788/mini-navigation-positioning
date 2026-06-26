#include <stdio.h>
#include <math.h>
#include "../include/uwb_types.h"
#include "../include/uwb_positioning.h"

int main(void) {
    /* 2D trilateration: 3 anchors forming a triangle, tag at center */
    uwb_pos3d_t anchors[3] = {
        {0.0, 0.0, 0.0},
        {10.0, 0.0, 0.0},
        {5.0, 8.660254, 0.0}  /* 8.66 = 10*sin(60deg), equilateral triangle */
    };

    /* Tag at (5.0, 2.887) - centroid */
    double tag_x = 5.0, tag_y = 2.886751;
    double ranges[3];
    int i;
    for (i = 0; i < 3; i++) {
        double dx = tag_x - anchors[i].x;
        double dy = tag_y - anchors[i].y;
        ranges[i] = sqrt(dx*dx + dy*dy);
    }

    uwb_pos2d_t result;
    if (trilateration_2d(anchors, ranges, &result)) {
        printf("Trilateration 2D Result:\n");
        printf("  Estimated: (%.4f, %.4f)\n", result.x, result.y);
        printf("  True:      (%.4f, %.4f)\n", tag_x, tag_y);
        printf("  Error:     %.4f m\n",
               sqrt((result.x-tag_x)*(result.x-tag_x) + (result.y-tag_y)*(result.y-tag_y)));
    } else {
        printf("Trilateration failed!\n");
    }

    /* 3D trilateration: 4 anchors */
    uwb_pos3d_t anchors3d[4] = {
        {0.0, 0.0, 0.0},
        {10.0, 0.0, 0.0},
        {0.0, 10.0, 0.0},
        {0.0, 0.0, 10.0}
    };
    double tag3d_x = 3.0, tag3d_y = 4.0, tag3d_z = 5.0;
    double ranges3d[4];
    for (i = 0; i < 4; i++) {
        double dx = tag3d_x - anchors3d[i].x;
        double dy = tag3d_y - anchors3d[i].y;
        double dz = tag3d_z - anchors3d[i].z;
        ranges3d[i] = sqrt(dx*dx + dy*dy + dz*dz);
    }

    uwb_pos3d_t result3d;
    if (trilateration_3d(anchors3d, ranges3d, &result3d)) {
        printf("\nTrilateration 3D Result:\n");
        printf("  Estimated: (%.4f, %.4f, %.4f)\n", result3d.x, result3d.y, result3d.z);
        printf("  True:      (%.4f, %.4f, %.4f)\n", tag3d_x, tag3d_y, tag3d_z);
        double err = sqrt((result3d.x-tag3d_x)*(result3d.x-tag3d_x) +
                          (result3d.y-tag3d_y)*(result3d.y-tag3d_y) +
                          (result3d.z-tag3d_z)*(result3d.z-tag3d_z));
        printf("  Error:     %.6f m\n", err);
    }

    return 0;
}
