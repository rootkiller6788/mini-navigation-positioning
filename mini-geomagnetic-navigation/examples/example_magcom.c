/**
 * example_magcom.c -- Magnetic Contour Matching demo
 */
#include "geomag_navigation.h"
#include "geomag_math.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void) {
    printf("=== MAGCOM Demo ===\n");
    MagneticMap map;
    map.origin.lat = 40.0; map.origin.lon = -105.0; map.origin.alt = 0.0;
    map.grid_spacing_lat = 0.01; map.grid_spacing_lon = 0.01;
    map.nlat = 20; map.nlon = 20;
    int np = map.nlat * map.nlon;
    map.total_field = (double *)malloc(np * sizeof(double));
    map.inclination = NULL; map.declination = NULL; map.anomaly = NULL;
    for (int i = 0; i < map.nlat; i++)
        for (int j = 0; j < map.nlon; j++) {
            double lat = map.origin.lat + i * map.grid_spacing_lat;
            double lon = map.origin.lon + j * map.grid_spacing_lon;
            map.total_field[i*map.nlon+j] = 52000.0 + lat*100.0 + lon*50.0;
        }
    int N = 15;
    double *meas = malloc(N*sizeof(double));
    double *tl = malloc(N*sizeof(double));
    double *to = malloc(N*sizeof(double));
    GeodeticCoord start = {40.05, -104.95, 0.0};
    generate_magnetic_profile(&map, &start, 45.0, 100.0, N, meas, tl, to);
    printf("Profile: %d points.\n", N);
    for (int i = 0; i < 3; i++)
        printf("  [%d] (%.4f,%.4f): %.1f nT\n", i, tl[i], to[i], meas[i]);
    double elat, elon, conf;
    magcom_correlation_match(&map, meas, tl, to, N, 0.1, 0.005, &elat, &elon, &conf);
    printf("NCC offset=(%.4f,%.4f) conf=%.3f\n", elat, elon, conf);
    double roughness;
    GeodeticCoord ctr = {40.05, -104.95, 0.0};
    magmap_roughness_metric(&map, &ctr, 0.05, &roughness);
    printf("Roughness: %.1f nT\n", roughness);
    free(map.total_field); free(meas); free(tl); free(to);
    printf("Done.\n");
    return 0;
}
