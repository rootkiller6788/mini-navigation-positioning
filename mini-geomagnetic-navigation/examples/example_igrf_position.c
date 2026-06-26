/**
 * example_igrf_position.c -- IGRF single-point positioning demo
 */
#include "geomag_model.h"
#include "geomag_navigation.h"
#include "geomag_kalman.h"
#include "geomag_math.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== IGRF Positioning Demo ===\n");
    IGRFModel model;
    if (igrf_init_model(13, &model) != 0) {
        printf("Failed to load IGRF-13.\n");
        return 1;
    }
    printf("IGRF-13 loaded: nmax=%d, %d coeffs\n", model.nmax, model.ncoeffs);
    GeodeticCoord truth = {40.015, -105.270, 1600.0};
    MagVector B; MagneticElements elem;
    igrf_compute_field(&model, &truth, &B);
    compute_magnetic_elements(&B, &elem);
    printf("Truth: (%.3f, %.3f)\n", truth.lat, truth.lon);
    printf("F=%.1f nT, D=%.2f deg, I=%.2f deg\n",
           elem.total_intensity, elem.declination, elem.inclination);
    double elat, elon, residual;
    igrf_single_point_position(&model, elem.total_intensity, elem.inclination,
                                1.0, 1600.0, 40.0, -100.0, 10.0,
                                &elat, &elon, &residual);
    printf("Estimated: (%.3f, %.3f), residual=%.2f\n", elat, elon, residual);
    printf("Error: %.1f km\n",
           great_circle_distance(&truth, &(GeodeticCoord){elat, elon, 0.0})/1000.0);
    double m = compute_dipole_moment(&model);
    printf("Dipole moment: %.2e A*m^2\n", m);
    double Fn[3];
    predict_magnetic_measurement(&model, truth.lat, truth.lon, truth.alt, Fn);
    printf("Predicted: F=%.1f D=%.2f I=%.2f\n", Fn[0], Fn[1], Fn[2]);
    igrf_free_model(&model);
    printf("Done.\n");
    return 0;
}
