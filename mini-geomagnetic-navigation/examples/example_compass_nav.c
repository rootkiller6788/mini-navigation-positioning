/**
 * example_compass_nav.c -- Magnetic compass dead-reckoning demo
 */
#include "geomag_core.h"
#include "geomag_model.h"
#include "geomag_navigation.h"
#include "geomag_math.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== Compass Navigation Demo ===\n");
    IGRFModel model;
    igrf_init_model(13, &model);
    NavSolution state;
    state.position.lat = 37.7749; state.position.lon = -122.4194;
    state.position.alt = 0.0; state.timestamp = 0.0;
    GeodeticCoord loc = {state.position.lat, state.position.lon, 0.0};
    MagVector B; MagneticElements elem;
    igrf_compute_field(&model, &loc, &B);
    compute_magnetic_elements(&B, &elem);
    printf("Start: (%.4f, %.4f), D=%.2f, I=%.2f, F=%.1f nT\n",
           state.position.lat, state.position.lon,
           elem.declination, elem.inclination, elem.total_intensity);
    double speeds[] = {5.0, 5.0, 5.0, 5.0};
    double hdgs[] = {90.0, 0.0, 270.0, 180.0};
    double durs[] = {200.0, 200.0, 200.0, 200.0};
    for (int i = 0; i < 4; i++) {
        mag_compass_dr_update(&state, hdgs[i], speeds[i], durs[i]);
        printf("Leg %d: hdg=%.0f -> (%.4f, %.4f)\n",
               i+1, hdgs[i], state.position.lat, state.position.lon);
    }
    double dlat = (state.position.lat - 37.7749) * 111000.0;
    double dlon = (state.position.lon - (-122.4194)) * 111000.0 * cos(37.7749*M_PI/180.0);
    printf("Closure error: %.1f m\n", sqrt(dlat*dlat + dlon*dlon));
    MagVector Bb = {20000.0, 2000.0, 35000.0};
    double lev = mag_heading_from_triaxial(&Bb, elem.declination);
    double tilt = tilt_compensated_heading(&Bb, 0.0, 10.0*M_PI/180.0, elem.declination);
    printf("Heading: level=%.2f, tilt=%.2f diff=%.2f deg\n", lev, tilt, lev-tilt);
    igrf_free_model(&model);
    printf("Done.\n");
    return 0;
}
