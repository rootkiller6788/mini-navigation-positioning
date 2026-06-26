/**
 * geomag_core.h -- Geomagnetic Navigation: Core Definitions
 *
 * L1: All fundamental definitions for geomagnetic field representation
 *   - Magnetic field vector components
 *   - Geodetic coordinates (WGS84)
 *   - Magnetic elements (declination, inclination, intensity)
 *   - Magnetometer measurement types
 *   - Navigation state vector
 *
 * Reference: Tsui, "Fundamentals of Global Positioning System Receivers" (2005)
 *            Merrill, "Magnetic Navigation" (2010)
 *            NOAA, "The World Magnetic Model" (2020)
 */

#ifndef GEOMAG_CORE_H
#define GEOMAG_CORE_H

#include <stddef.h>
#include <stdint.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================================
 * L1 Definition: Earth reference ellipsoid -- WGS84
 *
 * The WGS84 ellipsoid defines the shape of the Earth used for all geodetic
 * coordinate transforms. Constants from NIMA TR8350.2.
 * ============================================================================ */
#define WGS84_A     6378137.0       /* semi-major axis [m]                  */
#define WGS84_F     (1.0 / 298.257223563) /* flattening                     */
#define WGS84_B     (WGS84_A * (1.0 - WGS84_F))  /* semi-minor axis [m]    */
#define WGS84_E2    (1.0 - (WGS84_B * WGS84_B) / (WGS84_A * WGS84_A))
#define WGS84_E     sqrt(WGS84_E2)
#define WGS84_EP2   (WGS84_E2 / (1.0 - WGS84_E2))

#define DEG2RAD     (M_PI / 180.0)
#define RAD2DEG     (180.0 / M_PI)
#define ARCSEC2RAD  (M_PI / 648000.0)

/* L1: IGRF-13 reference radius for spherical harmonic expansion */
#define GEOMAG_REF_RADIUS  6371200.0    /* IGRF reference radius [m]       */
#define GEOMAG_EARTH_RADIUS 6371000.0   /* Mean Earth radius for mag [m]   */
#define EARTH_CORE_RADIUS   3485000.0   /* Approx outer core radius [m]    */

/* ============================================================================
 * L1 Definition: Geodetic coordinate (WGS84 ellipsoidal)
 *   lat: latitude [deg], -90 to +90
 *   lon: longitude [deg], -180 to +180
 *   alt: height above ellipsoid [m]
 * ============================================================================ */
typedef struct {
    double lat;
    double lon;
    double alt;
} GeodeticCoord;

/* L1: ECEF Cartesian coordinate */
typedef struct {
    double x, y, z;
} ECEFCoord;

/* L1: NED local tangent plane coordinate */
typedef struct {
    double n, e, d;
} NEDCoord;

/* ============================================================================
 * L1 Definition: Magnetic field vector [nT]
 *
 * Represents geomagnetic field. Units: nanoTesla (1 nT = 1e-9 T).
 * Surface range: ~22,000-68,000 nT total intensity.
 * ============================================================================ */
typedef struct {
    double bx;
    double by;
    double bz;
} MagVector;

/* ============================================================================
 * L1 Definition: Magnetic elements -- seven standard descriptors
 *
 * D = declination [deg], I = inclination [deg], F = total intensity [nT]
 * H = horizontal intensity [nT], X/Y/Z = NED components [nT]
 * ============================================================================ */
typedef struct {
    double declination;
    double inclination;
    double total_intensity;
    double horizontal;
    double north_component;
    double east_component;
    double vertical;
} MagneticElements;

/* L1: Spherical harmonic (Gauss) coefficient */
typedef struct {
    int     n;
    int     m;
    double  g_nm;       /* [nT] cosine coefficient                         */
    double  h_nm;       /* [nT] sine coefficient                           */
    double  dg_nm;      /* [nT/yr] secular variation of g                  */
    double  dh_nm;      /* [nT/yr] secular variation of h                  */
} GaussCoeff;

/* L1: IGRF model descriptor */
typedef struct {
    int     igrf_version;
    double  epoch;
    int     nmax;
    GaussCoeff *coeffs;
    int     ncoeffs;
    char    model_name[32];
} IGRFModel;

/* L1: Magnetometer types */
typedef enum {
    MAG_SENSOR_SCALAR,
    MAG_SENSOR_TRIAXIAL,
    MAG_SENSOR_GRADIOMETER
} MagnetometerType;

typedef struct {
    MagnetometerType type;
    double  measurement[3];
    double  timestamp;
    double  noise_std;
    double  bias[3];
    double  scale_factor[3];
    double  misalignment[9];
} MagMeasurement;

/* L1: Magnetic anomaly */
typedef struct {
    MagVector   total_field;
    MagVector   main_field;
    MagVector   anomaly;
    double      anomaly_magnitude;
    GeodeticCoord location;
} MagneticAnomaly;

/* L1: Navigation solution (9-DOF state) */
typedef struct {
    GeodeticCoord position;
    NEDCoord      velocity;
    double        roll;
    double        pitch;
    double        yaw;
    double        timestamp;
} NavSolution;

/* L1: Magnetic map grid for MAGCOM navigation */
typedef struct {
    GeodeticCoord origin;
    double        grid_spacing_lat;
    double        grid_spacing_lon;
    int           nlat;
    int           nlon;
    double       *total_field;
    double       *inclination;
    double       *declination;
    double       *anomaly;
} MagneticMap;

/* L1: Geomagnetic activity indices (Kp, Dst, AE) */
typedef struct {
    double kp;
    double dst;
    double ae;
    double timestamp;
} GeomagneticActivity;

/* L1: Legendre computation state */
typedef struct {
    int      nmax;
    double  *values;
    double  *derivatives;
    double   sin_theta;
    double   cos_theta;
} LegendreState;

/* ---- Core API ---- */

void geodetic_to_ecef(const GeodeticCoord *geo, ECEFCoord *ecef);
void ecef_to_geodetic(const ECEFCoord *ecef, GeodeticCoord *geo);
void ecef_to_ned_rotation(double lat, double lon, double R[9]);
void mag_ecef_to_ned(const MagVector *B_ecef, double lat, double lon, MagVector *B_ned);
void compute_magnetic_elements(const MagVector *B, MagneticElements *elem);
double mag_magnitude(const MagVector *B);
double mag_horizontal(const MagVector *B);
void magnetometer_calibrate(const double raw[3], const double bias[3],
                            const double scale_inv[9], double cal[3]);

#endif /* GEOMAG_CORE_H */
