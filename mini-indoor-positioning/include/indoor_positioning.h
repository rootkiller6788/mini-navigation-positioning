/**
 * @file indoor_positioning.h
 * @brief Core definitions for indoor positioning and navigation systems
 *
 * Knowledge Coverage:
 *   L1 - Definitions: coordinate systems, position types, measurement types,
 *        positioning modes, accuracy/precision metrics, reference frames
 *   L2 - Core Concepts: dead reckoning, trilateration, fingerprinting,
 *        sensor fusion, map matching, proximity detection
 *   L3 - Mathematical Structures: Euclidean/NED coordinates, quaternions,
 *        covariance matrices, state vectors
 *   L4 - Fundamental Laws: Free-space path loss, Bayesian estimation
 *
 * Reference: Tsui, "Fundamentals of Global Positioning System Receivers" (2005)
 *            Groves, "Principles of GNSS, Inertial, and Multisensor Navigation" (2013)
 *
 * Course Alignment:
 *   - MIT 6.003 (Signal Processing)
 *   - Stanford EE359 (Wireless)
 *   - ETH 227-0436 (Communications)
 *   - 清华 导航与定位 (Navigation & Positioning)
 */

#ifndef INDOOR_POSITIONING_H
#define INDOOR_POSITIONING_H

#define _USE_MATH_DEFINES
#include <stdint.h>
#include <stddef.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Speed of light constant shared across modules */
#ifndef SPEED_OF_LIGHT_MPS
#define SPEED_OF_LIGHT_MPS 299792458.0
#endif

/* ============================================================================
 * L1 - Definitions: Coordinate Systems
 * ============================================================================ */

/** Number of spatial dimensions (2D or 3D indoor positioning) */
#define IP_DIMENSIONS_2D  2
#define IP_DIMENSIONS_3D  3
#define IP_MAX_ANCHORS   32
#define IP_MAX_FINGERPRINTS 4096
#define IP_STATE_DIM      9   /**< Position(3)+Velocity(3)+Attitude(3) */
#define IP_MEASUREMENT_DIM 6  /**< Typical IMU: accel(3)+gyro(3) */

/**
 * @brief 2D position in Cartesian coordinates (east, north)
 *
 * L1 Definition: Indoor position is typically expressed in a local
 * tangent plane coordinate system (ENU: East-North-Up).
 */
typedef struct {
    double x;  /**< East coordinate in meters */
    double y;  /**< North coordinate in meters */
} position2d_t;

/**
 * @brief 3D position in ENU (East-North-Up) coordinates
 *
 * L1 Definition: 3D indoor positioning extends the 2D plane with
 * vertical (floor/altitude) information.
 * Units: meters.
 * Reference frame: local tangent plane, origin at building reference point.
 */
typedef struct {
    double x;  /**< East coordinate in meters */
    double y;  /**< North coordinate in meters */
    double z;  /**< Up coordinate in meters (height above reference floor) */
} position3d_t;

/**
 * @brief 2D velocity vector
 *
 * L1 Definition: Velocity is the time derivative of position.
 * Units: meters per second.
 */
typedef struct {
    double vx; /**< East velocity in m/s */
    double vy; /**< North velocity in m/s */
} velocity2d_t;

/**
 * @brief 3D velocity vector
 */
typedef struct {
    double vx; /**< East velocity in m/s */
    double vy; /**< North velocity in m/s */
    double vz; /**< Up velocity in m/s */
} velocity3d_t;

/**
 * @brief Quaternion for 3D attitude representation
 *
 * L1 Definition: Attitude (orientation) is represented as a unit quaternion
 * q = w + xi + yj + zk, with ||q|| = 1.
 * Quaternions avoid gimbal lock and are computationally efficient.
 */
typedef struct {
    double w;  /**< Scalar (real) part */
    double x;  /**< i component */
    double y;  /**< j component */
    double z;  /**< k component */
} quaternion_t;

/**
 * @brief Full navigation state vector (position, velocity, attitude)
 *
 * L2 Concept: The navigation state encodes all information needed
 * to describe the motion of a mobile agent in 3D space.
 *
 * State vector: [x, y, z, vx, vy, vz, roll, pitch, yaw]
 * Attitude stored as Euler angles for simplicity (internally may use quaternion).
 */
typedef struct {
    position3d_t pos;    /**< Position in ENU frame */
    velocity3d_t vel;    /**< Velocity in ENU frame */
    double roll;         /**< Roll angle in radians */
    double pitch;        /**< Pitch angle in radians */
    double yaw;          /**< Yaw (heading) angle in radians */
} navigation_state_t;

/**
 * @brief Error ellipse parameters for 2D positioning uncertainty
 *
 * L1 Definition: Positioning accuracy is characterized by the
 * error ellipse: semi-major axis (a), semi-minor axis (b), orientation (theta).
 * CEP = 0.59*(a+b) for circular error probable.
 */
typedef struct {
    double semi_major;   /**< Semi-major axis length in meters */
    double semi_minor;   /**< Semi-minor axis length in meters */
    double orientation;  /**< Orientation of major axis in radians (from east) */
    double confidence;   /**< Confidence level (e.g., 0.95 for 95%) */
} error_ellipse_t;

/* ============================================================================
 * L1 - Definitions: Measurement Types
 * ============================================================================ */

/** Measurement source type */
typedef enum {
    MEAS_SOURCE_RSSI      = 0,  /**< Received Signal Strength Indicator */
    MEAS_SOURCE_TOF       = 1,  /**< Time of Flight */
    MEAS_SOURCE_TDOA      = 2,  /**< Time Difference of Arrival */
    MEAS_SOURCE_AOA       = 3,  /**< Angle of Arrival */
    MEAS_SOURCE_IMU_ACCEL = 4,  /**< IMU accelerometer */
    MEAS_SOURCE_IMU_GYRO  = 5,  /**< IMU gyroscope */
    MEAS_SOURCE_MAG       = 6,  /**< Magnetometer */
    MEAS_SOURCE_BARO      = 7,  /**< Barometer (altitude) */
    MEAS_SOURCE_UWB       = 8,  /**< Ultra-Wideband ranging */
    MEAS_SOURCE_ODOM      = 9,  /**< Odometer / step counter */
    MEAS_SOURCE_CAMERA    = 10, /**< Visual odometry */
    MEAS_SOURCE_BLE       = 11, /**< Bluetooth Low Energy beacon */
    MEAS_SOURCE_WIFI_FP   = 12, /**< WiFi fingerprint */
    MEAS_SOURCE_MAG_FP    = 13, /**< Magnetic field fingerprint */
    MEAS_SOURCE_LIDAR     = 14, /**< LiDAR scan matching */
    MEAS_SOURCE_COUNT     = 15
} measurement_source_t;

/**
 * @brief Generic measurement structure
 *
 * L1 Definition: A measurement is an observation z = h(x) + v,
 * where h is the measurement function, x is the state, v is noise.
 */
typedef struct {
    measurement_source_t source;  /**< Source type */
    uint64_t timestamp_us;        /**< Measurement timestamp in microseconds */
    double value[4];              /**< Measurement values (up to 4 components) */
    double noise_std;             /**< Measurement noise standard deviation */
    uint8_t anchor_id;            /**< Anchor/beacon ID if applicable */
    uint8_t n_components;         /**< Number of valid value components */
} measurement_t;

/* ============================================================================
 * L2 - Core Concepts: Positioning Modes
 * ============================================================================ */

/** Positioning estimation mode */
typedef enum {
    POS_MODE_NONE           = 0,  /**< No position estimate */
    POS_MODE_PROXIMITY      = 1,  /**< Nearest beacon */
    POS_MODE_TRILATERATION  = 2,  /**< Distance-based trilateration */
    POS_MODE_MULTILATERATION = 3, /**< TDOA-based multilateration */
    POS_MODE_FINGERPRINT    = 4,  /**< RSSI/feature fingerprint matching */
    POS_MODE_DEAD_RECKONING = 5,  /**< Inertial dead reckoning */
    POS_MODE_KALMAN_FUSION  = 6,  /**< Kalman filter sensor fusion */
    POS_MODE_PARTICLE_FILTER = 7, /**< Particle filter fusion */
    POS_MODE_HYBRID         = 8,  /**< Hybrid multi-modal positioning */
} positioning_mode_t;

/**
 * @brief Map constraints for indoor positioning
 *
 * L2 Concept: Indoor positioning benefits from map information
 * (wall constraints, floor plans, forbidden regions).
 * Known as "map-aided" or "map-matched" positioning.
 */
typedef struct {
    double x_min;          /**< Building east boundary in meters */
    double x_max;          /**< Building west boundary in meters */
    double y_min;          /**< Building south boundary in meters */
    double y_max;          /**< Building north boundary in meters */
    double floor_height;   /**< Height per floor in meters */
    int    n_floors;       /**< Number of floors */
    double wall_points[IP_MAX_ANCHORS][4]; /**< Wall segments: (x0,y0,x1,y1) */
    int    n_walls;        /**< Number of wall segments */
} map_constraint_t;

/* ============================================================================
 * L1 - Accuracy and Precision Metrics
 * ============================================================================ */

/** Positioning accuracy metrics */
typedef struct {
    double mean_error_2d;        /**< Mean 2D error in meters */
    double mean_error_3d;        /**< Mean 3D error in meters */
    double rmse_2d;             /**< Root-mean-square error in 2D */
    double rmse_3d;             /**< Root-mean-square error in 3D */
    double cep50;               /**< Circular Error Probable 50% */
    double cep95;               /**< Circular Error Probable 95% */
    double cep99;               /**< Circular Error Probable 99% */
    double drms;                /**< Distance RMS error */
    double max_error;           /**< Maximum observed error */
    int    n_samples;           /**< Number of samples collected */
} positioning_accuracy_t;

/* ============================================================================
 * L3 - Mathematical Structures: Matrices and Vectors
 * ============================================================================ */

/**
 * @brief 2x2 matrix for 2D covariance
 */
typedef struct {
    double m[2][2];
} matrix2_t;

/**
 * @brief 3x3 matrix for 3D operations
 */
typedef struct {
    double m[3][3];
} matrix3_t;

/**
 * @brief 9x9 covariance matrix for navigation state
 *
 * L3 Structure: The covariance matrix encodes the uncertainty of
 * the state estimate. Diagonal elements are variances, off-diagonal
 * elements are covariances between state components.
 */
typedef struct {
    double m[IP_STATE_DIM][IP_STATE_DIM];
} covariance_matrix_t;

/**
 * @brief Dynamic measurement matrix (rows x columns, bounded)
 */
typedef struct {
    double data[IP_STATE_DIM][IP_STATE_DIM];
    int rows;
    int cols;
} matrix_t;

/* ============================================================================
 * L4 - Fundamental Laws: Free-Space Path Loss Model
 * ============================================================================ */

/**
 * @brief RSSI path loss model parameters
 *
 * L4 Law: Friis free-space equation adapted for indoor:
 *   RSSI(d) = RSSI_0 - 10 * n * log10(d / d_0) + X_sigma
 * where:
 *   RSSI_0 = reference RSSI at distance d_0 (typically 1m)
 *   n = path loss exponent (2.0 free space, 1.6-1.8 corridor, 4-6 office)
 *   X_sigma ~ N(0, sigma^2) shadow fading term
 */
typedef struct {
    double rssi_at_1m;    /**< RSSI at 1 meter reference distance (dBm) */
    double path_loss_exp; /**< Path loss exponent n */
    double shadow_std;    /**< Shadow fading standard deviation sigma (dB) */
    double frequency_mhz; /**< Carrier frequency in MHz */
} path_loss_model_t;

/**
 * @brief Convert RSSI to distance estimate using path loss model
 *
 * @param rssi Measured RSSI in dBm
 * @param model Path loss model parameters
 * @return Estimated distance in meters, or -1.0 on error
 *
 * Formula: d = d_0 * 10^((RSSI_0 - rssi) / (10 * n))
 *
 * L4: This is the Friis transmission equation adapted for
 * non-free-space environments through the path loss exponent.
 * Reference: Molisch, "Wireless Communications" (2011), Ch.4.
 *
 * Complexity: O(1), Time: 1 pow10 + 1 division
 */
double rssi_to_distance(double rssi, const path_loss_model_t *model);

/**
 * @brief Compute expected RSSI at given distance
 *
 * @param distance Distance in meters
 * @param model Path loss model parameters
 * @return Expected RSSI in dBm
 *
 * Complexity: O(1)
 */
double distance_to_rssi(double distance, const path_loss_model_t *model);

/* ============================================================================
 * L5 - Algorithms: Position Estimation
 * ============================================================================ */

/**
 * @brief Solve 2D trilateration problem using linear least squares
 *
 * Given distances to N anchors at known positions, estimate the
 * receiver position that minimizes squared distance errors.
 *
 * Algorithm: Linearize the N quadratic distance equations by
 * subtracting one anchor equation from all others, then solve
 * the linear system A*x = b via pseudo-inverse.
 *
 * Requires at least 3 anchors in 2D (3 equations -> 2 unknowns).
 *
 * @param anchor_positions Array of known anchor positions (at least 3)
 * @param distances Array of measured distances to each anchor
 * @param n_anchors Number of anchors (>= 3 for 2D)
 * @param[out] result Estimated position
 * @return 0 on success, -1 on numerical error or insufficient anchors
 *
 * L5: Linearized least-squares trilateration.
 * Reference: Langendoen & Reijers, "Distributed localization in wireless
 * sensor networks: a quantitative comparison" (2003).
 *
 * Complexity: O(N) for N anchors, Time: O(N)
 */
int trilateration_2d(const position2d_t *anchor_positions,
                     const double *distances,
                     int n_anchors,
                     position2d_t *result);

/**
 * @brief Solve 3D trilateration using non-linear least squares (Gauss-Newton)
 *
 * @param anchor_positions Array of known anchor positions (at least 4)
 * @param distances Array of measured distances to each anchor
 * @param n_anchors Number of anchors (>= 4 for 3D)
 * @param initial_guess Starting point for iterative solver
 * @param[out] result Estimated position
 * @param max_iterations Maximum Gauss-Newton iterations
 * @param tolerance Convergence tolerance
 * @return 0 on success, -1 on failure
 *
 * Complexity: O(K*N) where K iterations, N anchors
 */
int trilateration_3d(const position3d_t *anchor_positions,
                     const double *distances,
                     int n_anchors,
                     const position3d_t *initial_guess,
                     position3d_t *result,
                     int max_iterations,
                     double tolerance);

/**
 * @brief Solve TDOA multilateration using Chan's algorithm
 *
 * Given time-difference-of-arrival measurements with respect to a
 * reference anchor, estimate the receiver position.
 *
 * Algorithm: Chan's two-step weighted least squares (Chan & Ho, 1994).
 * Solves hyperbolic positioning equations.
 *
 * @param anchor_positions Anchor positions (reference anchor at index 0)
 * @param tdoa_measurements TDOA values in seconds relative to anchor[0]
 * @param n_anchors Total number of anchors (>= 4 for 3D TDOA)
 * @param speed_of_light Signal propagation speed in m/s
 * @param[out] result Estimated position
 * @return 0 on success, -1 on failure
 *
 * L5: Chan's algorithm for hyperbolic positioning.
 * Reference: Chan & Ho, "A simple and efficient estimator for hyperbolic
 * location," IEEE Trans. Signal Processing, 1994.
 *
 * Complexity: O(N^2) for N anchors
 */
int tdoa_multilateration(const position3d_t *anchor_positions,
                         const double *tdoa_measurements,
                         int n_anchors,
                         double speed_of_light,
                         position3d_t *result);

/**
 * @brief Weighted centroid localization
 *
 * Simple proximity-based positioning: weighted average of anchor
 * positions, with weights inversely proportional to distance.
 *
 * @param anchor_positions Known anchor positions
 * @param distances Measured distances to each anchor
 * @param n_anchors Number of anchors
 * @param[out] result Weighted centroid position
 * @return 0 on success
 *
 * L5: Weighted centroid localization.
 * Reference: Bulusu et al., "GPS-less low-cost outdoor localization
 * for very small devices," IEEE Personal Comm., 2000.
 *
 * Complexity: O(N)
 */
int weighted_centroid_2d(const position2d_t *anchor_positions,
                         const double *distances,
                         int n_anchors,
                         position2d_t *result);

/* ============================================================================
 * L1 - Coordinate Transformations
 * ============================================================================ */

/**
 * @brief Convert latitude/longitude/altitude to ENU coordinates
 *
 * Uses a reference point as the origin of the local tangent plane.
 *
 * @param lat Observer latitude in radians
 * @param lon Observer longitude in radians
 * @param alt Observer altitude in meters
 * @param ref_lat Reference point latitude in radians
 * @param ref_lon Reference point longitude in radians
 * @param ref_alt Reference point altitude in meters
 * @param[out] enu Result in East-North-Up meters
 *
 * L3: Coordinate transformation between geodetic and local
 * tangent plane coordinates. Uses WGS-84 ellipsoid parameters.
 */
void geodetic_to_enu(double lat, double lon, double alt,
                     double ref_lat, double ref_lon, double ref_alt,
                     position3d_t *enu);

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

/**
 * @brief Compute Euclidean distance between two 2D points
 *
 * @return Distance in meters
 * Complexity: O(1)
 */
static inline double distance_2d(position2d_t a, position2d_t b) {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    return sqrt(dx * dx + dy * dy);
}

/**
 * @brief Compute Euclidean distance between two 3D points
 *
 * @return Distance in meters
 * Complexity: O(1)
 */
static inline double distance_3d(position3d_t a, position3d_t b) {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    double dz = a.z - b.z;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

/**
 * @brief Normalize a quaternion to unit length
 *
 * @param q Quaternion to normalize (modified in place)
 */
static inline void quaternion_normalize(quaternion_t *q) {
    double norm = sqrt(q->w * q->w + q->x * q->x + q->y * q->y + q->z * q->z);
    if (norm > 1e-12) {
        q->w /= norm; q->x /= norm; q->y /= norm; q->z /= norm;
    }
}

/**
 * @brief Convert quaternion to Euler angles (roll, pitch, yaw)
 *
 * Assumes ZYX convention (yaw-pitch-roll).
 *
 * @param q Input quaternion
 * @param[out] roll Roll angle in radians
 * @param[out] pitch Pitch angle in radians
 * @param[out] yaw Yaw angle in radians
 */
void quaternion_to_euler(const quaternion_t *q,
                         double *roll, double *pitch, double *yaw);

/**
 * @brief Convert Euler angles to quaternion
 *
 * @param roll Roll angle in radians
 * @param pitch Pitch angle in radians
 * @param yaw Yaw angle in radians
 * @param[out] q Resulting quaternion
 */
void euler_to_quaternion(double roll, double pitch, double yaw, quaternion_t *q);

/**
 * @brief Rotate a 3D vector by a quaternion
 *
 * @param v Input vector (x, y, z components in array)
 * @param q Rotation quaternion
 * @param[out] result Rotated vector
 */
void quaternion_rotate_vector(const double v[3], const quaternion_t *q,
                              double result[3]);

#endif /* INDOOR_POSITIONING_H */
