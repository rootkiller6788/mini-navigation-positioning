/**
 * @file inertial_navigation.h
 * @brief Inertial navigation and dead reckoning for indoor positioning
 *
 * Knowledge Coverage:
 *   L1 - Definitions: IMU, accelerometer, gyroscope, magnetometer,
 *        specific force, angular rate, bias, scale factor, misalignment
 *   L2 - Core Concepts: Strapdown inertial navigation, dead reckoning,
 *        step detection, stride length estimation, heading estimation
 *   L3 - Mathematical Structures: Quaternion kinematics, rotation matrices,
 *        Euler angles, DCM (Direction Cosine Matrix)
 *   L4 - Fundamental Laws: Newton's laws in non-inertial frames,
 *        Coriolis effect, Schuler tuning, INS error dynamics
 *   L5 - Algorithms: Strapdown integration, ZUPT (Zero Velocity Update),
 *        step-and-heading dead reckoning, Madgwick/Mahony AHRS filters
 *   L6 - Canonical Problems: Pedestrian dead reckoning (PDR) in buildings,
 *        foot-mounted IMU with ZUPT for indoor positioning
 *
 * Reference: Groves, "Principles of GNSS, Inertial, and Multisensor
 *            Integrated Navigation Systems" (2013)
 *            Titterton & Weston, "Strapdown Inertial Navigation Technology" (2004)
 *            Madgwick, "An efficient orientation filter for IMU and MARG arrays" (2010)
 *
 * Course Alignment:
 *   - MIT 6.003 (Signal Processing) — integration/drift
 *   - Stanford EE267 (Virtual Reality) — IMU orientation
 *   - Berkeley EE16A/B (Circuits) — sensor interfacing
 *   - Michigan EECS 411 (Microwave) — sensor technology
 *   - 清华 导航与定位
 */

#ifndef INERTIAL_NAVIGATION_H
#define INERTIAL_NAVIGATION_H

#include "indoor_positioning.h"

/* ============================================================================
 * L1 - Definitions: IMU Data Structures
 * ============================================================================ */

/** IMU sample rate constants */
#define IMU_MAX_SAMPLES       10000
#define IMU_GYRO_SCALE_FACTOR 0.017453292519943295 /**< deg/s to rad/s */
#define IMU_GRAVITY_MSS       9.80665              /**< Standard gravity m/s^2 */

/**
 * @brief IMU sensor specifications
 *
 * L1 Definition: IMU (Inertial Measurement Unit) consists of:
 * - 3-axis accelerometer (measures specific force in m/s^2)
 * - 3-axis gyroscope (measures angular rate in rad/s)
 * - Optionally 3-axis magnetometer (measures magnetic field in uT)
 */
typedef struct {
    double gyro_noise_density;     /**< Gyro angular random walk (rad/s/sqrt(Hz)) */
    double gyro_bias_stability;    /**< Gyro bias instability (rad/s) */
    double accel_noise_density;    /**< Accel velocity random walk (m/s^2/sqrt(Hz)) */
    double accel_bias_stability;   /**< Accel bias instability (m/s^2) */
    double sample_rate_hz;         /**< Nominal sample rate in Hz */
    double gyro_full_scale_radps;  /**< Gyro full-scale range */
    double accel_full_scale_mss;   /**< Accelerometer full-scale range */
    int    has_magnetometer;       /**< 1 if magnetometer present */
    double mag_noise_density;      /**< Magnetometer noise (uT/sqrt(Hz)) */
} imu_specs_t;

/**
 * @brief Raw IMU measurement sample
 *
 * L1 Definition: A single IMU measurement comprises triaxial
 * accelerometer and gyroscope readings at a specific timestamp.
 */
typedef struct {
    double accel_x;       /**< X-axis specific force in m/s^2 (body frame) */
    double accel_y;       /**< Y-axis specific force in m/s^2 (body frame) */
    double accel_z;       /**< Z-axis specific force in m/s^2 (body frame) */
    double gyro_x;        /**< X-axis angular rate in rad/s (body frame) */
    double gyro_y;        /**< Y-axis angular rate in rad/s (body frame) */
    double gyro_z;        /**< Z-axis angular rate in rad/s (body frame) */
    double mag_x;         /**< X-axis magnetic field in uT (body frame, optional) */
    double mag_y;         /**< Y-axis magnetic field in uT (body frame, optional) */
    double mag_z;         /**< Z-axis magnetic field in uT (body frame, optional) */
    uint64_t timestamp_us; /**< Sample timestamp in microseconds */
    uint8_t mag_valid;     /**< 1 if magnetometer data is valid */
} imu_sample_t;

/**
 * @brief IMU calibration parameters
 *
 * L1 Definition: IMU calibration corrects for:
 * - Bias (offset): constant error in the measurement
 * - Scale factor: deviation from nominal sensitivity
 * - Misalignment: non-orthogonality between sensor axes
 * The 6-position or multi-position calibration estimates these.
 */
typedef struct {
    double accel_bias[3];         /**< Accelerometer bias [x,y,z] in m/s^2 */
    double accel_scale[3];        /**< Accelerometer scale factor [x,y,z] */
    double accel_misalign[3][3];  /**< Accelerometer misalignment matrix */
    double gyro_bias[3];          /**< Gyroscope bias [x,y,z] in rad/s */
    double gyro_scale[3];         /**< Gyroscope scale factor [x,y,z] */
    double gyro_misalign[3][3];   /**< Gyroscope misalignment matrix */
    double mag_bias[3];           /**< Magnetometer hard-iron bias in uT */
    double mag_scale[3];          /**< Magnetometer soft-iron scale */
    double mag_misalign[3][3];    /**< Magnetometer misalignment matrix */
} imu_calibration_t;

/**
 * @brief Apply calibration to an IMU sample
 *
 * @param raw Raw IMU sample
 * @param cal Calibration parameters
 * @param[out] corrected Calibrated IMU sample
 *
 * Corrected = scale * (misalign^-1 * (raw - bias))
 */
void imu_apply_calibration(const imu_sample_t *raw,
                           const imu_calibration_t *cal,
                           imu_sample_t *corrected);

/* ============================================================================
 * L2 / L5 - Strapdown Inertial Navigation System (INS)
 * ============================================================================ */

/**
 * @brief Strapdown INS state
 *
 * L2 Concept: Strapdown INS computes position, velocity, and attitude
 * by integrating IMU measurements in the navigation frame.
 * The mechanization equations propagate the state forward in time.
 */
typedef struct {
    navigation_state_t nav;        /**< Current navigation state */
    quaternion_t       quat;       /**< Attitude quaternion (body-to-nav) */
    imu_sample_t       last_imu;   /**< Previous IMU sample for integration */
    double             dt;         /**< Integration time step in seconds */
    uint64_t           last_time;  /**< Timestamp of last update */
    double             pos_cov[3]; /**< Position uncertainty estimates */
} ins_state_t;

/**
 * @brief Initialize a strapdown INS with known initial position and heading
 *
 * @param ins INS state to initialize
 * @param initial_pos Known initial position
 * @param initial_heading Initial heading in radians (from north clockwise)
 * @param sample_rate_hz IMU sample rate
 */
void ins_init(ins_state_t *ins, const position3d_t *initial_pos,
              double initial_heading, double sample_rate_hz);

/**
 * @brief Strapdown INS mechanization update step
 *
 * Propagates the navigation state forward by one IMU sample period.
 * This is the core INS algorithm.
 *
 * Steps:
 *   1. Attitude update: integrate gyroscope to update quaternion
 *   2. Transform specific force to navigation frame
 *   3. Velocity update: subtract gravity, integrate acceleration
 *   4. Position update: integrate velocity
 *
 * @param ins INS state (updated in place)
 * @param imu Current IMU sample (calibrated)
 *
 * L5: Strapdown INS mechanization equations.
 * Reference: Groves (2013), Ch.5 - "Inertial Navigation".
 *
 * Complexity: O(1) per sample
 */
void ins_mechanize(ins_state_t *ins, const imu_sample_t *imu);

/**
 * @brief Apply Zero-Velocity Update (ZUPT) to correct INS drift
 *
 * When the IMU is known to be stationary (e.g., foot flat on ground),
 * measured velocity should be zero. ZUPT uses this pseudo-measurement
 * to estimate and correct IMU biases, bounding the cubic-in-time
 * position drift of unaided INS.
 *
 * @param ins INS state to correct
 * @param zupt_detected 1 if zero-velocity condition is detected
 *
 * L5: ZUPT for foot-mounted pedestrian INS.
 * Reference: Foxlin, "Pedestrian tracking with shoe-mounted inertial
 * sensors," IEEE CG&A, 2005.
 *            Nilsson et al., "A note on the limitations of ZUPTs and
 * the implications on sensor error modeling," IEEE IPIN, 2012.
 */
void ins_apply_zupt(ins_state_t *ins, int zupt_detected);

/**
 * @brief Detect zero-velocity condition from IMU data
 *
 * Uses stance-phase detection: accelerometer magnitude near gravity,
 * gyroscope magnitude near zero.
 *
 * @param imu IMU sample to evaluate
 * @param accel_thresh Acceleration magnitude deviation threshold (m/s^2)
 * @param gyro_thresh Angular rate magnitude threshold (rad/s)
 * @return 1 if zero-velocity detected, 0 otherwise
 *
 * L5: Stance phase detection for pedestrian INS.
 */
int detect_zero_velocity(const imu_sample_t *imu,
                         double accel_thresh, double gyro_thresh);

/* ============================================================================
 * L5 - AHRS: Attitude and Heading Reference System
 * ============================================================================ */

/**
 * @brief Madgwick AHRS filter state
 *
 * L5: Madgwick's gradient descent orientation filter fuses
 * gyroscope, accelerometer, and magnetometer data to estimate
 * the sensor's orientation as a quaternion.
 *
 * Reference: Madgwick, "An efficient orientation filter for inertial
 * and inertial/magnetic sensor arrays" (2010).
 */
typedef struct {
    quaternion_t q;          /**< Estimated orientation quaternion */
    double       beta;       /**< Algorithm gain (filter divergence rate) */
    double       sample_period; /**< Sample period in seconds */
} madgwick_ahrs_t;

/**
 * @brief Initialize Madgwick AHRS filter
 *
 * @param ahrs Filter state to initialize
 * @param beta Algorithm gain (typical: 0.033 for IMU, 0.041 for MARG)
 * @param sample_period Sample period in seconds
 */
void madgwick_init(madgwick_ahrs_t *ahrs, double beta, double sample_period);

/**
 * @brief Update Madgwick AHRS with IMU data
 *
 * @param ahrs Filter state (updated in place)
 * @param gx Gyroscope x in rad/s
 * @param gy Gyroscope y in rad/s
 * @param gz Gyroscope z in rad/s
 * @param ax Accelerometer x in m/s^2
 * @param ay Accelerometer y in m/s^2
 * @param az Accelerometer z in m/s^2
 *
 * L5: Madgwick IMU orientation filter.
 * Complexity: O(1)
 */
void madgwick_update_imu(madgwick_ahrs_t *ahrs,
                         double gx, double gy, double gz,
                         double ax, double ay, double az);

/**
 * @brief Update Madgwick AHRS with IMU + magnetometer data (MARG mode)
 *
 * @param ahrs Filter state (updated in place)
 * @param gx,gy,gz Gyroscope angular rates in rad/s
 * @param ax,ay,az Accelerometer specific force in m/s^2
 * @param mx,my,mz Magnetometer field in uT
 *
 * L5: Madgwick MARG (Magnetic, Angular Rate, Gravity) orientation filter.
 * Complexity: O(1)
 */
void madgwick_update_marg(madgwick_ahrs_t *ahrs,
                          double gx, double gy, double gz,
                          double ax, double ay, double az,
                          double mx, double my, double mz);

/**
 * @brief Mahony AHRS filter state
 *
 * L5: Mahony's complementary filter uses a PI controller to fuse
 * gyroscope integration with accelerometer/magnetometer corrections.
 *
 * Reference: Mahony et al., "Nonlinear complementary filters on the
 * special orthogonal group," IEEE Trans. Automatic Control, 2008.
 */
typedef struct {
    quaternion_t q;      /**< Estimated orientation quaternion */
    double kp;           /**< Proportional gain */
    double ki;           /**< Integral gain */
    double integral_fb[3]; /**< Integral error feedback */
    double sample_period;  /**< Sample period in seconds */
} mahony_ahrs_t;

/**
 * @brief Initialize Mahony AHRS filter
 */
void mahony_init(mahony_ahrs_t *ahrs, double kp, double ki, double sample_period);

/**
 * @brief Update Mahony AHRS with IMU data
 */
void mahony_update_imu(mahony_ahrs_t *ahrs,
                       double gx, double gy, double gz,
                       double ax, double ay, double az);

/* ============================================================================
 * L6 - Pedestrian Dead Reckoning (PDR)
 * ============================================================================ */

/**
 * @brief Step detector state for pedestrian dead reckoning
 *
 * L6: PDR estimates position by detecting steps and estimating
 * step length and heading. It is a low-cost alternative to
 * full strapdown INS for smartphone-based indoor positioning.
 *
 * Reference: Weinberg, "Using the ADXL202 in pedometer and personal
 * navigation applications," Analog Devices AN-602, 2002.
 */
typedef struct {
    double step_threshold;     /**< Acceleration magnitude threshold */
    double step_timeout_ms;    /**< Minimum time between steps in ms */
    double last_step_time;     /**< Timestamp of last detected step */
    double accel_mag_prev;     /**< Previous acceleration magnitude */
    double stride_length;      /**< Estimated stride length in meters */
    int    step_count;         /**< Total steps detected */
} step_detector_t;

/**
 * @brief Pedestrian dead reckoning state
 */
typedef struct {
    position2d_t      position;      /**< Current 2D position */
    double            heading;       /**< Current heading in radians */
    step_detector_t   step_det;      /**< Step detector state */
    double            total_distance; /**< Total distance traveled in meters */
} pdr_state_t;

/**
 * @brief Initialize a PDR state
 *
 * @param pdr PDR state to initialize
 * @param initial_pos Starting position
 * @param initial_heading Starting heading in radians
 * @param stride_length Estimated stride length in meters
 */
void pdr_init(pdr_state_t *pdr, position2d_t initial_pos,
              double initial_heading, double stride_length);

/**
 * @brief Process an accelerometer sample for step detection and PDR update
 *
 * @param pdr PDR state (updated in place)
 * @param accel_magnitude Acceleration magnitude in m/s^2
 * @param heading Current heading estimate in radians
 * @param timestamp_ms Current timestamp in milliseconds
 * @return 1 if a step was detected, 0 otherwise
 *
 * L6: Step detection via acceleration magnitude threshold crossing.
 * Uses peak detection on the bandpass-filtered acceleration magnitude.
 *
 * Complexity: O(1)
 */
int pdr_process_accel(pdr_state_t *pdr, double accel_magnitude,
                      double heading, double timestamp_ms);

/**
 * @brief Tune the stride length model based on acceleration statistics
 *
 * Weinberg model: stride_length = K * (a_max - a_min)^(1/4)
 * where a_max, a_min are peak/trough vertical accelerations.
 *
 * @param accel_max Maximum vertical acceleration in recent window
 * @param accel_min Minimum vertical acceleration in recent window
 * @param K Calibration constant (typical: 0.4-0.5)
 * @return Estimated stride length in meters
 *
 * Reference: Weinberg (2002)
 */
double pdr_stride_length_weinberg(double accel_max, double accel_min, double K);

/**
 * @brief Alternative stride length model (Kim model)
 *
 * stride_length = K * (sum|a| / N)^(1/3)
 * where sum|a| is average absolute acceleration.
 *
 * @param avg_abs_accel Average absolute acceleration in m/s^2
 * @param K Calibration constant
 * @return Estimated stride length in meters
 */
double pdr_stride_length_kim(double avg_abs_accel, double K);

/* ============================================================================
 * L5 - Numerical Integration Methods for INS
 * ============================================================================ */

/**
 * @brief Update attitude using 1st-order quaternion integration
 *
 * q_{k+1} = q_k ⊗ [1, 0.5*w_x*dt, 0.5*w_y*dt, 0.5*w_z*dt]
 *
 * @param q Current quaternion (updated in place)
 * @param gx,gy,gz Angular rates in rad/s
 * @param dt Time step in seconds
 */
void quaternion_integrate_1st_order(quaternion_t *q,
                                   double gx, double gy, double gz, double dt);

/**
 * @brief Update attitude using 4th-order Runge-Kutta quaternion integration
 *
 * Higher accuracy than 1st-order for high-dynamic maneuvers.
 *
 * @param q Current quaternion (updated in place)
 * @param gx,gy,gz Angular rates in rad/s
 * @param dt Time step in seconds
 *
 * Complexity: O(1), 4 evaluations
 */
void quaternion_integrate_rk4(quaternion_t *q,
                              double gx, double gy, double gz, double dt);

/**
 * @brief Quaternion multiplication (Hamilton product): q1 ⊗ q2
 *
 * @param q1 First quaternion
 * @param q2 Second quaternion
 * @param[out] result q1 ⊗ q2
 */
void quaternion_multiply(const quaternion_t *q1, const quaternion_t *q2,
                         quaternion_t *result);

/**
 * @brief Compute the conjugate of a quaternion (inverse rotation)
 *
 * @param q Input quaternion
 * @param[out] result q* = (w, -x, -y, -z)
 */
void quaternion_conjugate(const quaternion_t *q, quaternion_t *result);

#endif /* INERTIAL_NAVIGATION_H */
