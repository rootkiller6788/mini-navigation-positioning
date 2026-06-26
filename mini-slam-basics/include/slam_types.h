#ifndef SLAM_TYPES_H
#define SLAM_TYPES_H

/**
 * @file    slam_types.h
 * @brief   SLAM core data structures: pose, landmark, map, covariance
 *
 * Simultaneous Localization and Mapping (SLAM) — the problem of building
 * a map of an unknown environment while simultaneously tracking the
 * robot's pose within that map.
 *
 * Reference: Durrant-Whyte & Bailey (2006), "Simultaneous Localization
 *            and Mapping: Part I", IEEE Robotics & Automation Magazine.
 *            Thrun, Burgard & Fox (2005), "Probabilistic Robotics", MIT Press.
 *
 * Knowledge Coverage:
 *   L1: Pose (SE(2)), Landmark, Map, Covariance, Odometry, Observation
 *   L2: SLAM probabilistic formulation, Online vs Full SLAM
 *   L3: SE(2) Lie algebra, Information matrix, Mahalanobis distance
 */

#include <stddef.h>
#include <stdint.h>
#include <math.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* =========================================================================
 * L1: Core Enumerations
 * ========================================================================= */

/** SLAM algorithm backend */
typedef enum {
    SLAM_BACKEND_EKF        = 0,
    SLAM_BACKEND_FASTSLAM   = 1,
    SLAM_BACKEND_GRAPH      = 2,
    SLAM_BACKEND_ISAM       = 3,
    SLAM_BACKEND_SEIF       = 4
} slam_backend_t;

/** Sensor type for observation models */
typedef enum {
    SLAM_SENSOR_RANGE_BEARING  = 0,
    SLAM_SENSOR_RANGE_ONLY     = 1,
    SLAM_SENSOR_BEARING_ONLY   = 2,
    SLAM_SENSOR_LIDAR_2D       = 3,
    SLAM_SENSOR_CAMERA_MONO    = 4,
    SLAM_SENSOR_CAMERA_STEREO  = 5,
    SLAM_SENSOR_RGBD           = 6,
    SLAM_SENSOR_SONAR          = 7
} slam_sensor_type_t;

/** Motion model type */
typedef enum {
    SLAM_MOTION_VELOCITY     = 0,
    SLAM_MOTION_ODOMETRY     = 1,
    SLAM_MOTION_DIFF_DRIVE   = 2,
    SLAM_MOTION_ACKERMANN    = 3
} slam_motion_type_t;

/** Data association method */
typedef enum {
    SLAM_DA_NEAREST_NEIGHBOR     = 0,
    SLAM_DA_MAHALANOBIS_GATE     = 1,
    SLAM_DA_JCBB                 = 2,
    SLAM_DA_ICP                   = 3,
    SLAM_DA_HUNGARIAN             = 4,
    SLAM_DA_GRAPH_MATCHING        = 5
} slam_da_method_t;

/** SLAM system operation status */
typedef enum {
    SLAM_STATUS_IDLE            = 0,
    SLAM_STATUS_INITIALIZING    = 1,
    SLAM_STATUS_RUNNING         = 2,
    SLAM_STATUS_LOOP_CLOSURE    = 3,
    SLAM_STATUS_RELOCALIZING    = 4,
    SLAM_STATUS_ERROR           = 5
} slam_status_t;

/** Landmark type classification */
typedef enum {
    SLAM_LM_POINT          = 0,
    SLAM_LM_LINE           = 1,
    SLAM_LM_PLANE          = 2,
    SLAM_LM_CYLINDER       = 3,
    SLAM_LM_SPHERE         = 4,
    SLAM_LM_ORB_FEATURE    = 5,
    SLAM_LM_SIFT_FEATURE   = 6
} slam_landmark_type_t;

/* =========================================================================
 * L1: Pose in SE(2) — 2D rigid-body transformation
 * ========================================================================= */

/**
 * @brief 2D robot pose in SE(2) Lie group
 *
 * SE(2) = R^2 × SO(2): (x, y, θ) where θ ∈ [-π, π)
 * The pose represents a rigid-body transformation from robot frame
 * to world frame:  p_world = R(θ) * p_robot + [x, y]^T
 *
 * Composition (⊕): pose_a ⊕ pose_b = pose_a ∘ pose_b
 *   x' = x_a + x_b*cos(θ_a) - y_b*sin(θ_a)
 *   y' = y_a + x_b*sin(θ_a) + y_b*cos(θ_a)
 *   θ' = normalize(θ_a + θ_b)
 */
typedef struct {
    double x;     /**< x translation [m] */
    double y;     /**< y translation [m] */
    double theta; /**< rotation angle [rad], normalized to [-π, π) */
} slam_pose2d_t;

/**
 * @brief 3D robot pose in SE(3): position + quaternion orientation
 *
 * SE(3) = R^3 × SO(3). Quaternion avoids gimbal lock.
 * Hamiltonian convention: q = [qw, qx, qy, qz], ‖q‖ = 1.
 */
typedef struct {
    double x, y, z;    /**< translation [m] */
    double qw, qx, qy, qz; /**< unit quaternion */
} slam_pose3d_t;

/**
 * @brief Velocity twist for differential drive robots
 *
 * Controls: translational velocity v [m/s], angular velocity ω [rad/s]
 * State propagation (discrete-time, Δt seconds):
 *   x_{t+1} = x_t + v·Δt·cos(θ_t + ω·Δt/2)
 *   y_{t+1} = y_t + v·Δt·sin(θ_t + ω·Δt/2)
 *   θ_{t+1} = θ_t + ω·Δt
 */
typedef struct {
    double v;    /**< linear velocity [m/s] */
    double omega; /**< angular velocity [rad/s] */
    double dt;   /**< time step [s] */
} slam_velocity_t;

/* =========================================================================
 * L1: Landmark
 * ========================================================================= */

/**
 * @brief 2D point landmark
 *
 * A landmark is a distinguishable feature in the environment whose
 * position is estimated jointly with the robot pose.
 *
 * In EKF-SLAM, each landmark adds 2 state dimensions (x, y).
 * The full state vector is: [pose(3) | lm1(2) | lm2(2) | ... | lmN(2)]
 *
 * Signature: unique descriptor for data association (e.g., SIFT/ORB vector).
 */
typedef struct {
    int32_t  id;         /**< unique landmark identifier */
    double   x;          /**< world-frame x position [m] */
    double   y;          /**< world-frame y position [m] */
    double   covariance[4]; /**< 2×2 covariance [σ_xx, σ_xy; σ_yx, σ_yy] */
    slam_landmark_type_t type;
    double   signature[32]; /**< feature descriptor (≤32 dims) */
    int32_t  signature_dim; /**< actual descriptor dimensionality */
    int32_t  observed_count; /**< number of times observed */
    bool     is_active;     /**< whether landmark is currently tracked */
} slam_landmark2d_t;

/**
 * @brief 3D point landmark with full covariance
 */
typedef struct {
    int32_t  id;
    double   x, y, z;     /**< world-frame position [m] */
    double   covariance[9]; /**< 3×3 covariance (row-major) */
    slam_landmark_type_t type;
    double   signature[64]; /**< 3D feature descriptor */
    int32_t  signature_dim;
    int32_t  observed_count;
    bool     is_active;
} slam_landmark3d_t;

/* =========================================================================
 * L1: Sensor Observation
 * ========================================================================= */

/**
 * @brief Range-bearing observation (2D LiDAR / sonar)
 *
 * Measurement model (sensor frame → world frame):
 *   z = [r, φ]^T
 *   h(x, m_j) = [ sqrt((m_jx − x)^2 + (m_jy − y)^2) − R,
 *                  atan2(m_jy − y, m_jx − x) − θ ]^T + noise
 *
 * where R is the sensor offset on the robot.
 * Noise: ε ∼ N(0, R_noise), R_noise = diag(σ_r², σ_φ²)
 */
typedef struct {
    double   range;        /**< range [m] */
    double   bearing;      /**< bearing angle [rad], relative to robot heading */
    int32_t  landmark_id;  /**< -1 if unknown (data association needed) */
    double   range_std;    /**< range measurement std dev [m] */
    double   bearing_std;  /**< bearing measurement std dev [rad] */
    uint64_t timestamp;    /**< measurement timestamp [μs] */
} slam_obs_rb_t;

/**
 * @brief 2D LiDAR scan — array of ranges at fixed angular increments
 */
typedef struct {
    double   angle_min;    /**< start angle [rad] */
    double   angle_max;    /**< end angle [rad] */
    double   angle_inc;    /**< angular increment [rad] */
    double   range_min;    /**< minimum valid range [m] */
    double   range_max;    /**< maximum valid range [m] */
    int32_t  num_ranges;   /**< number of range readings */
    double  *ranges;       /**< range readings [m], length = num_ranges */
    double  *intensities;  /**< optional intensities */
    uint64_t timestamp;    /**< scan timestamp [μs] */
    slam_pose2d_t sensor_pose; /**< sensor pose relative to robot frame */
} slam_lidar_scan_t;

/* =========================================================================
 * L1: Odometry / Control Input
 * ========================================================================= */

/**
 * @brief Odometry measurement between two time steps
 *
 * In wheeled robots, odometry measures the relative motion:
 *   δ_rot1 = atan2(y'−y, x'−x) − θ
 *   δ_trans = sqrt((x'−x)² + (y'−y)²)
 *   δ_rot2 = θ' − θ − δ_rot1
 *
 * Probabilistic motion model (Thrun p.136):
 *   δ_rot1_hat ∼ N(δ_rot1, α₁·|δ_rot1| + α₂·δ_trans)
 *   δ_trans_hat ∼ N(δ_trans, α₃·δ_trans + α₄·(|δ_rot1|+|δ_rot2|))
 *   δ_rot2_hat ∼ N(δ_rot2, α₁·|δ_rot2| + α₂·δ_trans)
 */
typedef struct {
    double   delta_rot1;  /**< first rotation [rad] */
    double   delta_trans; /**< translation [m] */
    double   delta_rot2;  /**< second rotation [rad] */
    double   alpha[4];    /**< noise parameters α₁,α₂,α₃,α₄ */
    uint64_t timestamp;
} slam_odometry_t;

/**
 * @brief Control input for velocity motion model
 */
typedef struct {
    double   v_cmd;       /**< commanded linear velocity [m/s] */
    double   omega_cmd;   /**< commanded angular velocity [rad/s] */
    double   alpha[6];    /**< noise params α₁..α₆ (Thrun p.128) */
    uint64_t timestamp;
} slam_control_t;

/* =========================================================================
 * L1: SLAM Map
 * ========================================================================= */

/**
 * @brief 2D feature-based map
 *
 * Stores all estimated landmark positions. In EKF-SLAM the map is
 * tightly coupled with the state vector. In FastSLAM/Graph SLAM
 * the map is maintained separately.
 */
typedef struct {
    int32_t             capacity;    /**< max number of landmarks */
    int32_t             count;       /**< current landmark count */
    slam_landmark2d_t  *landmarks;   /**< dynamic array of landmarks */
    uint64_t            last_update; /**< timestamp of last map change */
} slam_map2d_t;

/**
 * @brief Occupancy grid map (2D)
 *
 * Each cell stores log-odds of occupancy: l(x) = log(p(x)/(1−p(x)))
 * p(x) = 1 − 1/(1 + exp(l(x)))
 *
 * Sensor model: inverse sensor model for laser range finder
 * (Thrun, Chapter 9).
 */
typedef struct {
    double   origin_x;    /**< map origin x [m] */
    double   origin_y;    /**< map origin y [m] */
    double   resolution;  /**< cell side length [m] */
    int32_t  width;       /**< cells in x direction */
    int32_t  height;      /**< cells in y direction */
    double  *log_odds;    /**< log-odds values, size = width*height */
    double   lo_occ;      /**< log-odds for occupied cell */
    double   lo_free;     /**< log-odds for free cell */
    double   lo_min;      /**< clamp minimum */
    double   lo_max;      /**< clamp maximum */
} slam_occgrid_t;

/* =========================================================================
 * L1: SLAM System State
 * ========================================================================= */

/**
 * @brief Full SLAM system state (EKF-based)
 *
 * State vector layout for N landmarks:
 *   μ = [x, y, θ, m1x, m1y, m2x, m2y, ..., mNx, mNy]^T
 *   dimension = 3 + 2*N
 *
 * Covariance matrix Σ: (3+2N) × (3+2N) — dense for EKF-SLAM.
 * Key sub-blocks:
 *   Σ_rr (3×3): robot pose covariance
 *   Σ_rl (3×2N): robot-landmark cross-covariance
 *   Σ_ll (2N×2N): landmark-landmark covariance
 */
typedef struct {
    slam_pose2d_t  robot_pose;    /**< current robot pose estimate */
    int32_t        num_landmarks; /**< number of mapped landmarks */
    int32_t        state_dim;     /**< total state dimension = 3 + 2N */
    int32_t        cov_stride;    /**< allocated covariance row stride (same for all) */
    double        *state_mean;    /**< mean vector μ, length = state_dim */
    double        *covariance;    /**< covariance Σ, cov_stride × cov_stride, row-major */
    uint64_t       timestamp;     /**< last update timestamp */
    slam_status_t  status;
} slam_ekf_state_t;

/**
 * @brief FastSLAM particle
 *
 * Each particle maintains its own robot trajectory estimate and
 * a set of independent EKFs (one per landmark). This decouples
 * the landmark estimation problem.
 */
typedef struct {
    slam_pose2d_t   pose;          /**< robot pose sample */
    double          weight;        /**< importance weight */
    double          log_weight;    /**< log importance weight (numeric stability) */
    int32_t         num_landmarks; /**< landmarks in this particle's map */
    int32_t         lm_capacity;   /**< allocated landmark capacity */
    slam_landmark2d_t *landmarks;  /**< per-particle landmark estimates */
    double         *lm_covariances; /**< per-landmark covariances (num_lm × 4) */
    double           trajectory_x[1024]; /**< recent trajectory x history */
    double           trajectory_y[1024]; /**< recent trajectory y history */
    int32_t          traj_len;     /**< trajectory history length */
} slam_particle_t;

/**
 * @brief Pose graph node and edge for graph-based SLAM
 *
 * Graph SLAM formulates SLAM as a nonlinear least-squares problem
 * over a graph G = (V, E):
 *   V: robot poses (nodes), E: constraints between poses (edges)
 *
 * Objective: X* = argmin_X Σ e_{ij}^T Ω_{ij} e_{ij}
 * where e_{ij} = z_{ij} − h(x_i, x_j) is the constraint error
 */
typedef struct {
    int32_t         id;            /**< unique node id */
    slam_pose2d_t   pose;          /**< node pose estimate */
    bool            is_fixed;      /**< anchor node (fixed during optimization) */
} slam_graph_node_t;

typedef struct {
    int32_t         id_a;          /**< source node id */
    int32_t         id_b;          /**< target node id */
    slam_pose2d_t   constraint;    /**< relative pose constraint (measurement) */
    double          info_matrix[9]; /**< 3×3 information matrix (row-major) */
    bool            is_loop_closure; /**< true if this edge closes a loop */
    int32_t         inlier_count;  /**< number of inlier correspondences */
} slam_graph_edge_t;

/**
 * @brief Pose graph for graph-based SLAM
 */
typedef struct {
    int32_t             num_nodes;
    int32_t             node_capacity;
    slam_graph_node_t  *nodes;
    int32_t             num_edges;
    int32_t             edge_capacity;
    slam_graph_edge_t  *edges;
    double             *hessian;       /**< H = J^T Ω J, size (3N)×(3N) */
    double             *gradient;      /**< b = J^T Ω e, size (3N) */
    double             *chi2_history;  /**< optimization convergence history */
    int32_t             chi2_len;
} slam_pose_graph_t;

/* =========================================================================
 * L1: SLAM Configuration
 * ========================================================================= */

/**
 * @brief SLAM system configuration parameters
 */
typedef struct {
    slam_backend_t      backend;
    slam_sensor_type_t   sensor_type;
    slam_motion_type_t   motion_type;
    slam_da_method_t     da_method;

    /* EKF/FastSLAM parameters */
    int32_t   max_landmarks;
    int32_t   max_particles;       /**< FastSLAM: number of particles */
    double    sigma_r;             /**< range measurement noise std [m] */
    double    sigma_b;             /**< bearing measurement noise std [rad] */
    double    sigma_v;             /**< linear velocity noise std [m/s] */
    double    sigma_omega;         /**< angular velocity noise std [rad/s] */

    /* Data association */
    double    mahalanobis_gate;    /**< chi-squared gating threshold (e.g., 5.991 for 95%) */
    double    max_range;           /**< maximum sensor range [m] */

    /* Graph optimization */
    int32_t   max_graph_iterations; /**< GN/LM max iterations */
    double    convergence_thresh;  /**< chi2 change threshold */
    double    lm_lambda_init;      /**< initial LM damping factor */

    /* Loop closure */
    double    loop_closure_radius; /**< search radius for loop candidates [m] */
    int32_t   loop_closure_min_steps; /**< min steps before loop check */

    /* Low-level */
    bool      use_robust_kernel;   /**< enable Huber/Cauchy robust kernel */
    double    robust_kernel_delta; /**< robust kernel threshold */
} slam_config_t;

/* =========================================================================
 * L2: Core SLAM Metrics
 * ========================================================================= */

/**
 * @brief SLAM performance metrics
 *
 * ATE: Absolute Trajectory Error — RMSE between estimated and ground-truth poses
 * RPE: Relative Pose Error — drift per unit distance
 * Map quality: landmark position error statistics
 */
typedef struct {
    double   ate_rmse;         /**< absolute trajectory RMSE [m] */
    double   ate_mean;         /**< absolute trajectory mean error [m] */
    double   ate_max;          /**< absolute trajectory max error [m] */
    double   rpe_trans;        /**< relative translational error [m/m] */
    double   rpe_rot;          /**< relative rotational error [rad/m] */
    double   map_rmse;         /**< landmark position RMSE [m] */
    int32_t  num_loop_closures; /**< total loop closures detected */
    double   loop_closure_error; /**< mean loop closure constraint error */
    double   cpu_time_per_step; /**< computation time per SLAM step [ms] */
} slam_metrics_t;

/* =========================================================================
 * L2: SLAM Core API — Return codes
 * ========================================================================= */

/** SLAM function return codes */
typedef enum {
    SLAM_OK                = 0,
    SLAM_ERR_NULL_PTR      = -1,
    SLAM_ERR_DIM_MISMATCH  = -2,
    SLAM_ERR_SINGULAR      = -3,
    SLAM_ERR_NOT_INIT      = -4,
    SLAM_ERR_NO_ASSOC      = -5,
    SLAM_ERR_MEMORY        = -6,
    SLAM_ERR_DIVERGED      = -7,
    SLAM_ERR_INVALID_PARAM = -8
} slam_error_t;

#ifdef __cplusplus
}
#endif

#endif /* SLAM_TYPES_H */
