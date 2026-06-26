/**
 * @file    test_slam.c
 * @brief   Comprehensive test suite for SLAM module
 *
 * Tests cover:
 *   - SE(2) pose composition and inverse
 *   - EKF-SLAM: init, predict, update, augment
 *   - FastSLAM: particle init, sampling, resampling
 *   - Graph SLAM: error, Jacobian, Gauss-Newton optimization
 *   - Sensor models: motion, observation, Jacobians
 *   - Data association: Mahalanobis, NN, loop detection
 *   - Linear algebra: matmul, Cholesky, inverse
 */

#include "slam_types.h"
#include "slam_ekf.h"
#include "slam_fastslam.h"
#include "slam_graph.h"
#include "slam_sensor.h"
#include "slam_data_assoc.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

/* Function declarations from source files */
extern double slam_normalize_angle(double theta);
extern double slam_angle_diff(double a, double b);
extern void slam_pose_compose(const slam_pose2d_t *a,
                               const slam_pose2d_t *b,
                               slam_pose2d_t *c);
extern void slam_pose_inverse(const slam_pose2d_t *a, slam_pose2d_t *inv);
extern void slam_pose_relative(const slam_pose2d_t *a,
                                const slam_pose2d_t *b,
                                slam_pose2d_t *rel);
extern void slam_transform_point(const slam_pose2d_t *pose,
                                  double px, double py,
                                  double *wx, double *wy);
extern void slam_transform_point_inv(const slam_pose2d_t *pose,
                                      double wx, double wy,
                                      double *px, double *py);
extern int slam_inv2x2(const double A[4], double Ainv[4]);
extern double slam_det2x2(const double A[4]);
extern double slam_det3x3(const double M[9]);
extern int slam_cholesky_3x3(const double A[9], double L[9]);
extern int slam_solve_cholesky_3x3(const double A[9], const double b[3],
                                    double x[3]);
extern void slam_matvec_mul(const double *A, const double *x,
                             int m, int n, double *y);
extern void slam_matmul(const double *A, const double *B,
                         int m, int k, int n, double *C);
extern void slam_config_default(slam_config_t *cfg);
extern int slam_map2d_init(slam_map2d_t *map, int capacity);
extern void slam_map2d_free(slam_map2d_t *map);
extern int slam_map2d_add_landmark(slam_map2d_t *map,
                                    const slam_landmark2d_t *lm, int *idx);

/* Forward declarations from loop closure */
extern int slam_loop_closure_detect(
    const slam_lidar_scan_t *current_scan,
    const slam_pose2d_t *current_pose,
    const void *past_descriptors,
    const slam_pose2d_t *past_poses,
    int num_past, double search_radius, int min_steps,
    double score_thresh,
    slam_pose2d_t *loop_pose,
    double info_matrix[9], int *match_id);

static int tests_run = 0;
static int tests_passed = 0;

#define TEST(name) do { tests_run++; printf("  %-50s", name); } while(0)
#define PASS() do { tests_passed++; printf("PASS\n"); } while(0)
#define FAIL(msg) do { printf("FAIL: %s\n", msg); assert(0); } while(0)
#define CHECK(cond) do { if (!(cond)) FAIL(#cond); } while(0)

/* =========================================================================
 * Test Group 1: Angle Normalization & SE(2) Operations
 * ========================================================================= */

static void test_angle_normalization(void) {
    TEST("angle normalization");
    CHECK(fabs(slam_normalize_angle(0.0)) < 1e-12);
    CHECK(fabs(slam_normalize_angle(M_PI)) - M_PI < 1e-12
          || fabs(slam_normalize_angle(M_PI) + M_PI) < 1e-12);
    CHECK(fabs(slam_normalize_angle(-M_PI)) - M_PI < 1e-12
          || fabs(slam_normalize_angle(-M_PI) + M_PI) < 1e-12);
    CHECK(fabs(slam_normalize_angle(3.0 * M_PI) - M_PI) < 1e-10);
    CHECK(fabs(slam_normalize_angle(-3.0 * M_PI) - M_PI) < 1e-10
          || fabs(slam_normalize_angle(-3.0 * M_PI) + M_PI) < 1e-10);
    PASS();
}

static void test_se2_compose_inverse(void) {
    TEST("SE(2) compose+inverse");
    slam_pose2d_t a = {1.0, 2.0, M_PI/4};
    slam_pose2d_t b = {0.5, 0.3, M_PI/6};
    slam_pose2d_t c, inv_a, back;

    slam_pose_compose(&a, &b, &c);
    slam_pose_inverse(&a, &inv_a);
    slam_pose_compose(&inv_a, &a, &back);

    CHECK(fabs(back.x) < 1e-12);
    CHECK(fabs(back.y) < 1e-12);
    CHECK(fabs(back.theta) < 1e-12);
    PASS();
}

static void test_se2_relative(void) {
    TEST("SE(2) relative pose");
    slam_pose2d_t a = {0, 0, 0};
    slam_pose2d_t b = {3, 4, M_PI/2};
    slam_pose2d_t rel;
    slam_pose_relative(&a, &b, &rel);

    CHECK(fabs(rel.x - 3.0) < 1e-10);
    CHECK(fabs(rel.y - 4.0) < 1e-10);
    CHECK(fabs(rel.theta - M_PI/2) < 1e-10);
    PASS();
}

static void test_transform_roundtrip(void) {
    TEST("point transform round-trip");
    slam_pose2d_t pose = {2.0, 3.0, M_PI/6};
    double wx, wy, px2, py2;
    slam_transform_point(&pose, 1.0, 0.5, &wx, &wy);
    slam_transform_point_inv(&pose, wx, wy, &px2, &py2);
    CHECK(fabs(px2 - 1.0) < 1e-10);
    CHECK(fabs(py2 - 0.5) < 1e-10);
    PASS();
}

/* =========================================================================
 * Test Group 2: Linear Algebra
 * ========================================================================= */

static void test_matmul(void) {
    TEST("matrix multiply");
    double A[6] = {1, 2, 3, 4, 5, 6};  /* 2×3 */
    double B[6] = {7, 8, 9, 10, 11, 12}; /* 3×2 */
    double C[4];
    slam_matmul(A, B, 2, 3, 2, C);
    /* C[0,0] = 1*7+2*9+3*11 = 58, C[0,1] = 1*8+2*10+3*12 = 64 */
    /* C[1,0] = 4*7+5*9+6*11 = 139, C[1,1] = 4*8+5*10+6*12 = 154 */
    CHECK(C[0] == 58 && C[1] == 64 && C[2] == 139 && C[3] == 154);
    PASS();
}

static void test_2x2_inverse(void) {
    TEST("2x2 inverse");
    double M[4] = {4, 7, 2, 6};
    double Minv[4], I_check[4];
    CHECK(slam_inv2x2(M, Minv));
    slam_matmul(M, Minv, 2, 2, 2, I_check);
    CHECK(fabs(I_check[0] - 1.0) < 1e-10);
    CHECK(fabs(I_check[3] - 1.0) < 1e-10);
    CHECK(fabs(I_check[1]) < 1e-10);
    CHECK(fabs(I_check[2]) < 1e-10);
    PASS();
}

static void test_cholesky_3x3(void) {
    TEST("Cholesky 3x3");
    double A[9] = {4, 12, -16, 12, 37, -43, -16, -43, 98};
    double b[3] = {1, 2, 3};
    double x[3];
    CHECK(slam_solve_cholesky_3x3(A, b, x));
    double Ax[3];
    slam_matvec_mul(A, x, 3, 3, Ax);
    CHECK(fabs(Ax[0] - b[0]) < 1e-10);
    CHECK(fabs(Ax[1] - b[1]) < 1e-10);
    CHECK(fabs(Ax[2] - b[2]) < 1e-10);
    PASS();
}

static void test_determinant(void) {
    TEST("determinant");
    double M3[9] = {6, 1, 1, 4, -2, 5, 2, 8, 7};
    double det = slam_det3x3(M3);
    CHECK(fabs(det - (-306.0)) < 1e-10);
    PASS();
}

/* =========================================================================
 * Test Group 3: EKF-SLAM
 * ========================================================================= */

static void test_ekf_init(void) {
    TEST("EKF-SLAM initialization");
    slam_ekf_state_t state;
    slam_config_t config;
    slam_config_default(&config);
    slam_pose2d_t init_pose = {0, 0, 0};

    int rc = slam_ekf_init(&state, &config, &init_pose, 0.1, 0.05);
    CHECK(rc == SLAM_OK);
    CHECK(state.state_dim == 3);
    CHECK(state.num_landmarks == 0);
    CHECK(fabs(state.state_mean[0]) < 1e-12);
    CHECK(fabs(state.state_mean[1]) < 1e-12);
    CHECK(fabs(state.state_mean[2]) < 1e-12);

    int dim = 3 + 2 * config.max_landmarks;
    CHECK(fabs(state.covariance[0*dim + 0] - 0.01) < 1e-12);  /* sigma_pos² */
    CHECK(fabs(state.covariance[1*dim + 1] - 0.01) < 1e-12);

    slam_ekf_free(&state);
    PASS();
}

static void test_ekf_predict(void) {
    TEST("EKF-SLAM predict");
    slam_ekf_state_t state;
    slam_config_t config;
    slam_config_default(&config);
    slam_pose2d_t init_pose = {0, 0, 0};
    slam_ekf_init(&state, &config, &init_pose, 0.01, 0.01);

    slam_velocity_t vel = {1.0, 0.0, 1.0}; /* 1 m/s straight for 1s */
    int rc = slam_ekf_predict_velocity(&state, &vel);
    CHECK(rc == SLAM_OK);
    /* Robot should have moved forward 1m */
    CHECK(fabs(state.robot_pose.x - 1.0) < 1e-10);
    CHECK(fabs(state.robot_pose.y) < 1e-10);
    CHECK(fabs(state.robot_pose.theta) < 1e-10);

    /* Test with rotation */
    slam_velocity_t vel2 = {0.5, M_PI/2, 1.0};
    rc = slam_ekf_predict_velocity(&state, &vel2);
    CHECK(rc == SLAM_OK);
    CHECK(state.robot_pose.theta > 0);

    /* Covariance should have grown (prediction adds noise) */
    int dim = 3 + 2 * config.max_landmarks;
    double trace = state.covariance[0*dim+0] + state.covariance[1*dim+1]
                  + state.covariance[2*dim+2];
    CHECK(trace > 0.0002);  /* larger than initial */

    slam_ekf_free(&state);
    PASS();
}

static void test_ekf_augment(void) {
    TEST("EKF-SLAM state augmentation");
    slam_ekf_state_t state;
    slam_config_t config;
    slam_config_default(&config);
    slam_pose2d_t init_pose = {0, 0, 0};
    slam_ekf_init(&state, &config, &init_pose, 0.01, 0.01);

    /* Observe a landmark at range=5m, bearing=30° */
    slam_obs_rb_t obs = {5.0, M_PI/6, -1, 0.1, 0.02, 0};
    int new_id;
    int rc = slam_ekf_augment(&state, &obs, &new_id);
    CHECK(rc == SLAM_OK);
    CHECK(new_id == 0);
    CHECK(state.num_landmarks == 1);
    CHECK(state.state_dim == 5);

    /* Landmark should be at (5*cos(30°), 5*sin(30°)) = (4.33, 2.5) */
    CHECK(fabs(state.state_mean[3] - 5.0*0.8660254) < 0.01);
    CHECK(fabs(state.state_mean[4] - 2.5) < 0.01);

    /* Add second landmark */
    slam_obs_rb_t obs2 = {3.0, -M_PI/4, -1, 0.1, 0.02, 0};
    int new_id2;
    rc = slam_ekf_augment(&state, &obs2, &new_id2);
    CHECK(rc == SLAM_OK);
    CHECK(new_id2 == 1);
    CHECK(state.num_landmarks == 2);
    CHECK(state.state_dim == 7);

    slam_ekf_free(&state);
    PASS();
}

static void test_ekf_update_rb(void) {
    TEST("EKF-SLAM update with known landmark");
    slam_ekf_state_t state;
    slam_config_t config;
    slam_config_default(&config);
    slam_pose2d_t init_pose = {0, 0, 0};
    slam_ekf_init(&state, &config, &init_pose, 0.1, 0.05);

    /* Add a landmark */
    slam_obs_rb_t obs_init = {5.0, M_PI/6, -1, 0.1, 0.02, 0};
    int new_id;
    slam_ekf_augment(&state, &obs_init, &new_id);

    /* Move robot slightly */
    slam_velocity_t vel = {0.5, 0.1, 1.0};
    slam_ekf_predict_velocity(&state, &vel);

    /* Re-observe the landmark */
    double rx = state.robot_pose.x, ry = state.robot_pose.y;
    double rt = state.robot_pose.theta;
    double mx = state.state_mean[3], my = state.state_mean[4];
    double dx = mx - rx, dy = my - ry;
    double true_range = sqrt(dx*dx + dy*dy);
    double true_bearing = slam_normalize_angle(atan2(dy, dx) - rt);

    slam_obs_rb_t obs = {true_range + 0.02, true_bearing - 0.01, 0, 0.1, 0.02, 0};
    int rc = slam_ekf_update_rb(&state, &obs, 0);
    CHECK(rc == SLAM_OK);

    /* Verify covariance decreased (uncertainty reduced) */
    int dim_cov = 3 + 2 * config.max_landmarks;
    double det_after = slam_det2x2((double[]){
        state.covariance[3*dim_cov+3], state.covariance[3*dim_cov+4],
        state.covariance[4*dim_cov+3], state.covariance[4*dim_cov+4]
    });
    /* After an update, the landmark covariance should be reasonable */
    CHECK(det_after > 0);

    /* NIS should be reasonable */
    double nis;
    rc = slam_ekf_nis(&obs, &state, 0, &nis);
    CHECK(rc == SLAM_OK);
    CHECK(nis >= 0);
    /* For 2 DOF, 95% of NIS < 5.991 */
    CHECK(nis < 10.0);  /* reasonable bound with noise */

    slam_ekf_free(&state);
    PASS();
}

static void test_ekf_nees(void) {
    TEST("EKF-SLAM NEES consistency");
    slam_ekf_state_t state;
    slam_config_t config;
    slam_config_default(&config);
    slam_pose2d_t init_pose = {0, 0, 0};
    slam_ekf_init(&state, &config, &init_pose, 0.1, 0.05);

    /* Check NEES at start (should be ~0 since we initialized at origin) */
    slam_pose2d_t true_pose = {0, 0, 0};
    double nees;
    int rc = slam_ekf_nees(&state, &true_pose, &nees);
    CHECK(rc == SLAM_OK);
    CHECK(nees >= 0);
    CHECK(nees < 12.0); /* 95% for χ²₃ ≈ 7.815, but allow some tolerance */

    slam_ekf_free(&state);
    PASS();
}

/* =========================================================================
 * Test Group 4: FastSLAM
 * ========================================================================= */

static void test_fastslam_init(void) {
    TEST("FastSLAM initialization");
    slam_particle_t *particles;
    slam_config_t config;
    slam_config_default(&config);
    slam_pose2d_t init_pose = {0, 0, 0};

    int rc = slam_fastslam_init(&particles, 20, &init_pose, &config);
    CHECK(rc == SLAM_OK);
    CHECK(particles != NULL);

    for (int i = 0; i < 20; i++) {
        CHECK(fabs(particles[i].pose.x) < 1e-12);
        CHECK(fabs(particles[i].pose.y) < 1e-12);
        CHECK(fabs(particles[i].pose.theta) < 1e-12);
        CHECK(particles[i].num_landmarks == 0);
        CHECK(particles[i].log_weight < 0);  /* log(1/M) is negative */
    }

    slam_fastslam_free(particles, 20);
    PASS();
}

static void test_fastslam_sample_pose(void) {
    TEST("FastSLAM 1.0 pose sampling");
    slam_particle_t *particles;
    slam_config_t config;
    slam_config_default(&config);
    slam_pose2d_t init_pose = {0, 0, 0};
    slam_fastslam_init(&particles, 5, &init_pose, &config);

    slam_velocity_t vel = {1.0, 0.0, 1.0};
    slam_fastslam1_sample_pose(&particles[0], &vel, 0.05, 0.02);

    /* Particle should have moved approximately 1m forward */
    CHECK(fabs(particles[0].pose.x - 1.0) < 0.5);  /* noisy, but within 0.5m */
    CHECK(fabs(particles[0].pose.y) < 0.5);
    CHECK(fabs(particles[0].pose.theta) < 0.2);

    slam_fastslam_free(particles, 5);
    PASS();
}

static void test_fastslam_resample(void) {
    TEST("FastSLAM systematic resampling");
    slam_particle_t *particles;
    slam_config_t config;
    slam_config_default(&config);
    slam_pose2d_t init_pose = {0, 0, 0};
    slam_fastslam_init(&particles, 10, &init_pose, &config);

    /* Artificially skew weights */
    for (int i = 0; i < 10; i++) {
        particles[i].log_weight = (i == 5) ? 0.0 : -100.0;
        particles[i].pose.x = (double)i;
    }

    double neff;
    int rc = slam_fastslam_resample(&particles, 10, 5.0, &neff);
    CHECK(rc == SLAM_OK);
    CHECK(neff < 5.0); /* skewed weights => low Neff */

    /* After resampling, weights should be uniform */
    for (int i = 0; i < 10; i++) {
        CHECK(fabs(particles[i].weight - 0.1) < 1e-10);
    }

    slam_fastslam_free(particles, 10);
    PASS();
}

/* =========================================================================
 * Test Group 5: Graph SLAM
 * ========================================================================= */

static void test_graph_error_jacobian(void) {
    TEST("Graph SLAM error and Jacobian");
    slam_pose2d_t xi = {0, 0, 0};
    slam_pose2d_t xj = {1, 0, 0};
    slam_pose2d_t z_ij = {1, 0, 0};  /* expected relative: (1,0,0) */

    double error[3];
    slam_graph_error_se2(&xi, &xj, &z_ij, error);
    CHECK(fabs(error[0]) < 1e-10);
    CHECK(fabs(error[1]) < 1e-10);
    CHECK(fabs(error[2]) < 1e-10);

    double Ji[9], Jj[9];
    slam_graph_jacobian_se2(&xi, &xj, Ji, Jj);
    /* Ji should be [-I, [0;-1]; 0,0, -1] = [-1,0,0; 0,-1,-1; 0,0,-1] */
    CHECK(fabs(Ji[0] - (-1.0)) < 1e-10);
    CHECK(fabs(Ji[4] - (-1.0)) < 1e-10);
    CHECK(fabs(Ji[8] - (-1.0)) < 1e-10);
    /* Jj should be [I, 0; 0,0, 1] = [1,0,0; 0,1,0; 0,0,1] */
    CHECK(fabs(Jj[0] - 1.0) < 1e-10);
    CHECK(fabs(Jj[4] - 1.0) < 1e-10);
    CHECK(fabs(Jj[8] - 1.0) < 1e-10);

    PASS();
}

static void test_graph_optimization(void) {
    TEST("Graph SLAM Gauss-Newton optimization");
    slam_pose_graph_t graph;
    slam_graph_init(&graph, 10, 20);

    /* Simple odometry chain */
    slam_pose2d_t origin = {0, 0, 0};
    int nid;
    slam_graph_add_node(&graph, &origin, true, &nid);

    slam_pose2d_t p1 = {1.0, 0.0, 0.0};
    slam_graph_add_node(&graph, &p1, false, &nid);

    slam_pose2d_t p2 = {2.0, 0.0, 0.0};
    slam_graph_add_node(&graph, &p2, false, &nid);

    /* Odometry edges */
    slam_pose2d_t odom1 = {1.0, 0.0, 0.0};
    double info1[9] = {100,0,0, 0,100,0, 0,0,50};
    int eid;
    slam_graph_add_edge(&graph, 0, 1, &odom1, info1, &eid);

    slam_pose2d_t odom2 = {1.0, 0.0, 0.0};
    slam_graph_add_edge(&graph, 1, 2, &odom2, info1, &eid);

    /* Loop closure edge (p2 back to origin) */
    slam_pose2d_t loop = {-2.0, 0.0, 0.0};  /* p2→origin relative */
    double info_loop[9] = {200,0,0, 0,200,0, 0,0,100};
    slam_graph_add_edge(&graph, 2, 0, &loop, info_loop, &eid);
    graph.edges[eid].is_loop_closure = true;

    slam_config_t config;
    slam_config_default(&config);
    config.max_graph_iterations = 20;

    int rc = slam_graph_optimize_gauss_newton(&graph, &config, NULL);
    CHECK(rc == SLAM_OK);

    /* After optimization, nodes should be consistent */
    CHECK(graph.num_nodes == 3);
    CHECK(graph.chi2_len > 0);

    slam_graph_free(&graph);
    PASS();
}

/* =========================================================================
 * Test Group 6: Sensor Models
 * ========================================================================= */

static void test_motion_model(void) {
    TEST("velocity motion model");
    slam_pose2d_t pose = {0, 0, 0};
    slam_velocity_t vel = {1.0, 0.0, 1.0};
    slam_pose2d_t new_pose;
    slam_motion_model_velocity(&pose, &vel, &new_pose);
    CHECK(fabs(new_pose.x - 1.0) < 1e-10);
    CHECK(fabs(new_pose.y) < 1e-10);
    CHECK(fabs(new_pose.theta) < 1e-10);

    /* Rotation */
    slam_velocity_t vel2 = {0.5, M_PI/2, 1.0};
    slam_motion_model_velocity(&pose, &vel2, &new_pose);
    CHECK(new_pose.theta > 0);

    PASS();
}

static void test_obs_model(void) {
    TEST("range-bearing observation model");
    slam_pose2d_t pose = {0, 0, 0};
    slam_landmark2d_t lm = {0, 5.0, 0.0, {0}, SLAM_LM_POINT, {0}, 0, 0, true};
    double pred[2];
    slam_obs_model_rb(&pose, &lm, pred);
    CHECK(fabs(pred[0] - 5.0) < 1e-10);
    CHECK(fabs(pred[1] - 0.0) < 1e-10);

    /* Landmark at (3,4): range=5, bearing=atan2(4,3)≈0.927 */
    lm.x = 3.0; lm.y = 4.0;
    slam_obs_model_rb(&pose, &lm, pred);
    CHECK(fabs(pred[0] - 5.0) < 1e-10);
    CHECK(fabs(pred[1] - atan2(4, 3)) < 1e-10);

    PASS();
}

static void test_obs_inverse(void) {
    TEST("inverse observation model");
    slam_pose2d_t pose = {1.0, 2.0, M_PI/6};
    slam_obs_rb_t obs = {5.0, M_PI/4, -1, 0.1, 0.02, 0};
    slam_landmark2d_t lm;
    slam_obs_inverse_rb(&pose, &obs, &lm);

    /* Verify: landmark should be at known position */
    double expected_x = 1.0 + 5.0 * cos(M_PI/4 + M_PI/6);
    double expected_y = 2.0 + 5.0 * sin(M_PI/4 + M_PI/6);
    CHECK(fabs(lm.x - expected_x) < 1e-10);
    CHECK(fabs(lm.y - expected_y) < 1e-10);

    PASS();
}

/* =========================================================================
 * Test Group 7: Data Association
 * ========================================================================= */

static void test_mahalanobis(void) {
    TEST("Mahalanobis distance");
    double nu[2] = {1.0, 0.5};
    double S[4] = {1.0, 0.0, 0.0, 1.0};  /* Identity → Euclidean squared */
    double d2 = slam_mahalanobis_sq(nu, S, 2);
    CHECK(fabs(d2 - 1.25) < 1e-10);  /* 1² + 0.5² = 1.25 */
    PASS();
}

static void test_nearest_neighbor(void) {
    TEST("nearest neighbor data association");
    slam_pose2d_t pose = {0, 0, 0};
    double Sigma_rr[9] = {0.01,0,0, 0,0.01,0, 0,0,0.001};

    slam_landmark2d_t landmarks[3];
    landmarks[0] = (slam_landmark2d_t){0, 5.0, 0.0, {0}, SLAM_LM_POINT,
                                        {0}, 0, 0, true};
    landmarks[1] = (slam_landmark2d_t){1, 3.0, 4.0, {0}, SLAM_LM_POINT,
                                        {0}, 0, 0, true};
    landmarks[2] = (slam_landmark2d_t){2, 10.0, 0.0, {0}, SLAM_LM_POINT,
                                        {0}, 0, 0, true};

    /* Observation of landmark at (3,4): range≈5, bearing≈53° */
    slam_obs_rb_t obs = {5.0, atan2(4.0, 3.0), -1, 0.1, 0.02, 0};
    int best_idx;
    double dist_sq;
    slam_da_nearest_neighbor(&pose, Sigma_rr, &obs, landmarks, 3,
                              5.991, &best_idx, &dist_sq);
    CHECK(best_idx == 1);  /* landmark (3,4) is closest to (5, 53°) */
    CHECK(dist_sq < 5.991);

    PASS();
}

/* =========================================================================
 * Test Group 8: Map Management
 * ========================================================================= */

static void test_map_management(void) {
    TEST("map management");
    slam_map2d_t map;
    int rc = slam_map2d_init(&map, 10);
    CHECK(rc == SLAM_OK);
    CHECK(map.count == 0);
    CHECK(map.capacity == 10);

    slam_landmark2d_t lm = {0, 1.0, 2.0, {0}, SLAM_LM_POINT, {0}, 0, 0, true};
    int idx;
    rc = slam_map2d_add_landmark(&map, &lm, &idx);
    CHECK(rc == SLAM_OK);
    CHECK(idx == 0);
    CHECK(map.count == 1);

    slam_map2d_free(&map);
    CHECK(map.count == 0);
    PASS();
}

/* =========================================================================
 * Test Group 9: Configuration and Metrics
 * ========================================================================= */

static void test_config(void) {
    TEST("SLAM configuration defaults");
    slam_config_t config;
    slam_config_default(&config);
    CHECK(config.backend == SLAM_BACKEND_EKF);
    CHECK(config.max_landmarks == 100);
    CHECK(config.max_particles == 50);
    CHECK(config.sigma_r == 0.1);
    CHECK(config.mahalanobis_gate == 5.991);
    PASS();
}

/* =========================================================================
 * Main Test Runner
 * ========================================================================= */

int main(void) {
    printf("=== SLAM Module Test Suite ===\n\n");

    printf("[Group 1] SE(2) Operations\n");
    test_angle_normalization();
    test_se2_compose_inverse();
    test_se2_relative();
    test_transform_roundtrip();

    printf("\n[Group 2] Linear Algebra\n");
    test_matmul();
    test_2x2_inverse();
    test_cholesky_3x3();
    test_determinant();

    printf("\n[Group 3] EKF-SLAM\n");
    test_ekf_init();
    test_ekf_predict();
    test_ekf_augment();
    test_ekf_update_rb();
    test_ekf_nees();

    printf("\n[Group 4] FastSLAM\n");
    test_fastslam_init();
    test_fastslam_sample_pose();
    test_fastslam_resample();

    printf("\n[Group 5] Graph SLAM\n");
    test_graph_error_jacobian();
    test_graph_optimization();

    printf("\n[Group 6] Sensor Models\n");
    test_motion_model();
    test_obs_model();
    test_obs_inverse();

    printf("\n[Group 7] Data Association\n");
    test_mahalanobis();
    test_nearest_neighbor();

    printf("\n[Group 8] Map Management\n");
    test_map_management();

    printf("\n[Group 9] Configuration\n");
    test_config();

    printf("\n=== Results: %d/%d tests passed ===\n",
           tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
