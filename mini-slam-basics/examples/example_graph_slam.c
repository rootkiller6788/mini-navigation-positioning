/**
 * @file    example_graph_slam.c
 * @brief   Pose Graph SLAM: Loop Closure Demo
 *
 * Demonstrates graph-based SLAM with a robot driving a rectangular loop
 * and then detecting a loop closure when returning to the start.
 *
 * Scenario:
 *   - Robot drives: (0,0)→(10,0)→(10,10)→(0,10)→(0,0)
 *   - Odometry accumulates drift
 *   - Loop closure at (0,0) corrects the trajectory via graph optimization
 *
 * Reference: Lu & Milios (1997), Dellaert & Kaess (2006)
 */

#include "slam_types.h"
#include "slam_graph.h"
#include "slam_sensor.h"

extern void slam_config_default(slam_config_t *cfg);
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern double slam_normalize_angle(double theta);

#define NUM_NODES 5

int main(void) {
    printf("=== Pose Graph SLAM: Loop Closure Demo ===\n\n");

    /* 1. Initialize pose graph */
    slam_pose_graph_t graph;
    slam_graph_init(&graph, NUM_NODES, NUM_NODES + 1);

    /* 2. Add nodes (poses along a square path, with drift simulation) */
    /* Ground truth corners of 10×10 square */
    double gt_x[NUM_NODES] = {0, 10, 10, 0, 0};
    double gt_y[NUM_NODES] = {0,  0, 10, 10, 0};
    double gt_th[NUM_NODES] = {0, 0, M_PI/2, M_PI, -M_PI/2};

    /* Drifted poses (simulate odometry error) */
    double drift_x[NUM_NODES] = {0, 10.3, 10.5, 0.3, 0.5};
    double drift_y[NUM_NODES] = {0, -0.2, 9.7, 9.8, -0.3};
    double drift_th[NUM_NODES] = {0, -0.05, M_PI/2+0.03, M_PI-0.04, -M_PI/2-0.06};

    printf("Node  Ground Truth          Drifted (before optimization)\n");
    printf("----  ---------------------  ----------------------------\n");
    for (int i = 0; i < NUM_NODES; i++) {
        slam_pose2d_t p = {drift_x[i], drift_y[i], drift_th[i]};
        int nid;
        slam_graph_add_node(&graph, &p, (i == 0), &nid);
        printf("  %d   (%.1f, %.1f, %5.1f°)     (%.2f, %.2f, %5.1f°)\n",
               i, gt_x[i], gt_y[i], gt_th[i]*180/M_PI,
               drift_x[i], drift_y[i], drift_th[i]*180/M_PI);
    }

    /* 3. Add odometry edges between consecutive nodes */
    /* Information matrix: higher confidence in translation, lower in rotation */
    double odom_info[9] = {100, 0, 0,  0, 100, 0,  0, 0, 50};

    for (int i = 0; i < NUM_NODES - 1; i++) {
        /* Compute relative pose from ground truth (simulated perfect odometry) */
        double dx = gt_x[i+1] - gt_x[i];
        double dy = gt_y[i+1] - gt_y[i];
        double dth = slam_normalize_angle(gt_th[i+1] - gt_th[i]);

        /* Rotate into i's frame */
        double ci = cos(gt_th[i]), si = sin(gt_th[i]);
        double rel_x =  dx * ci + dy * si;
        double rel_y = -dx * si + dy * ci;

        slam_pose2d_t constraint = {rel_x, rel_y, dth};
        int eid;
        slam_graph_add_edge(&graph, i, i+1, &constraint, odom_info, &eid);

        printf("  Edge %d→%d: odometry constraint = (%.2f, %.2f, %.2f°)\n",
               i, i+1, rel_x, rel_y, dth*180/M_PI);
    }

    /* 4. Add loop closure edge (node 4 → node 0) */
    /* High-confidence loop closure */
    double loop_info[9] = {500, 0, 0,  0, 500, 0,  0, 0, 200};
    slam_pose2d_t loop_constraint = {0.0, 0.0, 0.0};  /* expected: back to origin */
    int eid;
    slam_graph_add_edge(&graph, NUM_NODES-1, 0, &loop_constraint, loop_info, &eid);
    graph.edges[eid].is_loop_closure = true;
    printf("  Edge %d→%d: LOOP CLOSURE constraint = (0, 0, 0°)\n",
           NUM_NODES-1, 0);

    /* 5. Run optimization */
    slam_config_t config;
    slam_config_default(&config);
    config.max_graph_iterations = 20;
    config.convergence_thresh = 1e-6;

    printf("\nOptimizing...\n");
    slam_graph_optimize_gauss_newton(&graph, &config, NULL);

    /* 6. Print results */
    printf("\nOptimized trajectory:\n");
    printf("Node  X        Y        Theta      Error(m)\n");
    printf("----  -------  -------  ---------  ---------\n");

    double total_err = 0.0;
    for (int i = 0; i < NUM_NODES; i++) {
        double err = sqrt((graph.nodes[i].pose.x - gt_x[i]) *
                          (graph.nodes[i].pose.x - gt_x[i])
                        + (graph.nodes[i].pose.y - gt_y[i]) *
                          (graph.nodes[i].pose.y - gt_y[i]));
        total_err += err;
        printf("  %d   %7.3f  %7.3f  %7.2f°  %9.3f\n",
               i, graph.nodes[i].pose.x, graph.nodes[i].pose.y,
               graph.nodes[i].pose.theta * 180.0 / M_PI, err);
    }

    printf("\nResults:\n");
    printf("  Average node error: %.3f m\n", total_err / NUM_NODES);
    printf("  Loop closure error: %.3f m\n",
           sqrt((graph.nodes[NUM_NODES-1].pose.x - 0.0) *
                (graph.nodes[NUM_NODES-1].pose.x - 0.0)
              + (graph.nodes[NUM_NODES-1].pose.y - 0.0) *
                (graph.nodes[NUM_NODES-1].pose.y - 0.0)));
    printf("  Chi2 iterations: %d\n", graph.chi2_len);

    slam_graph_free(&graph);
    printf("\nGraph SLAM demo complete.\n");
    return 0;
}
