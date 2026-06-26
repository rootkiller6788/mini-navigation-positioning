#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "../include/uwb_types.h"
#include "../include/uwb_ranging.h"

static int test_anchor_init(void) {
    uwb_anchor_t anchor;
    uwb_anchor_init(&anchor, 1, 10.0, 20.0, 0.0);
    assert(anchor.id == 1);
    assert(fabs(anchor.position.x - 10.0) < 1e-10);
    assert(fabs(anchor.position.y - 20.0) < 1e-10);
    assert(anchor.is_active == 1);
    printf("PASS: test_anchor_init\n");
    return 1;
}

static int test_distance_2d(void) {
    uwb_pos2d_t a = {0, 0}, b = {3, 4};
    double d = uwb_distance_2d(&a, &b);
    assert(fabs(d - 5.0) < 1e-10);
    printf("PASS: test_distance_2d\n");
    return 1;
}

static int test_distance_3d(void) {
    uwb_pos3d_t a = {0, 0, 0}, b = {1, 2, 2};
    double d = uwb_distance_3d(&a, &b);
    assert(fabs(d - 3.0) < 1e-10);
    printf("PASS: test_distance_3d\n");
    return 1;
}

static int test_tof_conversion(void) {
    double d = 10.0;
    double tof = uwb_distance_to_tof(d);
    double d2 = uwb_tof_to_distance(tof);
    assert(fabs(d - d2) < 1e-6);
    printf("PASS: test_tof_conversion\n");
    return 1;
}

static int test_crlb(void) {
    double crlb = uwb_crlb_distance(100.0, 500e6);
    assert(crlb > 0.0 && crlb < 0.1);
    printf("PASS: test_crlb (CRLB=%.6f m)\n", crlb);
    return 1;
}

static int test_ss_twr_tof(void) {
    twr_ss_timestamps_t ts;
    double tick_ps = 15.65; /* DW1000 tick period */
    ts.t_poll_tx = 1000;
    ts.t_poll_rx = 2000; /* 1000 ticks propagation + processing */
    ts.t_resp_tx = 3000;
    ts.t_resp_rx = 4000;
    double tof = twr_ss_compute_tof(&ts, tick_ps);
    assert(tof > 0.0);
    printf("PASS: test_ss_twr_tof (ToF=%.3f ns)\n", tof * 1e9);
    return 1;
}

static int test_ds_twr_tof(void) {
    twr_ds_timestamps_t ts;
    double tick_ps = 15.65;
    ts.t_poll_tx = 1000; ts.t_poll_rx = 1500;
    ts.t_resp_tx = 2000; ts.t_resp_rx = 2500;
    ts.t_final_tx = 3000; ts.t_final_rx = 3500;
    double tof = twr_ds_compute_tof(&ts, tick_ps);
    assert(tof > 0.0);
    printf("PASS: test_ds_twr_tof (ToF=%.3f ns)\n", tof * 1e9);
    return 1;
}

static int test_ranging_stats(void) {
    double dists[] = {5.0, 5.1, 5.2, 4.9, 5.0, 15.0, 5.1, 5.0};
    ranging_burst_stats_t stats;
    ranging_burst_compute_stats(dists, 8, &stats);
    assert(stats.count == 8);
    assert(fabs(stats.median - 5.05) < 0.1);
    printf("PASS: test_ranging_stats (median=%.3f, trimmed_mean=%.3f)\n",
           stats.median, stats.trimmed_mean);
    return 1;
}

static int test_nlos_score(void) {
    double los[] = {5.0, 5.05, 4.98, 5.02, 4.99, 5.01, 5.03, 4.97};
    double nlos[] = {5.0, 5.5, 6.0, 5.2, 7.0, 5.8, 6.5, 5.3};
    double score_los = ranging_burst_nlos_score(los, 8);
    double score_nlos = ranging_burst_nlos_score(nlos, 8);
    assert(score_nlos > score_los);
    printf("PASS: test_nlos_score (LOS=%.3f, NLOS=%.3f)\n", score_los, score_nlos);
    return 1;
}

static int test_error_budget(void) {
    double total = twr_error_budget(100.0, 500e6, 100.0, 0.2, 0.1, 0.02);
    assert(total > 0.0 && total < 1.0);
    printf("PASS: test_error_budget (total=%.4f m)\n", total);
    return 1;
}

static int test_ewma(void) {
    double hist[] = {5.0, 5.1, 5.05, 5.08, 5.03};
    double filt = ranging_ewma_filter(hist, 5, 0.3);
    assert(filt > 5.0 && filt < 5.2);
    printf("PASS: test_ewma (filtered=%.4f)\n", filt);
    return 1;
}

int main(void) {
    int passed = 0;
    passed += test_anchor_init();
    passed += test_distance_2d();
    passed += test_distance_3d();
    passed += test_tof_conversion();
    passed += test_crlb();
    passed += test_ss_twr_tof();
    passed += test_ds_twr_tof();
    passed += test_ranging_stats();
    passed += test_nlos_score();
    passed += test_error_budget();
    passed += test_ewma();
    printf("\n=== %d/11 tests passed ===\n", passed);
    return (passed == 11) ? 0 : 1;
}
