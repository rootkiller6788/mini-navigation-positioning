/**
 * @file example_ntp_time_transfer.c
 * @brief NTP client and time transfer example
 *
 * Demonstrates: NTP peer management, clock filter,
 * Marzullo intersection algorithm, clock combining,
 * Sagnac correction, TWSTFT delay model, White Rabbit link.
 *
 * L6 - Canonical Problem: NTP time synchronization
 * L7 - Application: Financial timestamping (MiFID II)
 * L8 - Advanced: White Rabbit sub-ns link model
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "timing_sync.h"
#include "ntp_client.h"
#include "time_transfer.h"

int main(void)
{
    printf("=== NTP Client and Time Transfer Example ===\n\n");

    /* ---- NTP Client Demo ---- */
    printf("--- NTP Client Synchronization ---\n");

    NtpClient client;
    ntp_client_init(&client);

    /* Add three NTP servers (simulating real servers) */
    ntp_add_peer(&client, 0xC0A80101, 1); /* Stratum 1 (GPS) */
    ntp_add_peer(&client, 0xC0A80102, 2); /* Stratum 2 */
    ntp_add_peer(&client, 0xC0A80103, 2); /* Stratum 2 */

    printf("Added 3 NTP peers.\n");

    /* Simulate receiving samples from peers */
    printf("\nSimulating NTP exchanges...\n");
    printf("  Exchange | Peer | Offset(ns) | Delay(ns)\n");
    printf("  ---------|------|-------------|----------\n");

    for (int ex = 0; ex < 5; ex++) {
        for (int peer = 0; peer < 3; peer++) {
            NtpSample sample;
            memset(&sample, 0, sizeof(sample));

            /* Simulate timestamps with different offsets per peer */
            double peer_true_offset = 100.0 + peer * 200.0
                                    + (double)(rand() % 21 - 10);
            double peer_true_delay = 10000.0 + peer * 5000.0;

            sample.orig.seconds = 2000 + ex;
            sample.orig.nanoseconds = 0;

            sample.recv = sample.orig;
            timing_timestamp_add_ns(&sample.recv,
                (int64_t)(peer_true_delay + peer_true_offset));

            sample.xmit = sample.recv;
            timing_timestamp_add_ns(&sample.xmit, 1000);

            sample.dest = sample.xmit;
            timing_timestamp_add_ns(&sample.dest,
                (int64_t)(peer_true_delay - peer_true_offset));

            double correction = ntp_client_update(&client, peer, &sample);

            printf("  %8d | %4d | %10.1f | %9.1f\n",
                   ex + 1, peer + 1,
                   client.peers[peer].peer_offset_ns,
                   client.peers[peer].peer_delay_ns);
        }
    }

    printf("\nFinal NTP status: ");
    switch (client.sync_status) {
    case SYNC_LOCKED:  printf("LOCKED\n"); break;
    default:          printf("UNLOCKED\n"); break;
    }
    printf("System offset: %.1f ns\n", client.system_offset_ns);

    /* MiFID II compliance */
    printf("\nMiFID II Financial Timestamping Compliance:\n");
    int mifid_std = ntp_mifid_compliance(client.system_offset_ns, 100000.0);
    int mifid_hft = ntp_mifid_compliance(client.system_offset_ns, 1000.0);
    printf("  Standard trading (<100us): %s\n",
           mifid_std ? "COMPLIANT" : "NON-COMPLIANT");
    printf("  HFT trading (<1us):       %s\n",
           mifid_hft ? "COMPLIANT" : "NON-COMPLIANT");

    /* ---- Time Transfer Demo ---- */
    printf("\n--- Time Transfer Methods ---\n");

    /* Sagnac effect between global sites */
    printf("\nSagnac Correction (Earth rotation effect):\n");
    double sites[][2] = {
        {35.7, 139.8},    /* Tokyo */
        {37.8, -122.4},   /* San Francisco */
        {48.9, 2.3},      /* Paris */
        {52.5, 13.4}      /* Berlin */
    };
    const char *site_names[] = {"Tokyo", "San Francisco", "Paris", "Berlin"};

    for (int i = 0; i < 4; i++) {
        for (int j = i + 1; j < 4; j++) {
            double sagnac = sagnac_correction_ns(
                sites[i][0], sites[i][1],
                sites[j][0], sites[j][1]);
            printf("  %s <-> %s: %+.1f ns\n",
                   site_names[i], site_names[j], sagnac);
        }
    }

    /* TWSTFT delay model */
    printf("\nTWSTFT Delay Model (GEO satellite at 100E):\n");
    double sagnac_corr, total_delay;
    twstft_delay_model(100.0,  /* Satellite longitude */
                       35.7, 139.8, 50.0,  /* Tokyo */
                       48.9, 2.3, 100.0,   /* Paris */
                       &sagnac_corr, &total_delay);
    printf("  Total path delay: %.1f ms\n", total_delay / 1.0e6);
    printf("  Sagnac correction: %.1f ns\n", sagnac_corr);

    /* White Rabbit link model */
    printf("\nWhite Rabbit Link Model (10 km fiber):\n");
    double fwd_delay, rev_delay;
    white_rabbit_link_model(10000.0, 1.467,  /* 10 km SMF-28 */
                            5.0, 5.0,         /* TX/RX delays */
                            &fwd_delay, &rev_delay);
    printf("  One-way fiber delay: %.1f ns\n", fwd_delay);
    printf("  Asymmetry (fwd-rev):  %.1f ns\n", fwd_delay - rev_delay);

    /* Phase tracking example */
    double phase_offset = white_rabbit_phase_track(125.0, 130.0, 125.0e6);
    printf("  Phase tracking offset: %.3f ps\n", phase_offset);

    printf("\n=== Example Complete ===\n");
    return 0;
}
