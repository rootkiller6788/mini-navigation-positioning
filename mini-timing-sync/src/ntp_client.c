/**
 * @file ntp_client.c
 * @brief NTP client implementation (RFC 5905)
 *
 * Implements: NTP offset/delay computation, clock filter,
 * Marzullo's intersection algorithm, clock combining,
 * NTP client update loop, MiFID II compliance.
 */

#include "ntp_client.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

/* ===================================================================
 * L4: NTP Offset and Delay Computation
 * =================================================================== */

int ntp_compute_offset_delay(const NtpSample *sample)
{
    if (!sample) return -1;

    /* RFC 5905 Section 8:
     *   offset = ((T2 - T1) + (T3 - T4)) / 2
     *   delay  = (T4 - T1) - (T3 - T2)
     *
     * We cast away const to write results back into the sample struct.
     */
    NtpSample *s = (NtpSample *)sample;

    return timing_ntp_offset_delay(&s->orig, &s->recv,
                                   &s->xmit, &s->dest,
                                   &s->offset_ns, &s->delay_ns);
}

/* ===================================================================
 * L4: NTP Clock Filter Algorithm (RFC 5905 Section 10)
 * =================================================================== */

int ntp_clock_filter(NtpClockFilter *filter, const NtpSample *new_sample)
{
    if (!filter || !new_sample) return 0;

    /* Insert new sample into circular buffer */
    filter->head = (filter->head + 1) % NTP_FILTER_SIZE;
    filter->samples[filter->head] = *new_sample;
    filter->samples[filter->head].valid = 1;

    if (filter->count < NTP_FILTER_SIZE) {
        filter->count++;
    }

    /* Find sample with minimum delay (most accurate offset estimate
     * according to NTP's minimum-delay principle) */
    int best_idx = -1;
    double min_delay = 1.0e18;

    for (int i = 0; i < filter->count; i++) {
        if (filter->samples[i].valid &&
            filter->samples[i].delay_ns < min_delay) {
            min_delay = filter->samples[i].delay_ns;
            best_idx = i;
        }
    }

    /* Compute filter dispersion as weighted sum of offset differences
     * from the minimum-delay sample */
    if (best_idx >= 0 && filter->count > 1) {
        double sum_disp = 0.0;
        double best_offset = filter->samples[best_idx].offset_ns;

        for (int i = 0; i < filter->count; i++) {
            if (filter->samples[i].valid) {
                double diff = filter->samples[i].offset_ns - best_offset;
                sum_disp += fabs(diff);
            }
        }
        filter->filter_dispersion_ns = sum_disp / (double)filter->count;
    }

    return (best_idx >= 0) ? 1 : 0;
}

/* ===================================================================
 * L5: Marzullo's Intersection Algorithm (NTP Truechimer Selection)
 * =================================================================== */

/* Helper structure for interval endpoints used in Marzullo's algorithm */
typedef struct {
    double value;    /* Endpoint value */
    int    type;     /* +1 for low endpoint, -1 for high endpoint */
} MarzulloEndpoint;

static int marzullo_compare_endpoints(const void *a, const void *b)
{
    const MarzulloEndpoint *ea = (const MarzulloEndpoint *)a;
    const MarzulloEndpoint *eb = (const MarzulloEndpoint *)b;

    if (ea->value < eb->value) return -1;
    if (ea->value > eb->value) return 1;
    /* Put high endpoints (-1) before low endpoints (+1) at same value */
    return ea->type - eb->type;
}

int ntp_intersection_algorithm(const double *offsets_ns,
                               const double *root_distances_ns,
                               int num_peers,
                               uint32_t *selection_mask,
                               double *consensus_offset,
                               double *consensus_error)
{
    if (!offsets_ns || !root_distances_ns || !selection_mask ||
        !consensus_offset || !consensus_error) {
        return 0;
    }

    *selection_mask = 0;
    *consensus_offset = 0.0;
    *consensus_error = 0.0;

    if (num_peers < 1) return 0;

    /* For each peer, define confidence interval:
     *   [offset - root_distance, offset + root_distance]
     *
     * Marzullo's algorithm: find the largest interval intersected
     * by the most peers.
     */

    int num_endpoints = 2 * num_peers;
    MarzulloEndpoint *endpoints =
        (MarzulloEndpoint *)malloc((size_t)num_endpoints *
                                    sizeof(MarzulloEndpoint));
    if (!endpoints) return 0;

    for (int i = 0; i < num_peers; i++) {
        endpoints[2*i].value = offsets_ns[i] - root_distances_ns[i];
        endpoints[2*i].type = 1; /* Low endpoint */

        endpoints[2*i + 1].value = offsets_ns[i] + root_distances_ns[i];
        endpoints[2*i + 1].type = -1; /* High endpoint */
    }

    /* Sort endpoints by value */
    qsort(endpoints, (size_t)num_endpoints, sizeof(MarzulloEndpoint),
          marzullo_compare_endpoints);

    /* Scan to find maximum intersection */
    int current_count = 0;
    int max_count = 0;
    double best_low = 0.0, best_high = 0.0;
    double low_mark = 0.0;
    int in_max_region = 0;

    for (int i = 0; i < num_endpoints; i++) {
        current_count += endpoints[i].type; /* +1 for low, -1 for high */

        if (current_count > max_count) {
            /* Entering a new maximum intersection region */
            max_count = current_count;
            low_mark = endpoints[i].value;
            in_max_region = 1;
        } else if (in_max_region && current_count < max_count) {
            /* Leaving the maximum intersection region */
            best_low = low_mark;
            best_high = endpoints[i].value;
            in_max_region = 0;
        }
    }

    /* If only one peer, use its interval directly */
    if (num_peers == 1) {
        max_count = 1;
        best_low = offsets_ns[0] - root_distances_ns[0];
        best_high = offsets_ns[0] + root_distances_ns[0];
    }

    /* Build selection mask: which peers intersect the best interval? */
    for (int i = 0; i < num_peers && i < 32; i++) {
        double low_i = offsets_ns[i] - root_distances_ns[i];
        double high_i = offsets_ns[i] + root_distances_ns[i];

        /* Peer interval intersects consensus interval? */
        if (low_i <= best_high && high_i >= best_low) {
            *selection_mask |= (1u << i);
        }
    }

    *consensus_offset = (best_low + best_high) / 2.0;
    *consensus_error = (best_high - best_low) / 2.0;

    free(endpoints);
    return max_count;
}

/* ===================================================================
 * L5: NTP Clock Combining Algorithm
 * =================================================================== */

double ntp_combine_clocks(const double *offsets_ns,
                          const double *root_distances_ns,
                          uint32_t selected_mask, int num_peers)
{
    if (!offsets_ns || !root_distances_ns) return 0.0;

    /* Weighted average by inverse of root distance */
    double sum_weight = 0.0;
    double weighted_sum = 0.0;

    for (int i = 0; i < num_peers && i < 32; i++) {
        if (selected_mask & (1u << i)) {
            double weight = 1.0 / (root_distances_ns[i] + 1.0); /* +1 to avoid div/0 */
            weighted_sum += offsets_ns[i] * weight;
            sum_weight += weight;
        }
    }

    if (sum_weight > 0.0) {
        return weighted_sum / sum_weight;
    }
    return 0.0;
}

/* ===================================================================
 * L6: NTP Client Initialization and Update
 * =================================================================== */

void ntp_client_init(NtpClient *client)
{
    if (!client) return;

    memset(client, 0, sizeof(NtpClient));
    client->selected_peer_index = -1;
    client->sync_status = SYNC_FREE_RUNNING;
    client->poll_interval_s = 64.0; /* Default: poll every 64 seconds */

    pi_servo_init(&client->servo_config, &client->servo_state);
}

int ntp_add_peer(NtpClient *client, uint32_t address, uint8_t stratum)
{
    if (!client) return -1;
    if (client->num_peers >= 10) return -1;

    int idx = client->num_peers;
    memset(&client->peers[idx], 0, sizeof(NtpPeer));
    client->peers[idx].peer_address = address;
    client->peers[idx].stratum = stratum;
    client->peers[idx].reachable = 1;
    client->peers[idx].selected = 0;

    client->num_peers++;
    return idx;
}

double ntp_client_update(NtpClient *client, int peer_index,
                         const NtpSample *sample)
{
    if (!client || !sample) return 0.0;
    if (peer_index < 0 || peer_index >= client->num_peers) return 0.0;

    NtpPeer *peer = &client->peers[peer_index];

    /* 1. Compute offset and delay */
    NtpSample working_sample = *sample;
    if (ntp_compute_offset_delay(&working_sample) != 0) {
        return 0.0;
    }

    /* 2. Insert into clock filter */
    ntp_clock_filter(&peer->filter, &working_sample);

    /* 3. Update peer statistics from filter's best sample */
    int best_idx = -1;
    double min_delay = 1.0e18;
    for (int i = 0; i < peer->filter.count; i++) {
        if (peer->filter.samples[i].valid &&
            peer->filter.samples[i].delay_ns < min_delay) {
            min_delay = peer->filter.samples[i].delay_ns;
            best_idx = i;
        }
    }
    if (best_idx >= 0) {
        peer->peer_offset_ns = peer->filter.samples[best_idx].offset_ns;
        peer->peer_delay_ns = peer->filter.samples[best_idx].delay_ns;
        peer->peer_dispersion_ns = peer->filter.filter_dispersion_ns;
        /* Root distance = delay/2 + root_dispersion + peer_dispersion */
        peer->root_distance_ns = peer->peer_delay_ns / 2.0
                               + peer->root_dispersion_ns
                               + peer->peer_dispersion_ns;
    }

    /* 4. Collect offsets and root distances from all reachable peers */
    double offsets[10], root_dists[10];
    int valid_count = 0;
    for (int i = 0; i < client->num_peers; i++) {
        if (client->peers[i].reachable && client->peers[i].filter.count > 0) {
            offsets[valid_count] = client->peers[i].peer_offset_ns;
            root_dists[valid_count] = client->peers[i].root_distance_ns;
            valid_count++;
        }
    }

    if (valid_count == 0) return 0.0;

    /* 5. Marzullo intersection algorithm */
    uint32_t selection_mask = 0;
    double consensus_offset, consensus_error;
    int n_selected = ntp_intersection_algorithm(offsets, root_dists,
                                                valid_count,
                                                &selection_mask,
                                                &consensus_offset,
                                                &consensus_error);

    /* 6. Combine selected clocks */
    double combined_offset = ntp_combine_clocks(offsets, root_dists,
                                                selection_mask, valid_count);

    /* 7. Run clock discipline (PI servo) */
    client->system_offset_ns = combined_offset;
    double correction = pi_servo_update(&client->servo_config,
                                        &client->servo_state,
                                        combined_offset);

    /* 8. Update sync status */
    if (n_selected >= 1 && fabs(combined_offset) < 1000.0) {
        client->sync_status = SYNC_LOCKED;
    } else if (n_selected >= 1) {
        client->sync_status = SYNC_ACQUIRING;
    } else {
        client->sync_status = SYNC_FREE_RUNNING;
    }

    return correction;
}

/* ===================================================================
 * L7: MiFID II Financial Timestamping Compliance
 * =================================================================== */

int ntp_mifid_compliance(double offset_ns, double max_divergence_ns)
{
    /* MiFID II RTS 25:
     * - Timestamp accuracy: max divergence from UTC
     * - High-Frequency Trading: < 1 us (1000 ns)
     * - Standard trading: < 100 us (100000 ns)
     *
     * NIST / US CAT NMS:
     * - Clock sync required every 1 second
     * - Maximum divergence: 50 ms
     */

    if (fabs(offset_ns) <= max_divergence_ns) {
        return 1; /* Compliant */
    }
    return 0; /* Non-compliant */
}

/* ===================================================================
 * L5: NTP Timestamp Format Conversion
 * =================================================================== */

double ntp_short_to_double(NtpTimestampShort nts)
{
    /* NTP short format: 32-bit seconds + 32-bit fraction
     * Fraction range: 0 to ~1 (2^32 units per second)
     */
    return (double)nts.seconds + (double)nts.fraction / 4294967296.0;
}

NtpTimestampShort ntp_double_to_short(double seconds)
{
    NtpTimestampShort nts;
    double int_part, frac_part;

    frac_part = modf(seconds, &int_part);

    nts.seconds = (uint32_t)(uint64_t)int_part;
    nts.fraction = (uint32_t)(frac_part * 4294967296.0);

    return nts;
}
