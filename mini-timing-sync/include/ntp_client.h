/**
 * @file ntp_client.h
 * @brief Network Time Protocol (NTP) client implementation (RFC 5905)
 *
 * L1: NTP packet format, stratum, reference clock identifiers
 * L2: NTP clock filter, clock selection, clock combining
 * L3: NTP clock discipline loop (FLL/PLL hybrid)
 * L4: NTP offset/delay equations, Marzullo's algorithm
 * L5: NTP clock mitigation algorithms, intersection algorithm
 * L6: NTP client-server synchronization
 * L7: Financial timestamping (MiFID II compliance)
 *
 * Reference: Mills, D.L. "Computer Network Time Synchronization" (2011)
 *            IETF RFC 5905 (NTPv4), RFC 5906 (Autokey)
 */

#ifndef NTP_CLIENT_H
#define NTP_CLIENT_H
#include <stdint.h>
#include "timing_sync.h"
#ifdef __cplusplus
extern "C" {
#endif

/* L1: NTP association modes (RFC 5905 Section 3) */
typedef enum {
    NTP_MODE_SYMMETRIC_ACTIVE   = 1,
    NTP_MODE_SYMMETRIC_PASSIVE  = 2,
    NTP_MODE_CLIENT             = 3,
    NTP_MODE_SERVER             = 4,
    NTP_MODE_BROADCAST          = 5,
    NTP_MODE_CONTROL            = 6,
} NtpMode;

/* L1: NTP leap indicator values */
typedef enum {
    NTP_LI_NOWARNING = 0,
    NTP_LI_LAST_MINUTE_61 = 1,
    NTP_LI_LAST_MINUTE_59 = 2,
    NTP_LI_ALARM = 3
} NtpLeapIndicator;

/* L1: NTP timestamp formats */
typedef struct {
    uint32_t seconds;
    uint32_t fraction;
} NtpTimestampShort;

/* L1: NTP reference clock identifiers (common) */
typedef enum {
    NTP_REFID_GOES   = 0x474F4553,  /* GOES */
    NTP_REFID_GPS    = 0x47505300,  /* GPS\0 */
    NTP_REFID_PTP    = 0x50545000,  /* PTP\0 */
    NTP_REFID_LOCL   = 0x4C4F434C,  /* LOCL (uncalibrated local) */
    NTP_REFID_CESM   = 0x4345534D,  /* CESM (cesium clock) */
    NTP_REFID_ATOM   = 0x41544F4D   /* ATOM (atomic clock) */
} NtpRefId;

/* L2: NTP peer (clock filter) sample */
typedef struct {
    Timestamp orig;           /* T1: originate timestamp */
    Timestamp recv;           /* T2: receive timestamp (server) */
    Timestamp xmit;           /* T3: transmit timestamp (server) */
    Timestamp dest;           /* T4: destination timestamp (client) */
    double     offset_ns;     /* Computed offset */
    double     delay_ns;      /* Computed round-trip delay */
    double     dispersion_ns; /* Filter dispersion */
    uint8_t    valid;         /* Sample validity flag */
} NtpSample;

/* L2: NTP clock filter (stores 8 most recent samples per peer) */
#define NTP_FILTER_SIZE 8
typedef struct {
    NtpSample samples[NTP_FILTER_SIZE];
    int       head;       /* Index of newest sample */
    int       count;      /* Number of valid samples */
    double    filter_dispersion_ns;
    int       poll_interval;
    uint8_t   reach;      /* Reachability register (8-bit shift) */
} NtpClockFilter;

/* L2: NTP peer association state */
typedef struct {
    uint32_t        peer_address;
    uint8_t         stratum;
    int8_t          precision;
    double          root_delay_ns;
    double          root_dispersion_ns;
    double          peer_dispersion_ns;
    double          peer_offset_ns;
    double          peer_delay_ns;
    double          root_distance_ns;
    NtpClockFilter  filter;
    uint8_t         reachable;
    uint8_t         selected;   /* Selected for synchronization */
} NtpPeer;

/* L2: NTP client state */
typedef struct {
    NtpPeer         peers[10];
    int             num_peers;
    int             selected_peer_index;
    double          system_offset_ns;
    double          system_delay_ns;
    double          system_dispersion_ns;
    SyncStatus      sync_status;
    PiServoConfig   servo_config;
    PiServoState    servo_state;
    double          poll_interval_s;
    double          freq_offset_ppm;
    double          root_distance_ns;
} NtpClient;

/**
 * L4: Compute NTP offset and delay from four timestamps.
 *
 * offset = ((T2 - T1) + (T3 - T4)) / 2
 * delay  = (T4 - T1) - (T3 - T2)
 */
int ntp_compute_offset_delay(const NtpSample *sample);

/**
 * L4: NTP clock filter algorithm (RFC 5905 Section 10).
 *
 * Maintains 8 most recent samples, selects the one with
 * minimum delay as the most accurate offset estimate.
 *
 * @param filter  Clock filter (updated in place)
 * @param sample  New sample to insert
 * @return        1 if filter produced a new best offset
 */
int ntp_clock_filter(NtpClockFilter *filter, const NtpSample *sample);

/**
 * L5: NTP intersection algorithm (Marzullo's algorithm).
 *
 * Given multiple peers with (offset, root_distance) pairs,
 * finds the largest intersection interval.
 * This is the key algorithm for selecting truechimers from falsetickers.
 *
 * @param offsets_ns        Array of peer offsets
 * @param root_distances_ns Array of peer root distances
 * @param num_peers         Number of peers
 * @param selection_mask    [out] Bitmask of selected peers (1 = selected)
 * @param consensus_offset  [out] Consensus offset estimate
 * @param consensus_error   [out] Error bound on consensus
 * @return Number of selected peers (0 if no intersection)
 */
int ntp_intersection_algorithm(const double *offsets_ns,
                               const double *root_distances_ns,
                               int num_peers,
                               uint32_t *selection_mask,
                               double *consensus_offset,
                               double *consensus_error);

/**
 * L5: NTP clock combining algorithm.
 *
 * Averages the offsets of all selected peers, weighted by
 * the inverse of their root distance.
 *
 * @param offsets_ns       Array of peer offsets
 * @param root_distances_ns Array of peer root distances
 * @param selected_mask    Bitmask of selected peers
 * @param num_peers        Number of peers
 * @return Combined offset estimate in nanoseconds
 */
double ntp_combine_clocks(const double *offsets_ns,
                          const double *root_distances_ns,
                          uint32_t selected_mask, int num_peers);

/**
 * L6: NTP client update: process new peer samples.
 *
 * Full NTP client processing:
 * 1. Insert sample into clock filter
 * 2. If new best offset from filter, update peer offset/delay
 * 3. Run intersection algorithm across all peers
 * 4. Combine selected peers
 * 5. Run clock discipline (PI servo)
 *
 * @param client         NTP client state
 * @param peer_index     Which peer received new sample
 * @param sample         New sample data
 * @return Clock correction in nanoseconds
 */
double ntp_client_update(NtpClient *client, int peer_index,
                         const NtpSample *sample);

/**
 * L6: Initialize NTP client with default parameters.
 */
void ntp_client_init(NtpClient *client);

/**
 * L6: Add or update a peer in the NTP client's peer list.
 *
 * @param client  NTP client state
 * @param address Peer identifier
 * @param stratum Peer's stratum
 * @return Peer index (or -1 if peer list full)
 */
int ntp_add_peer(NtpClient *client, uint32_t address, uint8_t stratum);

/**
 * L7: MiFID II / MiFIR timestamping compliance check.
 *
 * Financial regulations require timestamp accuracy:
 * - MiFID II RTS 25: max divergence from UTC < 100 us
 * - High-frequency trading: accuracy < 1 us
 * - CAT NMS (US): accuracy < 50 ms, clock sync every second
 *
 * @param offset_ns      Current UTC offset
 * @param max_divergence_ns  Maximum allowed divergence
 * @return 1 if compliant, 0 if not
 */
int ntp_mifid_compliance(double offset_ns, double max_divergence_ns);

/**
 * L5: Convert NTP short timestamp to double seconds.
 */
double ntp_short_to_double(NtpTimestampShort nts);

/**
 * L5: Convert double seconds to NTP short timestamp.
 */
NtpTimestampShort ntp_double_to_short(double seconds);

#ifdef __cplusplus
}
#endif
#endif
