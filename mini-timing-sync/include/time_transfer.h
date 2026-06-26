/**
 * @file time_transfer.h
 * @brief Time transfer methods: two-way, common-view, all-in-view
 *
 * L1: One-way, two-way, common-view time transfer definitions
 * L2: Asymmetric delay compensation
 * L3: Kalman-filtered time transfer
 * L4: Two-way time transfer equation with asymmetric delay
 * L5: Common-view GPS time transfer algorithm
 * L6: Time transfer uncertainty budget
 * L7: Two-way satellite time and frequency transfer (TWSTFT)
 * L8: White Rabbit precision time protocol extensions
 *
 * Reference: ITU-R TF.1153 (Two-way time transfer)
 *            CCTF Guidelines for time transfer
 *            White Rabbit Project (CERN, GSI)
 *
 * Course: Stanford EE359, ETH 227-0436, TU Munich
 */

#ifndef TIME_TRANSFER_H
#define TIME_TRANSFER_H
#include <stdint.h>
#include "timing_sync.h"
#ifdef __cplusplus
extern "C" {
#endif

/* L1: Time transfer method types */
typedef enum {
    TT_ONE_WAY          = 0,  /* One-way time transfer */
    TT_TWO_WAY          = 1,  /* Two-way time transfer (TWTT) */
    TT_COMMON_VIEW      = 2,  /* Common-view time transfer */
    TT_ALL_IN_VIEW      = 3,  /* All-in-view time transfer */
    TT_PPP              = 4,  /* Precise point positioning time transfer */
    TT_WHITE_RABBIT     = 5   /* White Rabbit sub-ns time transfer */
} TimeTransferMethod;

/* L1: Time transfer link model */
typedef struct {
    double forward_delay_ns;     /* One-way forward path delay */
    double reverse_delay_ns;     /* One-way reverse path delay */
    double asymmetry_ns;         /* (forward - reverse) / 2 */
    double delay_variance_ns2;   /* Delay measurement variance */
    double temperature_coeff_ps_per_c; /* Temperature coefficient */
} TimeTransferLink;

/* L2: Asymmetric delay model */
typedef struct {
    double static_asymmetry_ns;     /* Fixed asymmetry [ns] */
    double dynamic_asymmetry_ns;    /* Time-varying component [ns] */
    double asymmetry_uncertainty_ns; /* Uncertainty in asymmetry estimate */
    double asymmetry_drift_ns_per_hour; /* Drift rate of asymmetry */
} AsymmetryModel;

/* L2: Time transfer result */
typedef struct {
    double offset_ns;             /* Clock offset between two sites */
    double delay_ns;              /* One-way path delay */
    double uncertainty_ns;        /* Combined standard uncertainty */
    double link_instability_1s;   /* TDEV at 1 s [seconds] */
    Timestamp measurement_time;   /* Epoch of measurement */
    TimeTransferMethod method;    /* Method used */
    int    valid;                 /* 1 if result is valid */
} TimeTransferResult;

/* L1: Common-view GPS observation */
typedef struct {
    Timestamp gps_time;       /* GPS time of observation */
    int       satellite_prn;  /* GPS satellite PRN number */
    double    local_clock_offset_ns; /* Measured offset of local clock */
    double    elevation_deg;  /* Satellite elevation angle */
    double    iono_delay_ns;  /* Ionospheric delay estimate */
    double    tropo_delay_ns; /* Tropospheric delay estimate */
} GpsCommonViewObs;

/**
 * L4: Two-way time transfer with asymmetric delay compensation.
 *
 * offset = ((t2 - t1) - (t4 - t3)) / 2 + asymmetry
 * delay  = ((t2 - t1) + (t4 - t3)) / 2
 *
 * @param ts           PTP timestamps
 * @param asymmetry_ns Forward path - reverse path delay difference / 2 [ns]
 * @param offset       [out] Compensated offset [ns]
 * @param delay        [out] Mean one-way delay [ns]
 * @return 0 on success
 */
int twtt_with_asymmetry(const PtpTimestamps *ts, double asymmetry_ns,
                        double *offset, double *delay);

/**
 * L5: Common-view GPS time transfer.
 *
 * Two stations A and B observe the same GPS satellite at the
 * same time. The clock offset between A and B is:
 *   offset_AB = (obs_A - obs_B) + corrections_A - corrections_B
 *
 * This function processes a pair of common-view observations
 * to compute the time offset between stations.
 *
 * @param obs_A        Observation at station A
 * @param obs_B        Observation at station B
 * @param sta_A_pos_x  Station A ECEF X [m] (for geometric corrections)
 * @param sta_A_pos_y  Station A ECEF Y [m]
 * @param sta_A_pos_z  Station A ECEF Z [m]
 * @param sta_B_pos_x  Station B ECEF X [m]
 * @param sta_B_pos_y  Station B ECEF Y [m]
 * @param sta_B_pos_z  Station B ECEF Z [m]
 * @return Time transfer result
 */
TimeTransferResult common_view_transfer(const GpsCommonViewObs *obs_A,
                                        const GpsCommonViewObs *obs_B,
                                        double sta_A_pos_x, double sta_A_pos_y,
                                        double sta_A_pos_z,
                                        double sta_B_pos_x, double sta_B_pos_y,
                                        double sta_B_pos_z);

/**
 * L5: All-in-view GPS time transfer.
 *
 * Unlike common-view (same satellite), all-in-view uses different
 * satellites observed at the same epoch, linked through precise
 * satellite clock and orbit products.
 *
 * @param obs_A         Observations at station A
 * @param num_obs_A     Number of observations at A
 * @param obs_B         Observations at station B
 * @param num_obs_B     Number of observations at B
 * @param method        0 = simple mean, 1 = weighted by elevation
 * @return Time transfer result (offset between stations)
 */
TimeTransferResult all_in_view_transfer(const GpsCommonViewObs *obs_A,
                                        int num_obs_A,
                                        const GpsCommonViewObs *obs_B,
                                        int num_obs_B, int method);

/**
 * L5: Estimate time transfer uncertainty budget.
 *
 * Components:
 * - Timestamp resolution
 * - Path delay asymmetry uncertainty
 * - Calibration uncertainty
 * - Environmental (temperature) variation
 * - Multipath/reflections
 * - Instrumentation noise
 *
 * @param link       Link model parameters
 * @param asymmetry  Asymmetry model
 * @param temp_c     Ambient temperature [Celsius]
 * @param result     [out] Uncertainty budget (uncertainty_ns filled)
 */
void time_transfer_uncertainty(const TimeTransferLink *link,
                               const AsymmetryModel *asymmetry,
                               double temp_c,
                               TimeTransferResult *result);

/**
 * L6: TWSTFT (Two-Way Satellite Time and Frequency Transfer) delay model.
 *
 * Uses geostationary satellite relay:
 *   up_delay = d_up_1 + d_up_2 (uplinks to satellite)
 *   down_delay = d_down_1 + d_down_2 (downlinks)
 *   total_delay = up_delay + down_delay + sat_delay
 *
 * Sagnac correction needed due to Earth rotation (relativistic).
 *
 * @param sat_longitude_deg  Satellite longitude [degrees]
 * @param sta_A_lat_deg      Station A latitude [degrees]
 * @param sta_A_lon_deg      Station A longitude [degrees]
 * @param sta_A_alt_m        Station A altitude [meters]
 * @param sta_B_lat_deg      Station B latitude [degrees]
 * @param sta_B_lon_deg      Station B longitude [degrees]
 * @param sta_B_alt_m        Station B altitude [meters]
 * @param sagnac_correction  [out] Sagnac correction [ns]
 * @param total_delay        [out] Estimated total path delay [ns]
 * @return 0 on success
 */
int twstft_delay_model(double sat_longitude_deg,
                       double sta_A_lat_deg, double sta_A_lon_deg,
                       double sta_A_alt_m,
                       double sta_B_lat_deg, double sta_B_lon_deg,
                       double sta_B_alt_m,
                       double *sagnac_correction, double *total_delay);

/**
 * L6: Sagnac effect correction for time transfer.
 *
 * Due to Earth rotation, time transfer signals traveling eastward
 * experience different path lengths than westbound.
 *
 * Correction = (2 * omega_E * A) / c^2
 *
 * Where omega_E = 7.2921159e-5 rad/s (Earth rotation rate)
 * A = projected area on equatorial plane
 * c = 299792458 m/s
 *
 * @param sta_A_lat_deg  Station A latitude
 * @param sta_A_lon_deg  Station A longitude
 * @param sta_B_lat_deg  Station B latitude
 * @param sta_B_lon_deg  Station B longitude
 * @return Sagnac correction in nanoseconds
 */
double sagnac_correction_ns(double sta_A_lat_deg, double sta_A_lon_deg,
                            double sta_B_lat_deg, double sta_B_lon_deg);

/**
 * L8: White Rabbit link model calibration.
 *
 * White Rabbit achieves sub-nanosecond synchronization by:
 * 1. Precise hardware timestamps
 * 2. Link delay model (fiber asymmetry characterization)
 * 3. Phase tracking of recovered clock
 *
 * @param fiber_length_m        One-way fiber length [meters]
 * @param fiber_refractive_index Refractive index of fiber core
 * @param tx_delay_ns           Transmitter electronics delay
 * @param rx_delay_ns           Receiver electronics delay
 * @param forward_delay         [out] Forward path delay [ns]
 * @param reverse_delay         [out] Reverse path delay [ns]
 */
void white_rabbit_link_model(double fiber_length_m,
                             double fiber_refractive_index,
                             double tx_delay_ns, double rx_delay_ns,
                             double *forward_delay,
                             double *reverse_delay);

/**
 * L8: White Rabbit phase tracking for sub-ns alignment.
 *
 * Uses DDMTD (Digital Dual Mixer Time Difference) phase measurement.
 * Phase alignment = main_PLL_phase - helper_PLL_phase
 *
 * @param main_phase_ps    Main PLL phase [picoseconds]
 * @param helper_phase_ps  Helper PLL phase [picoseconds]
 * @param ref_freq_hz      Reference frequency [Hz]
 * @return Phase offset in picoseconds
 */
double white_rabbit_phase_track(double main_phase_ps,
                                double helper_phase_ps,
                                double ref_freq_hz);

#ifdef __cplusplus
}
#endif
#endif
