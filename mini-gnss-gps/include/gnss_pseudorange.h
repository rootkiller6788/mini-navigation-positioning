#ifndef GNSS_PSEUDORANGE_H
#define GNSS_PSEUDORANGE_H
#include "gnss_common.h"
#include "gnss_signal.h"

#ifdef __cplusplus
extern "C" {
#endif

/* =========================================================================
 * L1: Pseudorange Measurement
 *
 * P = ρ + c·(dt_r - dt_s) + I + T + M + ε
 *
 * P   : measured pseudorange [m]
 * ρ   : true geometric range [m]
 * dt_r: receiver clock bias [s]  (unknown, estimated)
 * dt_s: satellite clock bias [s] (from navigation message)
 * I   : ionospheric delay [m]   (frequency-dependent)
 * T   : tropospheric delay [m]  (frequency-independent)
 * M   : multipath error [m]
 * ε   : receiver thermal noise [m]
 * ========================================================================= */

/** Raw pseudorange measurement for one satellite */
typedef struct {
    gnss_satid_t    sat_id;
    gnss_gpstime_t  time_rx;      /* Receiver time of reception */
    double          pseudorange;   /* Raw pseudorange [m] */
    double          snr_db;        /* Carrier-to-noise ratio C/N₀ [dB-Hz] */
    double          sat_elevation; /* Elevation angle [rad] */
    double          sat_azimuth;   /* Azimuth angle [rad] */
    double          doppler_hz;    /* Measured Doppler shift [Hz] */
} gnss_pseudorange_t;

/* -------------------------------------------------------------------------
 * L1: Pseudorange Corrections Bundle
 * ------------------------------------------------------------------------- */

typedef struct {
    double sat_clock_corr;     /* Satellite clock correction [m] (c·dt_s) */
    double iono_corr;          /* Ionospheric correction [m] */
    double tropo_corr;         /* Tropospheric correction [m] */
    double sagnac_corr;        /* Sagnac (Earth rotation) correction [m] */
    double relativistic_corr;  /* Relativistic clock correction [m] */
    double total_corr;         /* Sum of all corrections [m] */
} gnss_range_corrections_t;

/* -------------------------------------------------------------------------
 * L1: Pseudorange Error Budget (UERE analysis)
 * ------------------------------------------------------------------------- */

typedef struct {
    double uere_sat_clock;     /* Satellite clock & ephemeris error [m] */
    double uere_iono_residual; /* Residual ionospheric error [m] */
    double uere_tropo_residual;/* Residual tropospheric error [m] */
    double uere_multipath;     /* Multipath error [m] */
    double uere_receiver_noise;/* Receiver noise [m] */
    double uere_total;         /* RSS total UERE [m] */
} gnss_error_budget_t;

/* -------------------------------------------------------------------------
 * L2: Ionospheric Model — Klobuchar (1987)
 * ------------------------------------------------------------------------- */

/**
 * @brief Klobuchar ionospheric model parameters
 *
 * Broadcast in GPS navigation message subframe 4.
 *   α₀, α₁, α₂, α₃: ionospheric amplitude coefficients [s]
 *   β₀, β₁, β₂, β₃: ionospheric period coefficients [s]
 *
 * The model approximates zenith ionospheric delay as a half-cosine
 * during daytime and constant during nighttime:
 *   I_z = A₁ + A₂·cos(2π·(t - A₃) / A₄)
 *   I   = c · F · I_z   (mapped to line-of-sight)
 *
 * where F = 1 + 16·(0.53 - E)³  is the obliquity factor.
 *
 * Reference: Klobuchar, J.A. (1987). "Ionospheric time-delay algorithm
 * for single-frequency GPS users." IEEE Trans. Aerospace.
 */
typedef struct {
    double alpha[4];   /* α₀,α₁,α₂,α₃ [s],[s/semi-circle],[s/semi-circle²],[s/semi-circle³] */
    double beta[4];    /* β₀,β₁,β₂,β₃ [s],[s/semi-circle],[s/semi-circle²],[s/semi-circle³] */
} gnss_klobuchar_params_t;

/* -------------------------------------------------------------------------
 * L2: Tropospheric Models
 * ------------------------------------------------------------------------- */

/**
 * @brief Tropospheric parameters for Saastamoinen model
 *
 * Requires: total pressure P [hPa], temperature T [K],
 *           water vapor partial pressure e [hPa].
 */
typedef struct {
    double pressure_hPa;       /* Total atmospheric pressure [hPa] */
    double temperature_K;      /* Surface temperature [K] */
    double humidity_hPa;       /* Water vapor partial pressure [hPa] */
} gnss_tropo_params_t;

/** Hopfield tropospheric model parameters (dry + wet components) */
typedef struct {
    double T0;                 /* Surface temperature [K] */
    double P0;                 /* Surface pressure [hPa] */
    double RH;                 /* Relative humidity [%] */
    double h0;                 /* Station height [m] */
} gnss_hopfield_params_t;

/* -------------------------------------------------------------------------
 * L2: Multipath Model
 * ------------------------------------------------------------------------- */

/** Single-reflector multipath */
typedef struct {
    double relative_amplitude; /* α = A_reflected / A_direct  (0 < α < 1) */
    double path_delay_m;       /* Extra path length [m] */
    double phase_shift_rad;    /* Relative carrier phase [rad] */
} gnss_multipath_t;

/* -------------------------------------------------------------------------
 * L5: Code-Minus-Carrier (CMC) for multipath analysis
 * ------------------------------------------------------------------------- */

/**
 * @brief CMC: Pcode - Lcarrier = 2·I + M_code + ε
 *
 * Since the ionosphere advances code and delays carrier by equal magnitude,
 * the CMC reveals code multipath (M_code) plus twice the ionospheric delay.
 *
 * Useful for detecting and characterizing multipath environment.
 */
typedef struct {
    size_t  n_epochs;
    double *code_minus_carrier; /* CMC [m] per epoch */
    double  mean;               /* Mean CMC bias [m] */
    double  std;                /* CMC standard deviation [m] */
    double  max_abs;            /* Maximum absolute CMC [m] */
} gnss_cmc_series_t;

/* -------------------------------------------------------------------------
 * API: Pseudorange correction models
 * ------------------------------------------------------------------------- */

/**
 * @brief Klobuchar ionospheric delay computation
 *
 * @param params  Alpha/Beta coefficients from nav message
 * @param t       GPS time
 * @param lla     Approximate user position (for geomagnetic mapping)
 * @param sat_pos Satellite ECEF position
 * @param freq_hz Carrier frequency [Hz] (L1=1575.42e6)
 * @return Ionospheric delay [m] (positive = group delay)
 *
 * Complexity: O(1), 7 transcendental functions.
 */
double gnss_iono_klobuchar(const gnss_klobuchar_params_t *params,
                            gnss_gpstime_t t,
                            gnss_lla_t lla,
                            gnss_ecef_t sat_pos,
                            double freq_hz);

/**
 * @brief Saastamoinen tropospheric model
 *
 * T = 0.002277 / cos(z') · [P + (1255/T + 0.05)·e - B·tan²(z')] + δR
 *
 * where z' is the apparent zenith angle corrected for curvature.
 *
 * @param params Atmospheric parameters (P, T, e)
 * @param elevation_rad Satellite elevation angle [rad]
 * @param alt_m User altitude above ellipsoid [m]
 * @return Tropospheric delay [m]
 *
 * Reference: Saastamoinen, J. (1972). "Contributions to the theory of
 * atmospheric refraction." Bulletin Géodésique.
 */
double gnss_tropo_saastamoinen(const gnss_tropo_params_t *params,
                                double elevation_rad, double alt_m);

/**
 * @brief Hopfield tropospheric model (dry + wet components)
 *
 * Separately models dry (hydrostatic) and wet delay.
 * Dry delay decays as (h_d - h)⁴; wet decays as (h_w - h)⁴.
 *
 * Reference: Hopfield, H.S. (1969). "Two-quartic tropospheric refractivity
 * profile." JGR.
 */
double gnss_tropo_hopfield_dry(const gnss_hopfield_params_t *params,
                                double elevation_rad);
double gnss_tropo_hopfield_wet(const gnss_hopfield_params_t *params,
                                double elevation_rad);
double gnss_tropo_hopfield_total(const gnss_hopfield_params_t *params,
                                  double elevation_rad);

/**
 * @brief Simple tropospheric model (UNB3m, used when no meteo data)
 *
 * T = 2.3 / (sin(el) + 0.15) [m]   for temperate regions.
 */
double gnss_tropo_simple(double elevation_rad);

/** @brief Compute all pseudorange corrections */
gnss_range_corrections_t gnss_correct_pseudorange(
    const gnss_pseudorange_t *raw,
    const gnss_ephemeris_t *eph,
    const gnss_klobuchar_params_t *iono,
    const gnss_tropo_params_t *tropo,
    gnss_ecef_t rx_pos);

/** @brief Compute UERE error budget */
gnss_error_budget_t gnss_uere_budget(double elevation_rad, int single_freq);

/** @brief Compute multipath-induced pseudorange error */
double gnss_multipath_error(const gnss_multipath_t *mp,
                             double chip_rate, double correlator_spacing);

/** @brief Compute code-minus-carrier series */
int gnss_cmc_compute(const double *code_range, const double *carrier_phase,
                      size_t n_epochs, double wavelength, gnss_cmc_series_t *cmc);
void gnss_cmc_free(gnss_cmc_series_t *cmc);

#ifdef __cplusplus
}
#endif
#endif /* GNSS_PSEUDORANGE_H */
