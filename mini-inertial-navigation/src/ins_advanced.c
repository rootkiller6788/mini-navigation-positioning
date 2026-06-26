/**
 * @file    ins_advanced.c
 * @brief   Advanced INS topics: particle filter, transfer alignment, coning/sculling analysis
 *
 * Knowledge Coverage:
 *   L8 (Advanced Topics): Particle filter for non-Gaussian INS errors,
 *       transfer alignment for in-flight calibration, Bortz coning analysis,
 *       higher-order sculling algorithms
 *
 * Reference:
 *   Arulampalam et al. (2002), "A Tutorial on Particle Filters",
 *     IEEE Trans. Signal Process., 50(2): 174-188.
 *   Groves (2013), Chapter 15, "Transfer Alignment".
 *   Savage (1998), "Strapdown Inertial Navigation Integration Algorithm
 *     Design", J. Guid. Control Dyn., 21(1): 19-28 and 21(2): 208-221.
 *
 * Course Mapping:
 *   MIT 6.437 - Inference and Information (particle filtering)
 *   Stanford AA272 - GPS (advanced integration)
 *   Caltech CDS 110 - Optimal Control (estimation)
 */

#include "ins_core.h"
#include "ins_attitude.h"
#include "ins_mechanization.h"
#include "ins_integration.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

/* =========================================================================
 * L8: Particle Filter for INS/GPS Integration
 *
 * Unlike the extended Kalman filter which assumes Gaussian statistics,
 * a particle filter (sequential Monte Carlo) can handle non-Gaussian
 * distributions and nonlinear dynamics without linearization.
 *
 * SIR (Sampling Importance Resampling) algorithm:
 *   1. Prediction: propagate each particle through the INS dynamics
 *   2. Update: weight each particle by GPS measurement likelihood
 *   3. Resampling: systematic resampling to avoid particle degeneracy
 *
 * Each particle represents a complete navigation state hypothesis.
 * ========================================================================= */

/** Single particle in the particle filter */
typedef struct {
    ins_nav_solution_t nav;   /** Navigation solution */
    double            weight; /** Importance weight */
} ins_particle_t;

/**
 * Initialize N particles around an initial navigation solution.
 *
 * Each particle is initialized by adding Gaussian noise to the
 * initial state: position, velocity, and attitude.
 *
 * @param particles    Output array of N particles
 * @param N            Number of particles
 * @param init         Initial navigation solution (mean)
 * @param sigma_pos    Position std dev for spreading [m]
 * @param sigma_vel    Velocity std dev for spreading [m/s]
 * @param sigma_att    Attitude std dev for spreading [rad]
 */
void ins_pf_init(ins_particle_t *particles, size_t N,
                  const ins_nav_solution_t *init,
                  double sigma_pos, double sigma_vel, double sigma_att) {
    if (!particles || !init || N == 0) return;

    for (size_t i = 0; i < N; i++) {
        memcpy(&particles[i].nav, init, sizeof(ins_nav_solution_t));

        /* Spread position (in meters to lat/lon) */
        double M = ins_meridian_radius(init->pos.lat);
        double N = ins_prime_vertical_radius(init->pos.lat);
        double cos_lat = cos(init->pos.lat);

        double r1 = ((double)rand() / RAND_MAX - 0.5) * 2.0;
        double r2 = ((double)rand() / RAND_MAX - 0.5) * 2.0;
        double r3 = ((double)rand() / RAND_MAX - 0.5) * 2.0;

        particles[i].nav.pos.lat += r1 * sigma_pos / (M + init->pos.alt);
        particles[i].nav.pos.lon += r2 * sigma_pos / ((N + init->pos.alt) * cos_lat);
        particles[i].nav.pos.alt += r3 * sigma_pos;

        /* Spread velocity */
        for (int j = 0; j < 3; j++) {
            double r = ((double)rand() / RAND_MAX - 0.5) * 2.0;
            particles[i].nav.vel_ned.x += r * sigma_vel;
            particles[i].nav.vel_ned.y += r * sigma_vel;
            particles[i].nav.vel_ned.z += r * sigma_vel;
        }

        /* Spread attitude */
        double rx = ((double)rand() / RAND_MAX - 0.5) * sigma_att;
        double ry = ((double)rand() / RAND_MAX - 0.5) * sigma_att;
        double rz = ((double)rand() / RAND_MAX - 0.5) * sigma_att;
        ins_quat_t q_init, q_delta, q_noisy;
        q_init.w = init->q[0]; q_init.x = init->q[1];
        q_init.y = init->q[2]; q_init.z = init->q[3];

        double angle = sqrt(rx * rx + ry * ry + rz * rz);
        if (angle > 1e-10) {
            q_delta.w = cos(angle);
            double s = sin(angle) / angle;
            q_delta.x = rx * s;
            q_delta.y = ry * s;
            q_delta.z = rz * s;
            ins_quat_mul(&q_init, &q_delta, &q_noisy);
        } else {
            ins_quat_copy(&q_init, &q_noisy);
        }
        particles[i].nav.q[0] = q_noisy.w;
        particles[i].nav.q[1] = q_noisy.x;
        particles[i].nav.q[2] = q_noisy.y;
        particles[i].nav.q[3] = q_noisy.z;
        particles[i].weight = 1.0 / (double)N;
    }
}

/**
 * Predict step: propagate each particle through INS dynamics.
 *
 * Each particle is propagated independently using the INS mechanization.
 * Process noise is added to simulate sensor random walk.
 *
 * @param particles  Particle array (in/out)
 * @param N          Number of particles
 * @param imu        IMU measurement for this step
 * @param Q_acc      Accelerometer process noise std [m/s^2]
 * @param Q_gyr      Gyroscope process noise std [rad/s]
 */
void ins_pf_predict(ins_particle_t *particles, size_t N,
                     const ins_imu_sample_t *imu,
                     double Q_acc, double Q_gyr) {
    if (!particles || !imu || N == 0) return;

    for (size_t i = 0; i < N; i++) {
        /* Build noisy IMU measurement */
        ins_imu_sample_t imu_noisy;
        imu_noisy.dt = imu->dt;

        /* Add Gaussian-like process noise (Box-Muller would be better,
           but for simplicity we use uniform noise scaled to std) */
        double s = sqrt(3.0);  /* Uniform -> std conversion */
        double ar = ((double)rand() / RAND_MAX - 0.5) * 2.0 * s * Q_acc;
        double ay = ((double)rand() / RAND_MAX - 0.5) * 2.0 * s * Q_acc;
        double az = ((double)rand() / RAND_MAX - 0.5) * 2.0 * s * Q_acc;
        double gr = ((double)rand() / RAND_MAX - 0.5) * 2.0 * s * Q_gyr;
        double gy = ((double)rand() / RAND_MAX - 0.5) * 2.0 * s * Q_gyr;
        double gz = ((double)rand() / RAND_MAX - 0.5) * 2.0 * s * Q_gyr;

        imu_noisy.accel.x = imu->accel.x + ar;
        imu_noisy.accel.y = imu->accel.y + ay;
        imu_noisy.accel.z = imu->accel.z + az;
        imu_noisy.gyro.x  = imu->gyro.x  + gr;
        imu_noisy.gyro.y  = imu->gyro.y  + gy;
        imu_noisy.gyro.z  = imu->gyro.z  + gz;

        /* Propagate particle navigation state */
        ins_mech_state_t ms;
        ins_mech_init(&ms, INS_MECH_CONING,
                      particles[i].nav.pos.lat,
                      particles[i].nav.pos.lon,
                      particles[i].nav.pos.alt);
        ms.quat.w = particles[i].nav.q[0];
        ms.quat.x = particles[i].nav.q[1];
        ms.quat.y = particles[i].nav.q[2];
        ms.quat.z = particles[i].nav.q[3];
        ins_vec3_copy(&particles[i].nav.vel_ned, &ms.vel_ned);

        ins_mech_step(&ms, &imu_noisy);
        ins_mech_get_solution(&ms, &particles[i].nav);
    }
}

/**
 * Update particle weights based on GPS measurement likelihood.
 *
 * Likelihood: p(z|x) = N(z; h(x), R)
 *
 * Weights are updated multiplicatively: w_i = w_i * p(z|x_i)
 * Then normalized.
 *
 * @param particles  Particle array (in/out)
 * @param N          Number of particles
 * @param gps        GPS measurement
 * @return           Sum of unnormalized weights (used for effective sample size)
 */
double ins_pf_update(ins_particle_t *particles, size_t N,
                      const ins_gps_measurement_t *gps) {
    if (!particles || !gps || N == 0) return 0.0;

    double weight_sum = 0.0;

    for (size_t i = 0; i < N; i++) {
        /* Position innovation */
        double M = ins_meridian_radius(particles[i].nav.pos.lat);
        double N_rad = ins_prime_vertical_radius(particles[i].nav.pos.lat);
        double cos_lat = cos(particles[i].nav.pos.lat);

        double dlat = (gps->pos.lat - particles[i].nav.pos.lat) * (M + particles[i].nav.pos.alt);
        double dlon = (gps->pos.lon - particles[i].nav.pos.lon) * (N_rad + particles[i].nav.pos.alt) * cos_lat;
        double dalt = gps->pos.alt - particles[i].nav.pos.alt;

        /* Velocity innovation */
        double dvn = gps->vel_ned.x - particles[i].nav.vel_ned.x;
        double dve = gps->vel_ned.y - particles[i].nav.vel_ned.y;
        double dvd = gps->vel_ned.z - particles[i].nav.vel_ned.z;

        /* Log Gaussian likelihood (unnormalized) */
        double log_lik = 0.0;
        log_lik -= 0.5 * dlat * dlat / (gps->pos_std[0] * gps->pos_std[0]);
        log_lik -= 0.5 * dlon * dlon / (gps->pos_std[1] * gps->pos_std[1]);
        log_lik -= 0.5 * dalt * dalt / (gps->pos_std[2] * gps->pos_std[2]);
        log_lik -= 0.5 * dvn * dvn / (gps->vel_std[0] * gps->vel_std[0]);
        log_lik -= 0.5 * dve * dve / (gps->vel_std[1] * gps->vel_std[1]);
        log_lik -= 0.5 * dvd * dvd / (gps->vel_std[2] * gps->vel_std[2]);

        /* Convert to weight (avoid zero weight) */
        double lik = exp(log_lik);
        if (lik < 1e-300) lik = 1e-300;

        particles[i].weight *= lik;
        weight_sum += particles[i].weight;
    }

    /* Normalize weights */
    if (weight_sum > 0) {
        for (size_t i = 0; i < N; i++) {
            particles[i].weight /= weight_sum;
        }
    }

    return weight_sum;
}

/**
 * Systematic resampling to avoid particle degeneracy.
 *
 * When the effective sample size N_eff = 1 / sum(w_i^2) is too low,
 * particles are resampled: low-weight particles are replaced by
 * copies of high-weight particles.
 *
 * This implementation uses systematic resampling (Kitagawa, 1996)
 * which has lower variance than multinomial resampling.
 *
 * Reference: Kitagawa (1996), J. Comput. Graph. Stat., 5(1): 1-25.
 *
 * @param particles Particle array (in/out, contents replaced by resampled set)
 * @param N         Number of particles
 * @param N_eff     Output: effective sample size
 */
void ins_pf_resample(ins_particle_t *particles, size_t N, double *N_eff) {
    if (!particles || N == 0) return;

    /* Compute cumulative sum of weights */
    double *cumsum = (double *)malloc(N * sizeof(double));
    if (!cumsum) return;

    cumsum[0] = particles[0].weight;
    double weight_sq_sum = particles[0].weight * particles[0].weight;
    for (size_t i = 1; i < N; i++) {
        cumsum[i] = cumsum[i - 1] + particles[i].weight;
        weight_sq_sum += particles[i].weight * particles[i].weight;
    }

    if (N_eff) {
        *N_eff = 1.0 / weight_sq_sum;
    }

    /* Systematic resampling */
    ins_particle_t *new_particles = (ins_particle_t *)malloc(N * sizeof(ins_particle_t));
    if (!new_particles) { free(cumsum); return; }

    double u0 = ((double)rand() / RAND_MAX) / (double)N;
    size_t j = 0;

    for (size_t i = 0; i < N; i++) {
        double u = u0 + (double)i / (double)N;
        while (u > cumsum[j] && j < N - 1) j++;
        memcpy(&new_particles[i], &particles[j], sizeof(ins_particle_t));
        new_particles[i].weight = 1.0 / (double)N;
    }

    memcpy(particles, new_particles, N * sizeof(ins_particle_t));
    free(new_particles);
    free(cumsum);
}

/* =========================================================================
 * L8: Transfer Alignment (Velocity Matching)
 *
 * Transfer alignment initializes a slave INS in-flight using a master
 * INS (e.g., aircraft INS initializing a missile INS before launch).
 *
 * Velocity matching is the simplest method:
 *   The slave INS is aligned by comparing its velocity output with
 *   the master INS velocity over a short time window (30-120 s).
 *
 * The angular difference between master and slave velocity
 * vectors reveals the heading misalignment.
 *
 * This function estimates the heading error between master and slave
 * from the velocity difference during a maneuver.
 * ========================================================================= */

/**
 * Estimate initial heading error from velocity matching.
 *
 * During a straight and level acceleration, the velocity vectors
 * of master and slave diverge if they have different heading.
 *
 * delta_psi ~ atan2(dV_east, dV_north) for small tilt errors.
 *
 * @param vel_master  Master INS velocity history [num_samples]
 * @param vel_slave   Slave INS velocity history [num_samples]
 * @param num_samples Number of measurements
 * @return            Estimated heading misalignment [rad]
 */
double ins_transfer_align_heading(const ins_vec3_t *vel_master,
                                   const ins_vec3_t *vel_slave,
                                   size_t num_samples) {
    if (!vel_master || !vel_slave || num_samples < 2) return 0.0;

    double sum_dn = 0.0, sum_de = 0.0;

    for (size_t i = 0; i < num_samples; i++) {
        double dv_n = vel_master[i].x - vel_slave[i].x;
        double dv_e = vel_master[i].y - vel_slave[i].y;
        sum_dn += dv_n;
        sum_de += dv_e;
    }

    return atan2(sum_de, sum_dn);
}

/* =========================================================================
 * L8: High-Order Sculling Algorithm (Savage, 1998)
 *
 * While the 2-sample sculling algorithm (in ins_attitude.c) provides
 * basic sculling compensation, higher-order algorithms with more
 * samples per update provide better accuracy in high-vibration
 * environments.
 *
 * 4-sample sculling algorithm:
 *   delta_v_scull = k1*(da1 x dv1) + k2*(da1 x dv2 + da2 x dv1)
 *                 + k3*(da1 x dv3 + da3 x dv1)
 *                 + k4*(da2 x dv2) + k5*(da2 x dv3 + da3 x dv2)
 *                 + k6*(da3 x dv3)
 *
 * where the coefficients k_i are optimized for coning frequency
 * response. The coefficients below use the optimal values from
 * Savage (1998) for 4-sample algorithm.
 * ========================================================================= */

/**
 * 4-sample sculling compensation (Savage, 1998 optimized).
 *
 * @param da          Array of 4 angular increment samples [rad]
 * @param dv          Array of 4 velocity increment samples [m/s]
 * @param dv_scull    Output sculling compensation [m/s]
 */
void ins_sculling_4sample(const ins_vec3_t da[4],
                           const ins_vec3_t dv[4],
                           ins_vec3_t *dv_scull) {
    if (!da || !dv || !dv_scull) return;

    /* Savage (1998) 4-sample optimized coefficients */
    const double k1 =  54.0 / 105.0;
    const double k2 =  92.0 / 105.0;
    const double k3 = 214.0 / 105.0;
    const double k4 =  49.0 / 210.0;
    const double k5 =  16.0 / 35.0;
    const double k6 =   9.0 / 35.0;

    ins_vec3_t term, sum;
    ins_vec3_zero(&sum);

    /* Term 1: k1 * (da0 x dv0) */
    ins_vec3_cross(&da[0], &dv[0], &term);
    ins_vec3_scale(&term, k1, &term);
    ins_vec3_add(&sum, &term, &sum);

    /* Term 2: k2 * (da0 x dv1) */
    ins_vec3_cross(&da[0], &dv[1], &term);
    ins_vec3_scale(&term, k2, &term);
    ins_vec3_add(&sum, &term, &sum);

    /* Term 3: k2 * (da1 x dv0) */
    ins_vec3_cross(&da[1], &dv[0], &term);
    ins_vec3_scale(&term, k2, &term);
    ins_vec3_add(&sum, &term, &sum);

    /* Term 4: k3 * (da0 x dv2) */
    ins_vec3_cross(&da[0], &dv[2], &term);
    ins_vec3_scale(&term, k3, &term);
    ins_vec3_add(&sum, &term, &sum);

    /* Term 5: k3 * (da2 x dv0) */
    ins_vec3_cross(&da[2], &dv[0], &term);
    ins_vec3_scale(&term, k3, &term);
    ins_vec3_add(&sum, &term, &sum);

    /* Term 6: k4 * (da1 x dv1) */
    ins_vec3_cross(&da[1], &dv[1], &term);
    ins_vec3_scale(&term, k4, &term);
    ins_vec3_add(&sum, &term, &sum);

    /* Term 7: k5 * (da1 x dv2) */
    ins_vec3_cross(&da[1], &dv[2], &term);
    ins_vec3_scale(&term, k5, &term);
    ins_vec3_add(&sum, &term, &sum);

    /* Term 8: k5 * (da2 x dv1) */
    ins_vec3_cross(&da[2], &dv[1], &term);
    ins_vec3_scale(&term, k5, &term);
    ins_vec3_add(&sum, &term, &sum);

    /* Term 9: k6 * (da2 x dv2) */
    ins_vec3_cross(&da[2], &dv[2], &term);
    ins_vec3_scale(&term, k6, &term);
    ins_vec3_add(&sum, &term, &sum);

    /* Term 10: da3 x dv0 + da0 x dv3, da3 x dv1 + da1 x dv3, da3 x dv2 + da2 x dv3, da3 x dv3 */
    /* These higher-order terms are available with 3 subsamples per update;
       for 4-sample, we include da3/dv3 cross with earlier samples */
    ins_vec3_cross(&da[3], &dv[0], &term);
    ins_vec3_scale(&term, k1, &term);
    ins_vec3_add(&sum, &term, &sum);

    ins_vec3_cross(&da[0], &dv[3], &term);
    ins_vec3_scale(&term, k1, &term);
    ins_vec3_add(&sum, &term, &sum);

    ins_vec3_cross(&da[3], &dv[1], &term);
    ins_vec3_scale(&term, k2, &term);
    ins_vec3_add(&sum, &term, &sum);

    ins_vec3_cross(&da[1], &dv[3], &term);
    ins_vec3_scale(&term, k2, &term);
    ins_vec3_add(&sum, &term, &sum);

    ins_vec3_cross(&da[3], &dv[2], &term);
    ins_vec3_scale(&term, k5, &term);
    ins_vec3_add(&sum, &term, &sum);

    ins_vec3_cross(&da[2], &dv[3], &term);
    ins_vec3_scale(&term, k5, &term);
    ins_vec3_add(&sum, &term, &sum);

    ins_vec3_cross(&da[3], &dv[3], &term);
    ins_vec3_scale(&term, k6, &term);
    ins_vec3_add(&sum, &term, &sum);

    ins_vec3_copy(&sum, dv_scull);
}

/* =========================================================================
 * L8: Coning Error Analysis
 *
 * The coning error is the dominant source of attitude drift in
 * strapdown INS under vibration. This function computes the
 * net coning error for a given coning motion.
 *
 * Coning motion: angular oscillation about two axes
 *   omega_x(t) = A * omega_c * cos(omega_c * t)
 *   omega_y(t) = A * omega_c * sin(omega_c * t)
 *
 * The net rotation per cycle (false rotation about z-axis) is:
 *   phi_z_per_cycle = pi * A^2  (for small A)
 *
 * Over N cycles: phi_z = N * pi * A^2
 *
 * Reference: Bortz (1971), IEEE Trans. Aerosp. Electron. Syst.
 * ========================================================================= */

/**
 * Compute coning error rate for given vibration parameters.
 *
 * @param amplitude    Coning half-angle amplitude [rad]
 * @param frequency    Coning frequency [rad/s]
 * @return             Net drift rate (false rotation) [rad/s]
 */
double ins_coning_error_rate(double amplitude, double frequency) {
    if (amplitude <= 0 || frequency <= 0) return 0.0;

    /* phi_z per cycle = pi * A^2 */
    double err_per_cycle = M_PI * amplitude * amplitude;
    double cycle_period = 2.0 * M_PI / frequency;

    return err_per_cycle / cycle_period;
}

/**
 * Compute the expected coning error for an uncompensated strapdown
 * algorithm over a given time duration.
 *
 * @param amplitude    Coning half-angle [rad]
 * @param frequency    Coning frequency [rad/s]
 * @param duration     Total time [s]
 * @return             Net coning angle error [rad]
 */
double ins_coning_error_accumulated(double amplitude, double frequency,
                                     double duration) {
    return ins_coning_error_rate(amplitude, frequency) * duration;
}

/* =========================================================================
 * L7: Boeing 747 Cruise Example — INS Accuracy Analysis
 *
 * A typical Boeing 747 cruises at ~900 km/h (250 m/s, Mach 0.85)
 * at altitude ~10,700 m (FL350). With a navigation-grade INS
 * (0.005 deg/hr gyro bias), the position drift during a 10-hour
 * transatlantic flight is:
 *
 *   P_drift_at_10hr = 10 * 0.005 = 0.05 nm/hr * 10 hr = 0.5 nm CEP
 *
 * With GPS and Kalman filter integration, the position accuracy
 * is bounded to < 5 m (95%) throughout the flight.
 * ========================================================================= */

/**
 * Boeing 747 cruise INS accuracy analysis.
 *
 * Computes time-dependent position, velocity, and attitude drift
 * for a navigation-grade INS during a long-haul flight.
 *
 * @param flight_time    Flight duration [s]
 * @param pos_drift      Output: 1-sigma position drift [m]
 * @param vel_drift      Output: 1-sigma velocity drift [m/s]
 * @param att_drift      Output: 1-sigma attitude drift [rad]
 */
void ins_boeing747_ins_analysis(double flight_time,
                                 double *pos_drift,
                                 double *vel_drift,
                                 double *att_drift) {
    /* Navigation-grade INS parameters (Boeing 747 typically uses Honeywell HG9900 or similar) */
    ins_imu_error_model_t model;

    /* Gyro errors: 0.005 deg/hr bias, 0.001 deg/sqrt(hr) ARW */
    double gyro_bias_rads = 0.005 * (M_PI / 180.0) / 3600.0;
    double gyro_arw_radsh = 0.001 * (M_PI / 180.0) / 60.0;

    model.gyro_x.bias_offset = gyro_bias_rads;
    model.gyro_y.bias_offset = gyro_bias_rads;
    model.gyro_z.bias_offset = gyro_bias_rads;
    model.gyro_x.white_noise_std = gyro_arw_radsh;
    model.gyro_y.white_noise_std = gyro_arw_radsh;
    model.gyro_z.white_noise_std = gyro_arw_radsh;

    /* Accel errors: 25 micro-g bias, 0.001 m/s/sqrt(hr) VRW */
    double accel_bias_ms2 = 25e-6 * INS_GRAVITY_EQUATOR;
    double accel_vrw = 0.001 / 60.0;

    model.accel_x.bias_offset = accel_bias_ms2;
    model.accel_y.bias_offset = accel_bias_ms2;
    model.accel_z.bias_offset = accel_bias_ms2;
    model.accel_x.white_noise_std = accel_vrw;
    model.accel_y.white_noise_std = accel_vrw;
    model.accel_z.white_noise_std = accel_vrw;

    /* Zero out the rest */
    memset(&model.accel_x.bias_instability, 0, sizeof(model) - 6 * sizeof(ins_axis_error_t));

    ins_error_predict_drift(&model, flight_time, pos_drift, vel_drift, att_drift);
}
