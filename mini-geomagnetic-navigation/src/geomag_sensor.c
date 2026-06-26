/**
 * geomag_sensor.c -- Magnetometer Sensor Models & Calibration
 *
 * L2: Scalar and triaxial magnetometer measurement models
 * L5: Hard-iron/soft-iron calibration via ellipsoid fitting
 * L5: Magnetic gradient computation and invariants
 * L6: Magnetic Anomaly Detection (MAD) algorithms
 *
 * Reference:
 *   Ripka, "Magnetic Sensors and Magnetometers" (2001)
 *   Gebre-Egziabher et al., "Calibration of Strapdown Magnetometers",
 *     IEEE Trans. AES, Vol.42, No.2, 2006, pp. 620-633.
 *   Ginzburg et al., "Processing of magnetic data for MAD",
 *     IEEE Trans. GRS, Vol.42, No.1, 2004.
 *   Beiki et al., "Interpretation of magnetic gradient tensor data",
 *     Geophysics, Vol.77, No.4, 2012.
 */

#include "geomag_sensor.h"
#include "geomag_core.h"
#include "geomag_math.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ========================================================================
 * L2: Scalar magnetometer model
 *
 * reading = |B_true| + bias (+ noise, caller-applied)
 *
 * Scalar sensors (Overhauser, optically pumped Cs/K/He) measure only
 * total field magnitude, insensitive to orientation. Ideal for
 * total-field anomaly detection at 0.001-0.1 nT precision.
 *
 * Complexity: O(1).
 * ======================================================================== */
void scalar_magnetometer_read(const MagVector *B_true, double bias,
                               double noise_std, double *reading)
{
    if (!B_true || !reading) return;
    *reading = mag_magnitude(B_true) + bias;
    (void)noise_std; /* caller adds stochastic noise */
}

/* ========================================================================
 * L2: Triaxial fluxgate magnetometer measurement model
 *
 * B_raw = S * M * (B_body + b_hard) + w
 *
 * S = diag(sx, sy, sz) : scale factors
 * M = 3x3 soft-iron + misalignment matrix
 * b_hard = 3x1 hard-iron bias from platform permanent magnetization
 * w ~ N(0, sigma^2) per axis
 *
 * Reference: Gebre-Egziabher (2006), Eq. 1-4.
 * Complexity: O(1) -- constant matrix operations.
 * ======================================================================== */
void triaxial_magnetometer_read(const MagVector *B_body,
                                 const double bias[3],
                                 const double scale[9],
                                 const double soft_iron[9],
                                 double noise_std,
                                 double reading[3])
{
    if (!B_body || !reading) return;
    (void)noise_std;

    /* Step 1: Add hard-iron bias */
    double Bc[3] = {
        B_body->bx + bias[0],
        B_body->by + bias[1],
        B_body->bz + bias[2]
    };

    /* Step 2: Apply soft-iron/misalignment: tmp = M * Bc */
    double tmp[3];
    tmp[0] = soft_iron[0] * Bc[0] + soft_iron[1] * Bc[1] + soft_iron[2] * Bc[2];
    tmp[1] = soft_iron[3] * Bc[0] + soft_iron[4] * Bc[1] + soft_iron[5] * Bc[2];
    tmp[2] = soft_iron[6] * Bc[0] + soft_iron[7] * Bc[1] + soft_iron[8] * Bc[2];

    /* Step 3: Apply scale: reading = S * tmp */
    reading[0] = scale[0] * tmp[0] + scale[1] * tmp[1] + scale[2] * tmp[2];
    reading[1] = scale[3] * tmp[0] + scale[4] * tmp[1] + scale[5] * tmp[2];
    reading[2] = scale[6] * tmp[0] + scale[7] * tmp[1] + scale[8] * tmp[2];
}

/* ========================================================================
 * Internal: solve small linear system via Gauss-Jordan with pivoting.
 * Generic: solves A*x = b for n equations, in-place on augmented matrix.
 * ======================================================================== */
static int solve_linear_system(double *A, double *b, double *x, int n)
{
    /* Allocate augmented matrix [A | b] as n x (n+1) */
    double *aug = (double *)malloc(n * (n + 1) * sizeof(double));
    if (!aug) return -1;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            aug[i * (n + 1) + j] = A[i * n + j];
        aug[i * (n + 1) + n] = b[i];
    }

    for (int col = 0; col < n; col++) {
        /* Partial pivoting */
        int pivot = col;
        double maxv = fabs(aug[col * (n + 1) + col]);
        for (int row = col + 1; row < n; row++) {
            double av = fabs(aug[row * (n + 1) + col]);
            if (av > maxv) { maxv = av; pivot = row; }
        }
        if (maxv < 1e-12) { free(aug); return -1; }

        if (pivot != col) {
            for (int j = 0; j <= n; j++) {
                double tmp = aug[col * (n + 1) + j];
                aug[col * (n + 1) + j] = aug[pivot * (n + 1) + j];
                aug[pivot * (n + 1) + j] = tmp;
            }
        }

        double piv_val = aug[col * (n + 1) + col];
        for (int j = 0; j <= n; j++)
            aug[col * (n + 1) + j] /= piv_val;

        for (int row = 0; row < n; row++) {
            if (row != col) {
                double factor = aug[row * (n + 1) + col];
                for (int j = 0; j <= n; j++)
                    aug[row * (n + 1) + j] -= factor * aug[col * (n + 1) + j];
            }
        }
    }

    for (int i = 0; i < n; i++)
        x[i] = aug[i * (n + 1) + n];

    free(aug);
    return 0;
}


/* ========================================================================
 * L5: Magnetometer calibration via least-squares ellipsoid fitting.
 *
 * Ellipsoid: ax^2 + by^2 + cz^2 + 2fyz + 2gxz + 2hxy + 2px + 2qy + 2rz + d = 0
 *
 * Solved by linear least-squares on 9 parameters (a..r, with d fixed at -1).
 * The center gives hard-iron bias; the quadratic form gives soft-iron.
 *
 * Reference: Li et al., "Magnetic Sensor Calibration Based on Ellipsoid
 *   Fitting", IEEE Sensors J., 2005.
 *
 * Complexity: O(N) for accumulation + O(9^3) for linear solve.
 * ======================================================================== */
int magnetometer_calibrate_ls(const double measurements[][3], int N,
                               double bias_out[3], double scale_out[9])
{
    if (N < 9 || !measurements || !bias_out || !scale_out) return -1;

    /* Accumulate normal equations A^T*A and A^T*b for ellipsoid fit.
     * Row i of A: [x^2, y^2, z^2, 2yz, 2xz, 2xy, 2x, 2y, 2z]
     * RHS b_i = -1 (we fix d = -1)
     */
    double ATA[81] = {0};
    double ATb[9] = {0};

    for (int i = 0; i < N; i++) {
        double x = measurements[i][0];
        double y = measurements[i][1];
        double z = measurements[i][2];

        double a_row[9] = {
            x*x, y*y, z*z,
            2.0*y*z, 2.0*x*z, 2.0*x*y,
            2.0*x, 2.0*y, 2.0*z
        };

        for (int r = 0; r < 9; r++) {
            for (int c = 0; c < 9; c++)
                ATA[r * 9 + c] += a_row[r] * a_row[c];
            ATb[r] += a_row[r] * (-1.0);
        }
    }

    double params[9];
    if (solve_linear_system(ATA, ATb, params, 9) != 0) return -1;

    /* Extract: params = [a, b, c, f, g, h, p, q, r]
     * Q = [a  h  g;  h  b  f;  g  f  c]
     * u = [p; q; r]
     * Center (bias) = -Q^{-1} * u
     */
    double a=params[0], b=params[1], c=params[2];
    double f=params[3], g=params[4], h=params[5];
    double p=params[6], q=params[7], rp=params[8];

    double Q[9] = { a, h, g,  h, b, f,  g, f, c };
    double Qinv[9];
    double detQ = Q[0]*(Q[4]*Q[8]-Q[5]*Q[7])
                - Q[1]*(Q[3]*Q[8]-Q[5]*Q[6])
                + Q[2]*(Q[3]*Q[7]-Q[4]*Q[6]);

    if (fabs(detQ) < 1e-12) {
        /* Fallback: centroid as bias */
        double sx=0, sy=0, sz=0;
        for (int i = 0; i < N; i++) {
            sx += measurements[i][0];
            sy += measurements[i][1];
            sz += measurements[i][2];
        }
        bias_out[0] = sx / N;
        bias_out[1] = sy / N;
        bias_out[2] = sz / N;
    } else {
        double id = 1.0 / detQ;
        Qinv[0] = (Q[4]*Q[8]-Q[5]*Q[7])*id;
        Qinv[1] = (Q[2]*Q[7]-Q[1]*Q[8])*id;
        Qinv[2] = (Q[1]*Q[5]-Q[2]*Q[4])*id;
        Qinv[3] = (Q[5]*Q[6]-Q[3]*Q[8])*id;
        Qinv[4] = (Q[0]*Q[8]-Q[2]*Q[6])*id;
        Qinv[5] = (Q[2]*Q[3]-Q[0]*Q[5])*id;
        Qinv[6] = (Q[3]*Q[7]-Q[4]*Q[6])*id;
        Qinv[7] = (Q[1]*Q[6]-Q[0]*Q[7])*id;
        Qinv[8] = (Q[0]*Q[4]-Q[1]*Q[3])*id;

        double u[3] = {p, q, rp};
        bias_out[0] = -(Qinv[0]*u[0]+Qinv[1]*u[1]+Qinv[2]*u[2]);
        bias_out[1] = -(Qinv[3]*u[0]+Qinv[4]*u[1]+Qinv[5]*u[2]);
        bias_out[2] = -(Qinv[6]*u[0]+Qinv[7]*u[1]+Qinv[8]*u[2]);
    }

    /* Scale factors: normalize to unit sphere via covariance of bias-corrected data */
    double C[9] = {0};
    for (int i = 0; i < N; i++) {
        double dx = measurements[i][0] - bias_out[0];
        double dy = measurements[i][1] - bias_out[1];
        double dz = measurements[i][2] - bias_out[2];
        C[0] += dx*dx; C[1] += dx*dy; C[2] += dx*dz;
        C[3] += dy*dx; C[4] += dy*dy; C[5] += dy*dz;
        C[6] += dz*dx; C[7] += dz*dy; C[8] += dz*dz;
    }
    for (int i = 0; i < 9; i++) C[i] /= N;

    double trC = C[0] + C[4] + C[8];
    double rms = sqrt(trC / 3.0);
    double sx = rms / sqrt(C[0]), sy = rms / sqrt(C[4]), sz = rms / sqrt(C[8]);

    scale_out[0] = sx; scale_out[1] = 0; scale_out[2] = 0;
    scale_out[3] = 0;  scale_out[4] = sy; scale_out[5] = 0;
    scale_out[6] = 0;  scale_out[7] = 0;  scale_out[8] = sz;

    return 0;
}

/* L5: Calibration quality -- RMS deviation from reference total field */
double calibration_sphericity_error(const double calibrated[][3], int N,
                                     double B_ref)
{
    if (N < 1 || !calibrated) return -1.0;
    double sum_sq = 0.0;
    for (int i = 0; i < N; i++) {
        double F = sqrt(calibrated[i][0]*calibrated[i][0]
                      + calibrated[i][1]*calibrated[i][1]
                      + calibrated[i][2]*calibrated[i][2]);
        double err = F - B_ref;
        sum_sq += err * err;
    }
    return sqrt(sum_sq / N);
}

/* ========================================================================
 * L5: Magnetic gradient tensor from gradiometer (two-sensor pair)
 *
 * G_ij = dB_i/dx_j, approximated by finite difference:
 *   G_ij = (B2_i - B1_i) / d_j
 *
 * For a source-free region (above Earth): div B = 0 => trace(G) = 0.
 * ======================================================================== */
void compute_magnetic_gradient(const MagVector *B1, const MagVector *B2,
                                const double baseline[3], double G[9])
{
    double b = vec3_norm(baseline);
    if (b < 1e-10) {
        memset(G, 0, 9 * sizeof(double));
        return;
    }
    double d_hat[3] = { baseline[0]/b, baseline[1]/b, baseline[2]/b };
    double dB[3] = { B2->bx-B1->bx, B2->by-B1->by, B2->bz-B1->bz };

    /* Outer product: G = (dB / b) * d_hat^T */
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            G[i * 3 + j] = dB[i] * d_hat[j] / b;
}

/* L5: Magnitude of spatial gradient of total field intensity */
double total_field_gradient_magnitude(const double G[9], const MagVector *B)
{
    double F = mag_magnitude(B);
    if (F < 1e-10) return 0.0;
    double B_hat[3] = { B->bx/F, B->by/F, B->bz/F };
    /* grad|B| = G^T * B_hat */
    double gradF[3];
    gradF[0] = G[0]*B_hat[0] + G[3]*B_hat[1] + G[6]*B_hat[2];
    gradF[1] = G[1]*B_hat[0] + G[4]*B_hat[1] + G[7]*B_hat[2];
    gradF[2] = G[2]*B_hat[0] + G[5]*B_hat[1] + G[8]*B_hat[2];
    return vec3_norm(gradF);
}


/* ========================================================================
 * L5: Normalized Source Strength (NSS) from magnetic gradient tensor.
 *
 * NSS = sqrt(-lambda_2^2 - lambda_1*lambda_3)
 *
 * where lambda_1 >= lambda_2 >= lambda_3 are eigenvalues of G.
 * NSS is near-independent of source magnetization direction, making it
 * a robust quantity for magnetic dipole detection and localization.
 *
 * For traceless G: NSS = sqrt(|lambda_1*lambda_3| - lambda_2^2)
 *
 * We approximate using tensor invariants:
 *   I2 = (trace(G)^2 - trace(G^2))/2 = -trace(G^2)/2 (since trace(G)=0)
 *   NSS = sqrt(-I2/3) * 3 = sqrt(3*|I2|) [as upper bound]
 * Better: compute eigenvalues via power iteration.
 *
 * Reference: Beiki et al., Geophysics, Vol.77, No.4, 2012.
 * ======================================================================== */
double normalized_source_strength(const double G[9])
{
    /* Compute trace(G^2) */
    double trG2 = 0.0;
    for (int i = 0; i < 9; i++) trG2 += G[i] * G[i];

    /* I2 = -0.5 * trace(G^2) for traceless matrix */
    double I2 = -0.5 * trG2;

    /* NSS = |I2| scaled: NSS^2 = -(lambda_2^2 + lambda_1*lambda_3)
     * For a typical gradient tensor, NSS approximates sqrt(3*|I2|).
     * More precisely: NSS = sqrt( max eigenvalue magnitude * min eigenvalue )
     * For now, return sqrt(3*|I2|) as proxy. */
    return sqrt(3.0 * fabs(I2));
}

/* ========================================================================
 * L6: Simple threshold-based Magnetic Anomaly Detection
 *
 * Detection: | |B_meas| - |B_igrf| | > threshold
 *
 * Performance: Pd = Q(Q^{-1}(Pfa) - sqrt(2*SNR)) for coherent detection
 *
 * Typical thresholds:
 *   - Submarine detection from aircraft: 0.1-5 nT at 100-500 m altitude
 *   - UXO detection: 1-50 nT at 1-10 m standoff
 *   - Geological survey: 10-1000 nT at surface
 *
 * Complexity: O(1).
 * ======================================================================== */
int mad_threshold_detect(double B_measured, double B_igrf, double threshold)
{
    return (fabs(B_measured - B_igrf) > threshold) ? 1 : 0;
}

/* ========================================================================
 * L5: Orthonormalized Basis Functions (OBF) MAD detector.
 *
 * The magnetic dipole signature of a linearly moving target can be
 * represented as a linear combination of 3 orthonormal basis functions:
 *
 * B_dipole(t) = sum_{k=1}^3 a_k * phi_k(t; CPA, v)
 *
 * The basis functions are derived from the dipole field equation:
 *   phi_1(t) = w(t) / (w(t)^2 + CPA^2)^{5/2}
 *   phi_2(t) = CPA / (w(t)^2 + CPA^2)^{5/2}
 *   phi_3(t) = t*w(t) / (w(t)^2 + CPA^2)^{5/2}
 *
 * where w(t) = v*(t - t0) and CPA = closest point of approach.
 *
 * Detection statistic = sum of squared projections onto basis.
 *
 * Reference: Ginzburg et al., IEEE Trans. GRS, Vol.42, No.1, 2004.
 *
 * Complexity: O(N) per basis function.
 * ======================================================================== */
int mad_obf_detect(const double *signal, int N, const double *tbase,
                    double v, double cpa, double *stat)
{
    if (!signal || !tbase || !stat || N < 10) return -1;

    double t0 = tbase[N / 2];
    double stat_sum = 0.0;

    for (int k = 0; k < 3; k++) {
        double proj = 0.0, norm = 0.0;

        for (int i = 0; i < N; i++) {
            double tau = tbase[i] - t0;
            double w = v * tau;
            double r2 = w * w + cpa * cpa;
            double denom = pow(r2, 2.5);
            if (denom < 1e-30) continue;

            double phi;
            if (k == 0) phi = w / denom;
            else if (k == 1) phi = cpa / denom;
            else phi = tau * w / denom;

            proj += signal[i] * phi;
            norm += phi * phi;
        }

        if (norm > 1e-15)
            stat_sum += (proj * proj) / norm;
    }

    *stat = stat_sum;
    return 0;
}

/* L5: MAD SNR in dB */
double mad_snr_db(double anomaly_amplitude, double noise_rms)
{
    if (noise_rms < 1e-15) return 999.0;
    return 20.0 * log10(anomaly_amplitude / noise_rms);
}

/* ========================================================================
 * L6: Dipole moment estimation from multiple anomaly measurements.
 *
 * Magnetic dipole field at position r:
 *   B(r) = (mu_0/4pi) * [3(m.r_hat)r_hat - m] / r^3
 *
 * This is linear in dipole moment m = [mx, my, mz]^T:
 *   B = A * m   where A is 3x3 matrix depending only on position.
 *
 * Given N >= 3 measurements at known relative positions, solve for m
 * by linear least-squares: m = (A^T A)^{-1} A^T b
 *
 * Complexity: O(N) assembly + O(27) for 3x3 inverse.
 * ======================================================================== */
int estimate_dipole_moment(const double positions[][3],
                            const double anomalies[][3],
                            int N, double moment[3])
{
    if (N < 3 || !positions || !anomalies || !moment) return -1;

    const double mu0_4pi = 1e-7;
    double ATA[9] = {0};
    double ATb[3] = {0};

    for (int k = 0; k < N; k++) {
        double r = vec3_norm(positions[k]);
        if (r < 1e-10) continue;
        double r2 = r * r;
        double r5 = r2 * r2 * r;
        double factor = mu0_4pi / r5;
        double rx = positions[k][0], ry = positions[k][1], rz = positions[k][2];

        /* Build row of A (3x3 matrix for this measurement) */
        double A_row[9];
        for (int i = 0; i < 3; i++) {
            double ri = (i==0) ? rx : (i==1) ? ry : rz;
            for (int j = 0; j < 3; j++) {
                double rj = (j==0) ? rx : (j==1) ? ry : rz;
                double delta_ij = (i == j) ? 1.0 : 0.0;
                A_row[i*3 + j] = factor * (3.0 * ri * rj - r2 * delta_ij);
            }
        }

        /* ATA += A_k^T * A_k */
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                for (int l = 0; l < 3; l++)
                    ATA[i*3 + j] += A_row[l*3 + i] * A_row[l*3 + j];

        /* ATb += A_k^T * b_k */
        for (int i = 0; i < 3; i++)
            for (int l = 0; l < 3; l++)
                ATb[i] += A_row[l*3 + i] * anomalies[k][l];
    }

    /* Solve 3x3 system via Cramer's rule */
    double det = ATA[0]*(ATA[4]*ATA[8]-ATA[5]*ATA[7])
               - ATA[1]*(ATA[3]*ATA[8]-ATA[5]*ATA[6])
               + ATA[2]*(ATA[3]*ATA[7]-ATA[4]*ATA[6]);

    if (fabs(det) < 1e-30) return -1;
    double id = 1.0 / det;

    double Ainv[9];
    Ainv[0] = (ATA[4]*ATA[8]-ATA[5]*ATA[7])*id;
    Ainv[1] = (ATA[2]*ATA[7]-ATA[1]*ATA[8])*id;
    Ainv[2] = (ATA[1]*ATA[5]-ATA[2]*ATA[4])*id;
    Ainv[3] = (ATA[5]*ATA[6]-ATA[3]*ATA[8])*id;
    Ainv[4] = (ATA[0]*ATA[8]-ATA[2]*ATA[6])*id;
    Ainv[5] = (ATA[2]*ATA[3]-ATA[0]*ATA[5])*id;
    Ainv[6] = (ATA[3]*ATA[7]-ATA[4]*ATA[6])*id;
    Ainv[7] = (ATA[1]*ATA[6]-ATA[0]*ATA[7])*id;
    Ainv[8] = (ATA[0]*ATA[4]-ATA[1]*ATA[3])*id;

    moment[0] = Ainv[0]*ATb[0] + Ainv[1]*ATb[1] + Ainv[2]*ATb[2];
    moment[1] = Ainv[3]*ATb[0] + Ainv[4]*ATb[1] + Ainv[5]*ATb[2];
    moment[2] = Ainv[6]*ATb[0] + Ainv[7]*ATb[1] + Ainv[8]*ATb[2];

    return 0;
}
