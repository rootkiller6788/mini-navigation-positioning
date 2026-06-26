# Knowledge Graph — mini-gnss-gps

## L1: Definitions

| Item | C Struct/Type | Lean Formalization | Status |
|------|---------------|-------------------|--------|
| WGS84 semi-major axis (a) | `GNSS_WGS84_A` macro | `wgs84_a_km : Int` | ✅ |
| WGS84 flattening (f) | `GNSS_WGS84_F` macro | `wgs84_b_km` comparison | ✅ |
| Eccentricity (e²) | `GNSS_WGS84_E2` macro | — | ✅ |
| Speed of light (c) | `GNSS_C_LIGHT` | `cLightMps : Nat` | ✅ |
| L1/L2/L5 frequencies | `GNSS_L1_FREQ` etc. | `l1FrequencyMHz : Nat` etc. | ✅ |
| ECEF coordinates | `gnss_ecef_t` | `Vec3` structure | ✅ |
| LLA coordinates | `gnss_lla_t` | — | ✅ |
| ENU coordinates | `gnss_enu_t` | — | ✅ |
| GPS time (week + SOW) | `gnss_gpstime_t` | `GPSTime` structure | ✅ |
| UTC time | `gnss_utctime_t` | — | ✅ |
| C/A code (1023 chips) | `gnss_ca_code_t` | `caCodeLength : Nat` | ✅ |
| PRN code | `gnss_lfsr_state_t` | C/A code properties | ✅ |
| Pseudorange (P) | `gnss_pseudorange_t` | `Pseudorange` structure | ✅ |
| Carrier phase (Φ) | `gnss_carrier_phase_t` | — | ✅ |
| Integer ambiguity (N) | `gnss_ambiguity_t` | — | ✅ |
| DOP (GDOP/PDOP/HDOP/VDOP/TDOP) | `gnss_dop_t` | `DOPRating` inductive | ✅ |
| Q-matrix | `gnss_qmatrix_t` | — | ✅ |
| Ephemeris (Keplerian) | `gnss_ephemeris_t` | — | ✅ |
| Almanac | `gnss_almanac_t` | — | ✅ |
| Ionospheric delay (I) | `gnss_klobuchar_params_t` | Ionospheric freq scaling | ✅ |
| Tropospheric delay (T) | `gnss_tropo_params_t` | — | ✅ |
| Multipath (M) | `gnss_multipath_t` | — | ✅ |
| UERE | `gnss_error_budget_t` | — | ✅ |
| Constellation types | `gnss_constellation_t` | — | ✅ |
| Frequency bands | `gnss_band_t` | — | ✅ |
| PVT solution | `gnss_pvt_solution_t` | — | ✅ |
| Weight strategy | `gnss_weight_strategy_t` | — | ✅ |
| Hatch filter | `gnss_hatch_filter_t` | Hatch steady-state | ✅ |
| DGPS correction | `gnss_dgps_correction_t` | — | ✅ |
| Doppler shift | `gnss_doppler_t` | — | ✅ |

**L1 Status: Complete ✅** (29 independent definitions)

## L2: Core Concepts

| Concept | Implementation | Status |
|---------|---------------|--------|
| Trilateration (TOA ranging) | `gnss_ls_position_solve()` | ✅ |
| Gold code generation (G1/G2 LFSR) | `gnss_ca_code_generate()` | ✅ |
| Satellite ephemeris computation | `gnss_satpos_from_ephemeris()` | ✅ |
| Kepler's equation (E - e·sinE = M) | `gnss_kepler_solve()` | ✅ |
| Ionospheric advance/delay | `gnss_iono_klobuchar()` | ✅ |
| Tropospheric delay | `gnss_tropo_saastamoinen()` | ✅ |
| Sagnac effect | `gnss_sagnac_correction()` | ✅ |
| Code-minus-carrier (CMC) | `gnss_cmc_compute()` | ✅ |
| Cycle slip detection | `gnss_slip_detect()` | ✅ |
| Double-differencing | `gnss_ambiguity_dd_t` | ✅ |
| RAIM integrity monitoring | `gnss_raim_fde()` | ✅ |
| Relativistic clock correction | `gnss_sat_clock_correction()` | ✅ |
| Doppler-based velocity | `gnss_velocity_solve()` | ✅ |
| Satellite elevation/azimuth | `gnss_sat_elevation()` | ✅ |
| Elevation-dependent weighting | `gnss_weight_strategy_t` | ✅ |
| Satellite geometry (DOP) | `gnss_dop_compute()` | ✅ |
| Ionosphere-free combination | `gnss_ionofree_combination()` | ✅ |
| Wide-lane combination | `gnss_mw_combination()` | ✅ |

**L2 Status: Complete ✅** (18 core concepts)

## L3: Mathematical Structures

| Structure | C Implementation | Lean Formalization | Status |
|-----------|-----------------|-------------------|--------|
| 3-vector operations | `gnss_vec3_*()` (7 functions) | `Vec3` + operations | ✅ |
| 3×3 matrix | `gnss_mat33_t` | — | ✅ |
| 4×4 matrix | `gnss_mat44_t` | — | ✅ |
| Cholesky decomposition | `gnss_mat44_cholesky()` | — | ✅ |
| Gauss-Newton iteration | `gnss_ls_position_solve()` | — | ✅ |
| Least squares normal eqns | `gnss_normal_eqn_solve()` | Normal equation counting | ✅ |
| Minkowski inner product | `gnss_bancroft_solve()` | Lorentz metric | ✅ |
| Rotation matrix (ECEF→ENU) | `gnss_ecef_to_enu()` | `rot2` rotation | ✅ |
| Coordinate transformations | ECEF↔LLA↔ENU | — | ✅ |
| Design matrix (Jacobian) | `gnss_design_matrix_t` | `DesignMatrix` structure | ✅ |
| Weighted least squares | `gnss_wls_position_solve()` | — | ✅ |
| Covariance propagation | DOP from Q-matrix | — | ✅ |

**L3 Status: Complete ✅** (12 mathematical structures)

## L4: Fundamental Laws/Theorems

| Law/Theorem | C Verification | Lean Statement | Status |
|-------------|---------------|---------------|--------|
| Kepler's Third Law (T²∝a³) | Satellite position computation | Orbital period approximates half-day | ✅ |
| Newton's law of gravitation | In `gnss_satpos_from_ephemeris()` | — | ✅ |
| Speed of light (c constant) | `GNSS_C_LIGHT` | `cLightMps` as constant | ✅ |
| Pseudorange equation | `gnss_range_corrections_t` | `Pseudorange` structure | ✅ |
| GPS timing (P=ρ+c·Δt) | Position solution residual | Timing precision lemma | ✅ |
| Ionospheric dispersion (∝1/f²) | Klobuchar model | `iono_scaling_f1_gt_f2` theorem | ✅ |
| DOP theory (σ=PDOP×σ_UERE) | `gnss_dop_t` | `DOPRating` classification | ✅ |
| Gold code correlation bound | Cross-correlation test | `gold_bound_lt_one` theorem | ✅ |
| Doppler shift (Δf=-f·v/c) | `gnss_doppler_compute()` | — | ✅ |
| Earth rotation (Sagnac) | `gnss_sagnac_correction()` | — | ✅ |
| Ambiguity integer property | `gnss_ambiguity_dd_t` | — | ✅ |

**L4 Status: Complete ✅** (11 fundamental laws/theorems with code verification)

## L5: Algorithms/Methods

| Algorithm | Implementation | Complexity | Status |
|-----------|---------------|------------|--------|
| Bancroft direct solution | `gnss_bancroft_solve()` | O(n) | ✅ |
| Iterative least squares PVT | `gnss_ls_position_solve()` | O(n·iter) | ✅ |
| Weighted least squares | `gnss_wls_position_solve()` | O(n·iter) | ✅ |
| Klobuchar ionospheric model | `gnss_iono_klobuchar()` | O(1) | ✅ |
| Saastamoinen tropospheric model | `gnss_tropo_saastamoinen()` | O(1) | ✅ |
| Hopfield tropospheric model | `gnss_tropo_hopfield_total()` | O(1) | ✅ |
| Hatch carrier smoothing | `gnss_hatch_smooth()` | O(1) per epoch | ✅ |
| Newton-Raphson (Kepler) | `gnss_kepler_solve()` | O(iter) | ✅ |
| C/A code generation (LFSR) | `gnss_ca_code_generate()` | O(1023) | ✅ |
| Satellite position (ephemeris) | `gnss_satpos_from_ephemeris()` | O(20) per sat | ✅ |
| RAIM FDE | `gnss_raim_fde()` | O(n²) | ✅ |
| Best subset selection | `gnss_select_best_subset()` | O(n²·inv) | ✅ |
| Cycle slip detection (GF+Doppler) | `gnss_slip_detect()` | O(1) per sat | ✅ |
| CMC analysis | `gnss_cmc_compute()` | O(n) | ✅ |
| Gauss-Jordan matrix inversion | `gnss_mat44_inverse()` | O(1) for 4×4 | ✅ |
| Cholesky decomposition | `gnss_mat44_cholesky()` | O(1) for 4×4 | ✅ |
| Bowring geodetic iteration | `gnss_ecef_to_lla()` | O(10) iters | ✅ |

**L5 Status: Complete ✅** (17 algorithms, all with full implementations)

## L6: Canonical Problems

| Problem | Example/Test | Status |
|---------|-------------|--------|
| Single-point positioning (SPP) | `example_position_solve.c` | ✅ |
| DOP analysis & constellation selection | `example_dop_calc.c` | ✅ |
| Carrier-smoothed pseudorange | `example_carrier_smoothing.c` | ✅ |
| C/A code acquisition (correlation) | `test_gnss.c` (test_ca_code) | ✅ |
| Bancroft initialization + LS refinement | `test_gnss.c` (test_bancroft + test_ls_position) | ✅ |
| Atmospheric correction application | `test_gnss.c` (test_klobuchar + test_saastamoinen) | ✅ |
| Cycle slip detection & recovery | `test_gnss.c` (test_cycle_slip_detection + example) | ✅ |

**L6 Status: Complete ✅** (7 canonical problems with full examples/tests)

## L7: Applications

| Application | Reference | Status |
|-------------|-----------|--------|
| GPS SPS accuracy (≤9m horizontal 95%) | `example_position_solve.c` | ✅ |
| Aviation approach (DOP prediction) | `example_dop_calc.c` | ✅ |
| High-precision surveying (carrier smoothing) | `example_carrier_smoothing.c` | ✅ |
| Smartphone A-GPS (TTFF reduction) | Lean `agps_faster` theorem | ✅ |
| Urban canyon simulation (multi-sat) | `example_position_solve.c` (8 sats at SFO) | ✅ |
| Multi-constellation GNSS | Lean `multi_constellation_advantage` | ✅ |

**L7 Status: Complete ✅** (6 real-world applications)

## L8: Advanced Topics

| Topic | Implementation | Status |
|-------|---------------|--------|
| Ionosphere-free dual-frequency | `gnss_ionofree_combination()` | ✅ |
| Melbourne-Wübbena wide-lane | `gnss_mw_combination()` | ✅ |
| RAIM integrity monitoring | `gnss_raim_fde()` | ✅ |
| Best constellation subset selection | `gnss_select_best_subset()` | ✅ |
| DGPS differential corrections | `gnss_dgps_apply()` | ✅ |

**L8 Status: Partial ✅** (5/10 advanced topics)

## L9: Research Frontiers

| Topic | Documentation | Implementation | Status |
|-------|--------------|----------------|--------|
| Multi-constellation GNSS fusion | Referenced in Lean | Partial | Partial |
| 6G positioning integration | Documented in gap-report | Not implemented | Partial |
| Quantum navigation concepts | Documented | Not implemented | Partial |

**L9 Status: Partial** (documented, minimal implementation)
