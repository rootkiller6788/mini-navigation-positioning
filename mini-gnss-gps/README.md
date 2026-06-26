# mini-gnss-gps

Global Navigation Satellite System (GNSS) & GPS receiver fundamentals:
signal structure, pseudorange modeling, position solution algorithms,
DOP analysis, carrier phase processing, and atmospheric corrections.
Implemented in C with Lean 4 formalization.

## Module Status: COMPLETE ✅

- **L1-L6**: Complete
- **L7**: Complete (6 applications: SPS accuracy, aviation DOP, surveying,
  urban canyon, A-GPS, multi-GNSS)
- **L8**: Partial (5/10 advanced topics: IF combination, MW wide-lane,
  RAIM, subset selection, DGPS)
- **L9**: Partial (documented: LEO-PNT, quantum positioning, 6G positioning)

| Level | Status | Score |
|-------|--------|-------|
| L1 Definitions | Complete | 2 |
| L2 Core Concepts | Complete | 2 |
| L3 Math Structures | Complete | 2 |
| L4 Fundamental Laws | Complete | 2 |
| L5 Algorithms/Methods | Complete | 2 |
| L6 Canonical Problems | Complete | 2 |
| L7 Applications | Complete | 2 |
| L8 Advanced Topics | Partial | 1 |
| L9 Research Frontiers | Partial | 1 |
| **Total** | | **16/18** |

## Core Definitions

### Coordinate Systems
| System | Symbol | Description |
|--------|--------|-------------|
| ECEF | (x, y, z) | Earth-Centered Earth-Fixed Cartesian [m] |
| LLA | (φ, λ, h) | Geodetic Latitude, Longitude, Altitude |
| ENU | (e, n, u) | East-North-Up local tangent plane |

### GPS Signal
| Parameter | Value |
|-----------|-------|
| C/A code length | 1023 chips |
| C/A code period | 1 ms |
| Chip rate | 1.023 Mcps |
| L1 frequency | 1575.42 MHz |
| L2 frequency | 1227.60 MHz |
| L5 frequency | 1176.45 MHz |

### DOP Types
| DOP | Formula | Meaning |
|-----|---------|---------|
| GDOP | √(q₁₁+q₂₂+q₃₃+q₄₄) | Total geometric effect |
| PDOP | √(q₁₁+q₂₂+q₃₃) | 3D position effect |
| HDOP | √(c₁₁+c₂₂) in ENU | Horizontal effect |
| VDOP | √(c₃₃) in ENU | Vertical effect |
| TDOP | √(q₄₄) | Clock bias effect |

## Core Theorems

### Pseudorange Equation
```
P = ρ + c·(dt_r - dt_s) + I + T + M + ε
```
where ρ = geometric range, c = speed of light, dt_r/dt_s = clock biases,
I = ionospheric delay, T = tropospheric delay, M = multipath, ε = noise.

### DOP-Accuracy Relationship
```
σ_position = PDOP × σ_UERE
σ_horizontal = HDOP × σ_UERE
σ_vertical = VDOP × σ_UERE
```
PDOP × UERE translates ranging error to position error.

### Kepler's Equation
```
M = E - e·sin(E)
```
Solved via Newton-Raphson: E_{k+1} = E_k - (E_k - e·sin(E_k) - M)/(1 - e·cos(E_k))

### Klobuchar Ionospheric Model
```
I_z = A₁ + A₂·cos(2π·(t - A₃)/A₄)   [seconds]
I = c · F · I_z · (f₁²/f²)            [meters]
```

### Hatch Carrier Smoothing
```
P_sm[n] = (1/N)·P[n] + (1 - 1/N)·(P_sm[n-1] + Φ[n] - Φ[n-1])
```
Noise reduction: σ_smooth = σ_raw / √N

## Core Algorithms

| Algorithm | Function | Complexity | Reference |
|-----------|----------|------------|-----------|
| Bancroft direct solution | `gnss_bancroft_solve()` | O(n_sats) | Bancroft (1985) |
| Iterative LS PVT | `gnss_ls_position_solve()` | O(n·iter) | IS-GPS-200 |
| Weighted LS | `gnss_wls_position_solve()` | O(n·iter) | Strang & Borre |
| Klobuchar ionosphere | `gnss_iono_klobuchar()` | O(1) | Klobuchar (1987) |
| Saastamoinen troposphere | `gnss_tropo_saastamoinen()` | O(1) | Saastamoinen (1972) |
| Hatch carrier smoothing | `gnss_hatch_smooth()` | O(1)/epoch | Hatch (1982) |
| Kepler solver (Newton) | `gnss_kepler_solve()` | O(iter) | — |
| C/A code generation | `gnss_ca_code_generate()` | O(1023) | IS-GPS-200 |
| Satellite position | `gnss_satpos_from_ephemeris()` | O(20) | IS-GPS-200 |
| RAIM FDE | `gnss_raim_fde()` | O(n²) | Brown (1992) |
| Best subset selection | `gnss_select_best_subset()` | O(n²) | — |
| Cycle slip detection | `gnss_slip_detect()` | O(1)/sat | Blewitt (1990) |

## Classic Problems

1. **Single-Point Positioning** — Solve for (x,y,z,Δt) from ≥4 pseudorange measurements
2. **DOP Analysis** — Predict positioning accuracy from satellite geometry
3. **Carrier-Smoothed Pseudorange** — Reduce code noise using carrier phase
4. **Bancroft Initialization** — Direct algebraic solution for LS seed
5. **Atmospheric Correction** — Apply Klobuchar + Saastamoinen models
6. **Cycle Slip Detection** — GF + Doppler phase-rate comparison
7. **RAIM Integrity** — Fault detection from residual analysis

## Course Mapping

| School | Courses | Topics |
|--------|---------|--------|
| MIT | 6.003, 16.36 | GPS signal structure, correlation |
| Stanford | AA272C, EE359 | GPS position solution, DOP, Bancroft |
| Berkeley | EE123, EE221A | Optimal estimation, LS positioning |
| Illinois | ECE 310, ECE 459 | Spread spectrum, GPS receivers |
| Michigan | EECS 351, EECS 455 | GPS signal processing |
| Georgia Tech | ECE 4270, ECE 6601 | GNSS receiver design |
| TU Munich | Navigation, Satellite Comm | GNSS architecture |
| ETH Zurich | 227-0427, 227-0436 | Precise positioning, ambiguity |
| Tsinghua | 卫星导航原理, 通信原理 | BD/GPS positioning algorithms |

## Build & Test

```bash
make          # Build library and all targets
make test     # Build and run 22-test suite
make examples # Build 3 end-to-end examples
make bench    # Build and run benchmarks
make clean    # Remove build artifacts
```

## File Structure

```
mini-gnss-gps/
├── Makefile
├── README.md                          ← This file (COMPLETE ✅)
├── include/
│   ├── gnss_common.h                  Core types, WGS84, coordinates, matrices
│   ├── gnss_signal.h                  C/A code, Doppler, ephemeris, Kepler
│   ├── gnss_pseudorange.h             Pseudorange model, Klobuchar, Saastamoinen
│   ├── gnss_position.h                Bancroft, LS, WLS, RAIM position solvers
│   ├── gnss_dop.h                     DOP theory, geometry analysis
│   └── gnss_carrier.h                 Carrier phase, smoothing, DGPS, combinations
├── src/
│   ├── gnss_common.c                  ECEF↔LLA↔ENU, 3×3/4×4 matrix algebra, time
│   ├── gnss_signal.c                  C/A Gold codes, satellite position, Doppler
│   ├── gnss_pseudorange.c             Klobuchar, Saastamoinen, Hopfield, CMC, UERE
│   ├── gnss_position.c                Bancroft, iterative LS, weighted LS, RAIM
│   ├── gnss_dop.c                     DOP computation, constellation analysis
│   ├── gnss_carrier.c                 Hatch filter, cycle slip, IF/MW/GF combos
│   └── gnss_formal.lean               Lean 4: GPS time, WGS84, DOP, pseudorange
├── tests/
│   └── test_gnss.c                    22-test suite with mathematical assertions
├── examples/
│   ├── example_position_solve.c       SPP at SFO with 8 satellites
│   ├── example_dop_calc.c             DOP analysis & subset selection
│   └── example_carrier_smoothing.c    Carrier smoothing with cycle slips
├── docs/
│   ├── knowledge-graph.md             L1-L9 complete knowledge map
│   ├── coverage-report.md             Coverage assessment (16/18)
│   ├── gap-report.md                  Missing L8/L9 items + priority
│   ├── course-alignment.md            Nine-school curriculum mapping
│   └── course-tree.md                 Prerequisite dependency tree
└── benches/
    └── bench_gnss.c                   Micro benchmarks
```

## References

- Tsui, J.B.Y. (2005). *Fundamentals of Global Positioning System Receivers*, 2nd ed. Wiley.
- Parkinson, B.W. & Spilker, J.J. (1996). *GPS: Theory and Applications*. AIAA.
- Misra, P. & Enge, P. (2011). *Global Positioning System*, 2nd ed. Ganga-Jamuna Press.
- Hofmann-Wellenhof, B. et al. (2007). *GNSS — GPS, GLONASS, Galileo*, 2nd ed. Springer.
- Teunissen, P.J.G. & Montenbruck, O. (2017). *Springer Handbook of GNSS*. Springer.
- Strang, G. & Borre, K. (1997). *Linear Algebra, Geodesy, and GPS*. Wellesley-Cambridge.
- Bancroft, S. (1985). "An Algebraic Solution of the GPS Equations." IEEE Trans. AES, 21(1).
- Hatch, R. (1982). "The Synergism of GPS Code and Carrier Measurements." Proc. 3rd Int'l Geodetic Symp.
- Klobuchar, J.A. (1987). "Ionospheric time-delay algorithm for single-frequency GPS users." IEEE Trans. AES.
- Saastamoinen, J. (1972). "Contributions to the theory of atmospheric refraction." Bulletin Géodésique.
- IS-GPS-200, *NAVSTAR GPS Space Segment / Navigation User Interface Specification*.
