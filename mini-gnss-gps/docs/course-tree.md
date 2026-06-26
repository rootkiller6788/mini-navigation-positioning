# Course Tree — mini-gnss-gps

## Prerequisite Dependencies

```
mini-gnss-gps
├── mini-signal-basis (L2: correlation, convolution for GPS acquisition)
├── mini-fourier-analysis (L3: FFT for parallel code-phase search)
├── mini-filter-theory (L4: Hatch smoothing as low-pass filter)
├── mini-system-analysis (L4: least squares estimation theory)
├── mini-coordinate-systems (L3: ECEF/LLA/ENU transforms)
├── mini-communication-theory (L4: spread spectrum, Gold codes)
├── mini-noise-theory (L4: SNR, UERE error budget)
├── mini-wireless-propagation (L7: ionospheric/tropospheric models)
├── mini-adaptive-filter (L8: Kalman filter for PVT)
└── mini-inertial-navigation (L8: GNSS/INS coupling)
```

## Downstream Dependencies

```
mini-gnss-gps
├── mini-integrated-navigation (sensor fusion with GPS)
├── mini-slam-basics (GPS for loop closure)
├── mini-autonomous-vehicle (localization)
├── mini-timing-sync (GPS-disciplined oscillators)
└── mini-geomagnetic-navigation (backup to GNSS)
```

## Knowledge Dependency Graph (L1-L9)

```
L1: WGS84, time, signals  ───┐
L2: Trilateration, codes ────┤
L3: Linear algebra, rotation─┤
L4: DOP, Kepler, pseudorange─┼── converge to L6 canonical SPP
L5: LS, Bancroft, atmosphere─┘
                                    │
L6: SPP ──── L7: Aviation, urban canyon
       │
       └────── L8: RTK, RAIM, multi-GNSS
                    │
                    └─────── L9: LEO-PNT, 6G positioning
```

## Learning Sequence (recommended)

1. L1: Understand GPS signal structure, WGS84 datum, time systems
2. L2: Study trilateration geometry, C/A code, satellite orbits
3. L3: Master linear algebra for positioning: design matrix, normal equations
4. L4: Apply Kepler's laws, pseudorange equation, DOP theory
5. L5: Implement Bancroft, LS, atmospheric models, smoothing
6. L6: Solve full SPP problem end-to-end
7. L7: Apply to real scenarios (aviation, urban, surveying)
8. L8: Explore advanced topics (RTK, RAIM, EKF)
9. L9: Research directions (LEO-PNT, quantum, 6G)
