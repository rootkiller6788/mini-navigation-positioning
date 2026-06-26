# Course Dependency Tree: UWB Localization
## Prerequisites
- Linear Algebra -> Optimization -> Multilateration -> Positioning
- Statistics -> CRLB -> Error Budgets -> Quality Assessment
- Signal Processing -> CIR Analysis -> NLOS Detection -> Mitigation
- Electromagnetics -> Channel Models -> Path Loss -> Link Budget
- Kalman Theory -> EKF -> Tracking -> RTS Smoother
## Sub-Module Dependencies
- uwb_types.h -> (all)
- uwb_mathematics.h -> uwb_positioning.c
- uwb_ranging.h -> uwb_positioning.c, uwb_tracking.c
- uwb_channel.h -> uwb_nlos.c
- uwb_nlos.h -> uwb_ranging.c
- uwb_positioning.h -> uwb_tracking.c
