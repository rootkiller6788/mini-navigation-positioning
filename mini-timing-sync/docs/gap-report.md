# Gap Report — mini-timing-sync

## Missing Knowledge Points

### L8 Advanced Topics (not yet implemented)
| # | Topic | Priority | Effort |
|---|-------|----------|--------|
| 1 | Kalman-filtered multi-sensor clock ensemble | Medium | High |
| 2 | AI/ML-based network delay prediction | Low | High |
| 3 | Monte Carlo uncertainty propagation for time transfer | Medium | Medium |

### L9 Research Frontiers (documented, not implemented)
| # | Topic | Reference | Feasibility |
|---|-------|-----------|-------------|
| 1 | Quantum clock synchronization (entanglement-based) | Komar et al., Nature Physics 2014 | 5-10 years |
| 2 | Optical clock networks (Sr/Yb lattice clocks) | Boulder Atomic Clock, NIST 2020 | 3-5 years |
| 3 | Chip-scale atomic clocks (CSAC) for portable timing | Microsemi SA.45s, NIST on-chip | Available now |
| 4 | 6G integrated sensing and timing (ISAC) | 3GPP Release 20+ | 5+ years |
| 5 | Free-space optical time transfer (satellite-to-ground) | T2T2 mission concept | 5-10 years |

## Improvement Recommendations
1. Add multi-master clock ensemble with Kalman fusion (L8)
2. Add Monte Carlo time transfer uncertainty propagation (L8)
3. Extend White Rabbit to full WR-PTP profile (L8)
4. Add CSAC drift model for long-holdover applications (L7/L8)
