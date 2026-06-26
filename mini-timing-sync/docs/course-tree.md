# Prerequisite Dependency Tree — mini-timing-sync

## Upstream Dependencies (What this module needs)
```
Digital Signal Processing (Fourier, filtering, estimation)
  -> mini-digital-signal-process (DSP foundations)
Communication Principles (modulation, demodulation, timing recovery)
  -> mini-communication-principle
Navigation/Positioning context (GPS, GNSS)
  -> mini-navigation-positioning (parent module)
```

## Internal Dependency Tree
```
L1: Definitions
  -> L2: Core Concepts (PI servo, two-way, holdover, BMCA)
  --> L3: Mathematical Structures (Kalman, regression, state-space)
  ----> L4: Fundamental Laws (offset/delay equations, Allan)
  -------> L5: Algorithms (Allan computation, BMCA, Marzullo, noise fit)
  ----------> L6: Canonical Problems (PTP sync, oscillator char, NTP)
  -------------> L7: Applications (5G, power grid, GPSDO, MiFID)
  ----------------> L8: Advanced (3-cornered hat, White Rabbit, TWSTFT)
  -------------------> L9: Frontiers (quantum, optical, CSAC)
```

## Downstream Dependents (What needs this module)
```
Wireless/Mobile Communications (5G fronthaul timing)
  -> mini-wireless-mobile-comm
Radar remote sensing (pulse timing, coherent integration)
  -> mini-radar-remote-sensing
Optical fiber communications (clock distribution)
  -> mini-optical-fiber-comm
Industrial fieldbus (TSN, EtherCAT distributed clocks)
  -> mini-industrial-fieldbus
Audio/Video engineering (genlock, word clock)
  -> mini-audio-video-eng
IoT edge computing (time-sensitive networking)
  -> mini-iot-edge-computing
```
