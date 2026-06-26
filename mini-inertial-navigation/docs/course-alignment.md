# Course Alignment - mini-inertial-navigation

## MIT
- 6.832 Underactuated Robotics: Quaternion kinematics, 3D rotations
  -> ins_attitude.c (quaternion update, rotation)
- 6.003 Signals and Systems: Coordinate transforms, vector ops
  -> ins_core.c (vector/matrix ops)

## Stanford
- AA272 Global Positioning Systems: INS mechanization, GPS/INS KF
  -> ins_mechanization.c, ins_integration.c
- EE267 Virtual Reality: IMU orientation tracking
  -> example_ahrs.c (complementary filter)

## Berkeley
- EE C128 Mechatronics: IMU sensor characterization
  -> ins_calibration.c (6-position, multi-position calib)
- EE123 DSP: Allan variance, sensor noise analysis
  -> ins_errors.c (Allan variance computation)

## Michigan
- EECS 455 Communications: Kalman filter estimation
  -> ins_integration.c (15-state error-state KF)
- EECS 461 Embedded Control: Sensor fusion
  -> ins_integration.c (ZUPT, GPS/INS loop)

## Georgia Tech
- ECE 6601 Communications: Navigation systems
  -> ins_integration.c (GNSS outage budget)

## TU Munich
- Navigation (Inertial Systems): Strapdown INS lab
  -> ins_mechanization.c, example_strapdown.c

## ETH Zurich
- 151-0854 Trajectory Generation: Quaternion splines
  -> ins_attitude.c (quaternion algebra)

## Tsinghua University
- Inertial Navigation: Full strapdown simulation
  -> ins_mechanization.c, ins_advanced.c

## Caltech
- CDS 110 Optimal Control: State estimation
  -> example_ahrs.c, ins_advanced.c (particle filter)
