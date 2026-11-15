# GPS-Denied Navigation for Spin-Stabilized Projectiles

MATLAB simulation of a multi-sensor Kalman-filter navigation system for a
155 mm spin-stabilized artillery projectile operating without GPS. Developed
during an internship at SpaceLabs Analytics and Dynamics, Trivandrum, based on
the architecture in US Patent 9,702,674 B2.

## Overview
Spin-stabilized smart munitions rely on GPS for precision guidance, but GPS can
be jammed, spoofed, or blocked. This project builds a contingency navigation
capability that uses only onboard inertial and magnetic sensors, and evaluates
it in a full closed-loop MATLAB simulation.

## What it does
- 6-DOF truth simulation (aerodynamics, Magnus effect, coning)
- Realistic sensor models: triaxial magnetometer, accelerometer, rate gyro, intermittent GPS
- Magnetometer de-biasing (Stage C) and roll/spin extraction (Stage D)
- 8-state gyro-augmented Extended Kalman Filter for attitude and spin rate (Stage E)
- Accelerometer normal-force estimator (Layer 6A)
- Two independent position/velocity tracks (flight-dynamic and inertial)
- Bar-Shalom track-to-track fusion with GPS-reset handling (Layer 6C)
- Open-loop and closed-loop evaluation with a 15-trial Monte Carlo analysis

## Key findings
- A single 3-axis magnetometer gives a rank-2 observation matrix (confirmed by SVD),
  so full 3D attitude cannot be recovered from it alone.
- A 19:1 pitch/yaw Jacobian sensitivity ratio makes yaw weakly observable,
  causing systematic East-axis position drift during long GPS outages.
- Open-loop pitch RMS error: 0.227°. Closed-loop: 5.098°.

## Simulation setup
- Projectile: 155 mm M795 (43.2 kg), muzzle velocity 827 m/s, 45° elevation
- Initial spin: 50 Hz, decaying to about 38 Hz over the flight
- Trajectory: about 8.8 km apogee, about 26 km range
- Environment: ISA atmosphere model, Earth magnetic field in NED frame

## Repository contents
This repository contains the MATLAB scripts developed for the project,
along with the full internship report.

## Report
The full derivations, simulation results, and discussion are in the
internship report (see the `docs/` folder).

## Reference
S. H. Hall et al., "Navigation system and method for high-spin projectiles,"
US Patent 9,702,674 B2, 2017.
