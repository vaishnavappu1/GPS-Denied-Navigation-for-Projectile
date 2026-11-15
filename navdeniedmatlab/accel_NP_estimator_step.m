function [x_np_rad, P_np_rad, x_np_ax, P_np_ax, accel_clean] = ...
    accel_NP_estimator_step(x_np_rad, P_np_rad, x_np_ax, P_np_ax, ...
    accel_meas, f_aero_body, pb_est, Q_np_rad, R_np_rad, Q_np_ax, R_np_ax)
% ACCELEROMETER NP ESTIMATOR — physical 4-state radial model
% Patent US 9,702,674 B2, FIG.9
%
% IMPROVEMENT OVER PREVIOUS VERSION:
% Previous state: x = [NP_j; NP_k]  (total contamination per channel)
%   Problem: NP = b + r*pb^2 changes as pb decays, forcing the filter
%   to chase a moving target throughout the flight.
%
% This version: x = [b_j; b_k; r_j; r_k]
%   b_j, b_k = sensor biases [m/s^2]  -- CONSTANT (fixed sensor property)
%   r_j, r_k = CG offsets [m]         -- CONSTANT (fixed geometry)
%   Centripetal term computed analytically: a_c = r * pb^2
%   using pb_est from Stage E (already accurate to ~0.7 Hz RMS)
%
% Why this is better:
%   The filter now estimates two fixed physical constants instead of
%   a changing acceleration. As pb decays 200->150 Hz over flight,
%   the old filter had to continuously readjust. This filter's states
%   stay constant -- only the H matrix changes each step (via pb^2).
%   This is more accurate, especially in the second half of flight.
%
% MEASUREMENT MODEL (linear in states):
%   z_j = accel_j - f_aero_j = b_j + r_j * pb^2
%   z_k = accel_k - f_aero_k = b_k + r_k * pb^2
%
%   H = [1  0  pb^2   0  ]
%       [0  1   0    pb^2]
%
% NOTE ON CONDITIONING: H contains pb^2 ~ 1.58e6 at 200 Hz. This creates
% large scale differences between bias states (~1) and offset states (~1e6).
% Joseph-form covariance update is MANDATORY here for exactly the same
% reason it was necessary in Stage E. P is also symmetrized each step.
%
% AXIAL: unchanged from previous version (already at noise floor 0.05 m/s^2)
%
% INPUTS:
%   x_np_rad = [b_j; b_k; r_j; r_k] radial NP state (4x1)
%              b in [m/s^2], r in [m]
%   P_np_rad = 4x4 radial NP covariance
%   x_np_ax  = scalar axial NP state [m/s^2]
%   P_np_ax  = scalar axial NP variance
%   accel_meas   = 3x1 raw body accel [m/s^2], order [X; Y; Z]
%   f_aero_body  = 3x1 predicted aero specific force [m/s^2], [X; Y; Z]
%   pb_est       = scalar spin rate estimate from Stage E [rad/s]
%   Q_np_rad = 4x4 process noise (tiny -- states are physical constants)
%   R_np_rad = 2x2 measurement noise (large during balloting, small after)
%   Q_np_ax  = scalar axial process noise
%   R_np_ax  = scalar axial measurement noise
%
% OUTPUTS:
%   x_np_rad, P_np_rad = updated radial state and covariance
%   x_np_ax,  P_np_ax  = updated axial state and variance
%   accel_clean = 3x1 NP-corrected accel [m/s^2]

% ================================================
% RADIAL NP ESTIMATOR
% ================================================

% Predict (identity process -- both b and r are constants)
x_rad_pred = x_np_rad;
P_rad_pred = P_np_rad + Q_np_rad;

% H matrix: depends on current pb^2
pb2 = pb_est^2;
H_rad = [1, 0, pb2, 0;
         0, 1, 0,   pb2];

% Measurement
z_rad = [accel_meas(2) - f_aero_body(2);
         accel_meas(3) - f_aero_body(3)];
y_rad = z_rad - H_rad * x_rad_pred;

% Kalman gain
S_rad = H_rad * P_rad_pred * H_rad' + R_np_rad;
K_rad = P_rad_pred * H_rad' / S_rad;

% Update
x_np_rad = x_rad_pred + K_rad * y_rad;

% Joseph form -- mandatory due to scale mismatch between b and r states
I4 = eye(4);
IKH_rad  = I4 - K_rad * H_rad;
P_np_rad = IKH_rad * P_rad_pred * IKH_rad' + K_rad * R_np_rad * K_rad';
P_np_rad = (P_np_rad + P_np_rad') / 2;   % symmetry guard

% ================================================
% AXIAL NP ESTIMATOR (unchanged -- already at noise floor)
% ================================================

x_ax_pred = x_np_ax;
P_ax_pred = P_np_ax + Q_np_ax;

z_ax = accel_meas(1) - f_aero_body(1);
y_ax = z_ax - x_ax_pred;   % H = 1

S_ax = P_ax_pred + R_np_ax;
K_ax = P_ax_pred / S_ax;

x_np_ax = x_ax_pred + K_ax * y_ax;

IKH_ax  = 1 - K_ax;
P_np_ax = IKH_ax * P_ax_pred * IKH_ax + K_ax * R_np_ax * K_ax;

% ================================================
% COMPUTE NP AND CLEAN ACCEL OUTPUT
% ================================================
% Reconstruct total NP from physical states + current pb
NP_j = x_np_rad(1) + x_np_rad(3) * pb2;   % b_j + r_j * pb^2
NP_k = x_np_rad(2) + x_np_rad(4) * pb2;   % b_k + r_k * pb^2
NP_i = x_np_ax;

accel_clean = accel_meas - [NP_i; NP_j; NP_k];

end