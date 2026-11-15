function [NP_y, NP_z, x_np, P_np] = radial_NP_estimator_step(x_np, P_np, ...
    accel_y_meas, accel_z_meas, ay_pred, az_pred, Q_np, R_np)
% RADIAL ACCELEROMETER NP ESTIMATOR
% Patent US 9,702,674 B2, FIG.9, block 910.
%
% WHAT THIS DOES: estimates and removes the nuisance parameter (NP)
% contamination from the two radial (body Y,Z) accelerometer channels.
% The contamination = constant bias + cpb2 term (centripetal offset).
% In the body frame, this contamination appears as a slowly-varying
% offset on each channel independently.
%
% WHY THIS IS SIMPLE: the patent notes "from an integration standpoint,
% these NP errors will produce no measurable effect on the velocity and
% position estimates" -- meaning the NP doesn't need to be tracked with
% a complex sinusoid model for navigation purposes. What we need is just
% the CLEAN aerodynamic signal. Since the predicted aerodynamic signal
% (ay_pred, az_pred) comes from eom.m's force equations (already
% validated), the difference between measured and predicted is a direct
% observation of NP -- forming a simple, linear, scalar-per-channel
% measurement equation.
%
% STATE: x_np = [NP_y; NP_z] (2x1, one NP per radial channel)
% PROCESS: NP_dot = 0 (slowly varying, small Q)
% MEASUREMENT: z = accel_meas - accel_pred = NP (direct observation)
%
% INPUTS:
% x_np, P_np = previous state (2x1) and covariance (2x2)
% accel_y_meas = raw Y accelerometer reading [m/s^2]
% accel_z_meas = raw Z accelerometer reading [m/s^2]
% ay_pred = predicted true Y aerodynamic acceleration [m/s^2]
% from eom.m's force equations in body frame
% az_pred = predicted true Z aerodynamic acceleration [m/s^2]
% Q_np = 2x2 process noise covariance
% R_np = 2x2 measurement noise covariance
%
% OUTPUTS:
% NP_y, NP_z = estimated NP on each radial channel [m/s^2]
% x_np, P_np = updated state and covariance

% --- PREDICT ---
x_pred = x_np; % identity process model
P_pred = P_np + Q_np;

% --- UPDATE (both channels simultaneously, H = eye(2)) ---
H = eye(2);
z = [accel_y_meas - ay_pred;
    accel_z_meas - az_pred];
y = z - H*x_pred;

S = H*P_pred*H' + R_np;
K = P_pred*H'/S;

x_np = x_pred + K*y;

% Joseph form for numerical stability (proven necessary in Stage E)
IKH = eye(2) - K*H;
P_np = IKH*P_pred*IKH' + K*R_np*K';

NP_y = x_np(1);
NP_z = x_np(2);
end
