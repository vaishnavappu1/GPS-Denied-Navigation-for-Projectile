function [pos_new, vel_new, Phi_a, P_a] = track2_inertial_step(...
    pos_prev, vel_prev, P_prev_combined, ...
    accel_clean_prev, accel_clean_curr, ...
    phi_prev, theta_prev, psi_prev, ...
    phi_curr, theta_curr, psi_curr, ...
    dt, p, Q_a)
% TRACK 2: Accelerometer/Inertial Position/Velocity Predictor
% Patent US 9,702,674 B2, FIG.6 block 650, FIG.10
%
% PURPOSE:
% Propagates position and velocity forward one timestep by rotating
% NP-corrected body-frame accelerometer measurements to inertial frame
% using Stage E's attitude estimate, adding gravity, and trapezoidal
% integrating (patent eq.38, FIG.10 blocks 1010/1020/1030).
%
% STATE TRANSITION Phi_a (EXACT, not numeric):
% Unlike Track 1, Track 2's "acceleration" is a MEASURED quantity
% (the NP-corrected accelerometer reading), not a function of the
% [pos; vel] state being propagated. There is no force-Jacobian term
% because the integration is purely kinematic given the measurement:
%   pos_new = pos_prev + dt*vel_prev  (to first order)
%   vel_new = vel_prev + dt*accel_measured (accel is external input)
% So:
%   Phi_a = [I,    dt*I]
%           [0,    I   ]
% This is the same kinematic block used in Track 1's upper-left/
% upper-right (also exact there), but Track 2's lower-right is
% exactly I (no aero Jacobian term) since acceleration doesn't depend
% on the propagated velocity state at all -- it comes from the sensor.
%
% PROCESS NOISE Q_a:
% Reflects uncertainty in the NP-corrected accelerometer signal
% (residual NP estimation error + sensor noise), propagated through
% the integration. This is the source of Track 2's unbounded drift --
% there is no physical model to bound it, unlike Track 1's aero model.
% Should be set larger than Track 1's Q_d to reflect this.
%
% INPUTS:
%   pos_prev, vel_prev   = 3x1 previous position [m], velocity [m/s]
%   P_prev_combined      = 6x6 combined covariance from previous step
%                          (shared prior with Track 1)
%   accel_clean_prev     = 3x1 NP-corrected body accel at previous step [m/s^2]
%   accel_clean_curr     = 3x1 NP-corrected body accel at current step [m/s^2]
%   phi_prev,theta_prev,psi_prev = Euler angles at previous step [rad]
%                                  (phi from Stage D, theta/psi from Stage E)
%   phi_curr,theta_curr,psi_curr = Euler angles at current step [rad]
%   dt   = timestep [s]
%   p    = projectile_params struct
%   Q_a  = 6x6 process noise covariance for inertial track
%
% OUTPUTS:
%   pos_new = 3x1 predicted position [m]
%   vel_new = 3x1 predicted velocity [m/s]
%   Phi_a   = 6x6 state transition matrix (exact, kinematic)
%   P_a     = 6x6 predicted covariance for Track 2

% ================================================
% ROTATE BODY ACCEL TO INERTIAL FRAME, ADD GRAVITY
% (patent FIG.10 blocks 1010, 1020)
% ================================================
g_inertial = [0; 0; p.g];

R_b2i_prev = build_R_b2i(phi_prev, theta_prev, psi_prev);
R_b2i_curr = build_R_b2i(phi_curr, theta_curr, psi_curr);

a_prev = R_b2i_prev * accel_clean_prev + g_inertial;
a_curr = R_b2i_curr * accel_clean_curr + g_inertial;

% ================================================
% TRAPEZOIDAL INTEGRATION (patent eq.38, FIG.10 block 1030)
% ================================================
vel_new = vel_prev + 0.5*dt*(a_prev + a_curr);
pos_new = pos_prev + 0.5*dt*(vel_prev + vel_new);

% ================================================
% EXACT (NOT NUMERIC) PHI_A
% Acceleration is a measured input, independent of [pos;vel] state,
% so the lower-right block is exactly I, no Jacobian needed.
% ================================================
Phi_a = [eye(3),    dt*eye(3);
    zeros(3,3), eye(3)  ];

% ================================================
% COVARIANCE PROPAGATION (patent eq.40, FIG.10 block 1040)
% ================================================
P_a = Phi_a * P_prev_combined * Phi_a' + Q_a;
P_a = (P_a + P_a') / 2;  % symmetry guard

end


function R_b2i = build_R_b2i(phi, theta, psi)
% Body to inertial rotation, full Euler angle composition.
Ry = [cos(theta), 0, -sin(theta); 0, 1, 0; sin(theta), 0, cos(theta)];
Rz = [cos(psi), sin(psi), 0; -sin(psi), cos(psi), 0; 0, 0, 1];
Rx = [1, 0, 0; 0, cos(phi), -sin(phi); 0, sin(phi), cos(phi)];
R_b2i = (Rx * Ry * Rz)';  % R_i2b transpose = body to inertial
end