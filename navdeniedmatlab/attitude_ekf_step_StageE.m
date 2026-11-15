function [x_new, P_new, gate_tripped] = attitude_ekf_step_StageE(x, P, mag_raw, phi_roll, vel_est, p, qbar, V, dt, Q, R_axial, pb_meas, R_pb, gyro_meas)
% STAGE E v8: main attitude Kalman filter, EXPANDED to 8 states
% (US 9,702,674 B2, FIG.7, block 730) -- adds direct gyro rate
% measurement update (q, r) with spin-rate-adaptive R, per Phase 1/2.
%
% *** THIS IS AN EXTENSION OF THE VALIDATED 6-STATE VERSION ***
% All 6-state logic (predict, axial mag update, pb update, clamps,
% Joseph-form covariance, NaN guard, innovation gate) is UNCHANGED and
% copied verbatim. Only additions:
% 1. State vector grows from 6 to 8: appends [b_q; b_r] (gyro rate
% biases), per Phase 1's gyro_h_StageE design.
% 2. Q, P must be 8x8 by the caller (not validated/reshaped here --
% caller's responsibility, flagged explicitly below).
% 3. A THIRD measurement update is added: the gyro rate channel,
% using gyro_h_StageE (Phase 1) and gyro_R_scaling (Phase 2).
%
% CALL-SITE RESPONSIBILITY (not enforced inside this function):
% - x must be 8x1: [dtheta_S; dpsi_S; q; r; pb; axial_bias; b_q; b_r]
% - P, Q must be 8x8. The caller must have grown these from the old
% 6x6 versions -- e.g. by embedding the old 6x6 P/Q in the top-left
% block and adding diagonal entries for b_q, b_r (suggested initial
% variance and process noise below, but these are ASSUMPTIONS same
% as Phase 1/2 -- no datasheet backing).
% - gyro_meas is 2x1: [q_meas; r_meas], the raw (biased) gyro rate
% reading, in the SAME fixed-plane sense as q, r (see Phase 1's
% documented simplification: no body-to-fixed-plane rotation is
% applied here).
%
% STATE: x = [dtheta_S; dpsi_S; q; r; pb; axial_bias; b_q; b_r] (8 states)
%
% INPUTS (all as 6-state version, PLUS):
% gyro_meas = 2x1 raw gyro rate measurement [q_meas; r_meas] [rad/s]
%
% OUTPUTS: unchanged (x_new, P_new now 8x1/8x8; gate_tripped as before,
% refers only to the axial-mag gate, same as the 6-state version).

% ================================================
% PREDICT (states 1-6 identical dynamics to 6-state version;
% states 7-8, b_q/b_r, modeled as constant -- random walk via Q only,
% same convention as axial_bias_dot = 0 above)
% ================================================
dtheta_S = x(1); dpsi_S = x(2); q = x(3); r = x(4); pb = x(5); axial_bias = x(6);
b_q = x(7); b_r = x(8); %#ok<NASGU> -- unused in dynamics, present for clarity

L_B = qbar * p.A * p.d * (p.Clp * pb * p.d / (2*V));
M_B = qbar * p.A * p.d * (p.Cma*dtheta_S + p.Cmq*q*p.d/(2*V) + p.Cnpa*pb*p.d/(2*V)*dpsi_S);
N_B = qbar * p.A * p.d * (-p.Cma*dpsi_S + p.Cmq*r*p.d/(2*V) + p.Cnpa*pb*p.d/(2*V)*dtheta_S);

dtheta_S_dot = q;
dpsi_S_dot = r;
pb_dot = L_B / p.Ix;
q_dot = pb*r*(p.Ix/p.Iy) + M_B/p.Iy;
r_dot = -pb*q*(p.Ix/p.Iy) + N_B/p.Iy;
axial_bias_dot = 0;
b_q_dot = 0;
b_r_dot = 0;

x_pred = x + dt * [dtheta_S_dot; dpsi_S_dot; q_dot; r_dot; pb_dot; axial_bias_dot; b_q_dot; b_r_dot];
F = numeric_F_StageE8(x, p, qbar, V, dt);
P_pred = F*P*F' + Q;

% ================================================
% UPDATE 1 -- axial magnetometer channel (unchanged logic, H extended
% to 8 columns; new b_q/b_r columns are exactly zero since the axial
% mag model has no dependence on gyro states)
% ================================================
psi_v = atan2(vel_est(2), vel_est(1));
theta_v = atan2(-vel_est(3), sqrt(vel_est(1)^2+vel_est(2)^2));

if any(isnan(vel_est)) || isnan(phi_roll) || any(isnan(x_pred)) || isnan(qbar) || isnan(V) || any(isnan(gyro_meas))
warning('attitude_ekf_step_StageE_v8: NaN detected. vel_est=[%g %g %g] phi_roll=%g qbar=%g V=%g x_pred=[%s]', ...
vel_est(1), vel_est(2), vel_est(3), phi_roll, qbar, V, mat2str(x_pred'));
x_new = x;
P_new = P;
gate_tripped = true;
return;
end

h0 = mag_axial_h_StageE8(x_pred, phi_roll, theta_v, psi_v, p.B_earth);
H = mag_axial_jacobian_numeric8(x_pred, phi_roll, theta_v, psi_v, p.B_earth);
z = mag_raw(1);
y = z - h0;

innovation_gate_sigma = 6;
S_pred = H*P_pred*H' + R_axial;
axial_gate_tripped = abs(y) > innovation_gate_sigma * sqrt(S_pred);

if axial_gate_tripped
x_new = x_pred;
P_new = P_pred;
else
S = S_pred;
K = P_pred*H' / S;
x_new = x_pred + K*y;
IKH = eye(8) - K*H;
P_new = IKH*P_pred*IKH' + K*R_axial*K';
end
gate_tripped = axial_gate_tripped;

% ================================================
% UPDATE 2 -- pb measurement (Stage 2 zero-crossing), unchanged logic,
% H_pb extended to 8 columns (still only touches index 5)
% ================================================
H_pb = [0 0 0 0 1 0 0 0];
y_pb = pb_meas - x_new(5);
S_pb = H_pb*P_new*H_pb' + R_pb;
K_pb = P_new*H_pb' / S_pb;
x_new = x_new + K_pb*y_pb;
IKH_pb = eye(8) - K_pb*H_pb;
P_new = IKH_pb*P_new*IKH_pb' + K_pb*R_pb*K_pb';

% ================================================
% UPDATE 3 (NEW) -- gyro rate channel [q; r], adaptive R via spin rate
% Uses Phase 1's gyro_h_StageE / gyro measurement model, and Phase 2's
% gyro_R_scaling, evaluated at the PREVIOUS pb estimate (x(5), i.e.
% pb BEFORE this function's own pb update above) to avoid same-step
% self-reference -- consistent with the causal-delay convention used
% for vel_c_log(i-1,:) elsewhere in this project.
% ================================================
pb_for_R = x(5); % pre-update pb, i.e. the value entering this function
R_gyro_i = gyro_R_scaling(pb_for_R, p);
R_gyro_diag = diag([R_gyro_i, R_gyro_i]);

% --- CORRECTED (Phase 1b): gyro reads BODY frame; state's q,r are
% S-frame (fixed-plane). Rotate via phi_roll (external, Stage D) before
% comparing -- same Rx(phi_roll) convention used by the axial mag model.
H_gyro = gyro_jacobian_numeric8(x_new, phi_roll);
h_gyro = gyro_h_StageE8(x_new, phi_roll);
y_gyro = gyro_meas - h_gyro;
S_gyro = H_gyro*P_new*H_gyro' + R_gyro_diag;
K_gyro = P_new*H_gyro' / S_gyro;
x_new = x_new + K_gyro*y_gyro;
IKH_gyro = eye(8) - K_gyro*H_gyro;
P_new = IKH_gyro*P_new*IKH_gyro' + K_gyro*R_gyro_diag*K_gyro';

% --- Hard physical-plausibility clamp (unchanged 6-state logic,
% extended with a sane bound on b_q/b_r so a runaway gyro-bias
% estimate can't silently poison future steps the way the unclamped
% q/r divergence did before the original fix) ---
if abs(x_new(1)) >= deg2rad(5) || abs(x_new(2)) >= deg2rad(5)
x_new(3) = max(min(x_new(3), deg2rad(50)), deg2rad(-50)); % q
x_new(4) = max(min(x_new(4), deg2rad(50)), deg2rad(-50)); % r
end
x_new(1) = max(min(x_new(1), deg2rad(5)), deg2rad(-5)); % dtheta_S
x_new(2) = max(min(x_new(2), deg2rad(5)), deg2rad(-5)); % dpsi_S
x_new(5) = max(min(x_new(5), 2*pi*250), 0); % pb: 0-250 Hz
% Gyro bias clamp: placeholder bound, ASSUMPTION not datasheet-derived.
% Set wide enough not to interfere in normal operation; tighten once you
% have real bias-instability numbers for your gyro.
x_new(7) = max(min(x_new(7), deg2rad(20)), deg2rad(-20)); % b_q
x_new(8) = max(min(x_new(8), deg2rad(20)), deg2rad(-20)); % b_r

end

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function h_axial = mag_axial_h_StageE8(x, phi_roll, theta_v, psi_v, B_earth)
% Identical to the 6-state mag_axial_h_StageE -- axial mag model has no
% dependence on states 7-8 (b_q, b_r), so this is unchanged except that
% x is now 8x1 (only indices 1,2,6 are read).
dtheta_S = x(1);
dpsi_S = x(2);
axial_bias = x(6);

theta_pred = theta_v + dtheta_S;
psi_pred = psi_v + dpsi_S;

Ry = [cos(theta_pred), 0, -sin(theta_pred); 0, 1, 0; sin(theta_pred), 0, cos(theta_pred)];
Rz = [cos(psi_pred), sin(psi_pred), 0; -sin(psi_pred), cos(psi_pred), 0; 0, 0, 1];
Rx = [1, 0, 0; 0, cos(phi_roll), -sin(phi_roll); 0, sin(phi_roll), cos(phi_roll)];

R_i2b_pred = Rx * Ry * Rz;
B_body_pred = R_i2b_pred * B_earth;
h_axial = B_body_pred(1) + axial_bias;
end

function H = mag_axial_jacobian_numeric8(x, phi_roll, theta_v, psi_v, B_earth)
% Numeric Jacobian of mag_axial_h_StageE8 w.r.t. all 8 states.
% Columns 7,8 (b_q, b_r) are expected to come out exactly 0, since the
% axial mag model doesn't touch them -- this is the intended structural
% independence, mirrored from Phase 1's cross-state independence test.
n = length(x);
H = zeros(1, n);
eps_ang = 1e-6;
eps_bias = 1e-3;
eps_rate = 1e-6;
eps_gyro_bias = 1e-6;

h0 = mag_axial_h_StageE8(x, phi_roll, theta_v, psi_v, B_earth);
for i = 1:n
dx = zeros(n,1);
if i == 6
dx(i) = eps_bias;
elseif i == 1 || i == 2
dx(i) = eps_ang;
elseif i == 7 || i == 8
dx(i) = eps_gyro_bias;
else
dx(i) = eps_rate; % q, r, pb -- don't enter h directly, exactly 0
end
h1 = mag_axial_h_StageE8(x + dx, phi_roll, theta_v, psi_v, B_earth);
H(i) = (h1 - h0) / dx(i);
end
end

function F = numeric_F_StageE8(x, p, qbar, V, dt)
n = length(x);
F = eye(n);
eps_val = 1e-6;
f0 = stateE_dot8(x, p, qbar, V);
for i = 1:n
dx = zeros(n,1);
dx(i) = eps_val;
f1 = stateE_dot8(x + dx, p, qbar, V);
F(:,i) = F(:,i) + dt * (f1 - f0) / eps_val;
end
end

function xdot = stateE_dot8(x, p, qbar, V)
dtheta_S = x(1); dpsi_S = x(2); q = x(3); r = x(4); pb = x(5);
% axial_bias (x6), b_q (x7), b_r (x8) are all constant -> dot = 0
L_B = qbar * p.A * p.d * (p.Clp * pb * p.d / (2*V));
M_B = qbar * p.A * p.d * (p.Cma*dtheta_S + p.Cmq*q*p.d/(2*V) + p.Cnpa*pb*p.d/(2*V)*dpsi_S);
N_B = qbar * p.A * p.d * (-p.Cma*dpsi_S + p.Cmq*r*p.d/(2*V) + p.Cnpa*pb*p.d/(2*V)*dtheta_S);
xdot = [q; r; pb*r*(p.Ix/p.Iy)+M_B/p.Iy; -pb*q*(p.Ix/p.Iy)+N_B/p.Iy; L_B/p.Ix; 0; 0; 0];
end

function h_gyro = gyro_h_StageE8(x, phi_roll)
% CORRECTED (Phase 1b): gyro measures BODY-frame rates; state's q, r
% (indices 3,4) are S-frame (fixed-plane) rates. S-frame and body frame
% share the same x/roll axis and differ only by rotation through
% phi_roll (Stage D, external input) -- same Rx(phi_roll) convention
% used in mag_axial_h_StageE8 / Layer 6A / Layer 6C's R_i2b.
%
% [q_body] [ cos(phi_roll) -sin(phi_roll) ] [q_S]
% [r_body] = [ sin(phi_roll) cos(phi_roll) ] [r_S]
% h_gyro = [q_body_pred + b_q; r_body_pred + b_r]
%
% At phi_roll=0 this collapses exactly to the old (Phase 1) unrotated
% model -- validated in phase1b_gyro_h_corrected.m, Test 1.
q_S = x(3);
r_S = x(4);
b_q = x(7);
b_r = x(8);

q_body_pred = cos(phi_roll)*q_S - sin(phi_roll)*r_S;
r_body_pred = sin(phi_roll)*q_S + cos(phi_roll)*r_S;

h_gyro = [q_body_pred + b_q; r_body_pred + b_r];
end

function H = gyro_jacobian_numeric8(x, phi_roll)
% Numeric Jacobian of the corrected gyro_h_StageE8 w.r.t. all 8 states,
% at fixed (externally supplied) phi_roll -- same treatment as
% mag_axial_jacobian_numeric8 holding theta_v/psi_v fixed. Validated in
% phase1b_gyro_h_corrected.m (Test 4: q/r block matches the exact
% rotation matrix; other columns remain exactly zero; bias columns
% remain clean unity/zero).
n = length(x);
H = zeros(2, n);
eps_ang = 1e-6;
eps_rate = 1e-6;
eps_bias_mag = 1e-3;
eps_bias_gyro = 1e-6;

h0 = gyro_h_StageE8(x, phi_roll);
for i = 1:n
dx = zeros(n,1);
if i == 1 || i == 2
dx(i) = eps_ang;
elseif i == 3 || i == 4 || i == 5
dx(i) = eps_rate;
elseif i == 6
dx(i) = eps_bias_mag;
else
dx(i) = eps_bias_gyro;
end
h1 = gyro_h_StageE8(x + dx, phi_roll);
H(:,i) = (h1 - h0) / dx(i);
end
end

function R_gyro_i = gyro_R_scaling(pb_i, p)
% Same as Phase 2's gyro_R_scaling -- restated here so this file is
% self-contained. See Phase 2 file for the full modeling-assumption
% documentation (NOT datasheet-derived, tune R_floor/R_scale/n).
R_floor = p.sigma_gyro^2;
R_scale = p.sigma_gyro_spin^2;
pb_ref = p.p0;
n = p.gyro_R_exponent;
R_gyro_i = R_floor + R_scale * (abs(pb_i)/pb_ref)^n;
end
