% =========================================================================
% UNIFIED CAUSAL PIPELINE
% Patent US 9,702,674 B2 -- Spin-Stabilized Projectile Navigation
%
% Merges the previously-separate Layer 7/6A script and Layer 6C fusion
% script into a single, time-interleaved, causal loop:
%
% [Stage E Step] -> [Layer 6A Step] -> [Track1 & Track2 Propagation]
% -> [Layer 6C Combiner Update] -> [GPS Correction Anchor Layer]
%
% GROUND-TRUTH LEAKAGE FIX:
% Previously, Stage E computed its flight-geometry reference angles
% (theta_v, psi_v) from true_vel_log(i,:) (or GPS-corrected truth).
% This is now replaced with the FUSED velocity from the PREVIOUS
% timestep, vel_c_log(i-1,:), closing the loop causally via a unit
% delay (z^-1). At t=0 (i=1) the loop bootstraps from the honest,
% gun-derived initial velocity [u0, v0, w0].
%
% KNOWN, DELIBERATE, FLAGGED LIMITATIONS (Phase 2 scope):
% - Stage C (mag de-biasing geometry: R_field_i) and Stage D
% (build_S_frame's chi_i) STILL use true_vel_log for their internal
% geometry construction. This is an intentional, explicit scope
% limitation for this pass -- not an oversight. A fully causal
% version would need Stage C/D to also consume vel_c_log(i-1,:),
% but that reworking is deferred.
% - accel_clean_log(1,:) is never populated (loop starts at i=2, same
% as the prior separate-script version) -- Track2's first inertial
% step at i=2 therefore integrates from a zero accel sample at i=1,
% identical behavior to before. Not introduced by this merge.
% =========================================================================

clc; clear all; close all;
p = projectile_params();

% ================================================
% RUN HIFI TRUTH MODEL (unchanged)
% ================================================
elev_rad = deg2rad(p.elev);
az_rad = deg2rad(p.az);
u0 = p.V0*cos(elev_rad)*cos(az_rad);
v0 = p.V0*cos(elev_rad)*sin(az_rad);
w0 = -p.V0*sin(elev_rad);
phi0 = 0; theta0 = elev_rad; psi0 = az_rad;

q0i = cos(theta0/2); q1i = 0; q2i = sin(theta0/2); q3i = 0;
state_hifi = [0;0;0; u0;v0;w0; q0i;q1i;q2i;q3i; p.p0;0;0];
phi_c = 0;

N = round(p.t_max/p.dt) + 1;
time_log = zeros(N,1);
true_log = zeros(N,4);
true_vel_log = zeros(N,3);
true_alt_log = zeros(N,1);
mag_log = zeros(N,3);
accel_log = zeros(N,3);
accel_true_log = zeros(N,3);
accel_NP_log = zeros(N,3);
gps_vel_log = NaN(N,3);
true_pos = zeros(N,3);

k = 1; t = 0;
true_log(1,:) = [phi0, theta0, psi0, p.p0];
true_vel_log(1,:) = [u0, v0, w0];
true_alt_log(1) = 0;
true_pos(1,:) = [0, 0, 0];

while t < p.t_max
altitude = -state_hifi(3);
if t > 0.1 && altitude <= 0; break; end
[rho,~] = atmosphere(max(altitude,0));
k1 = eom_hifi(state_hifi, p, rho, phi_c);
k2 = eom_hifi(state_hifi+0.5*p.dt*k1, p, rho, phi_c);
k3 = eom_hifi(state_hifi+0.5*p.dt*k2, p, rho, phi_c);
k4 = eom_hifi(state_hifi+p.dt*k3, p, rho, phi_c);
state_hifi = state_hifi + (p.dt/6)*(k1+2*k2+2*k3+k4);
qn = norm(state_hifi(7:10));
if qn > 0; state_hifi(7:10) = state_hifi(7:10)/qn; end
q0=state_hifi(7); q1=state_hifi(8); q2=state_hifi(9); q3=state_hifi(10);
phi_e = atan2(2*(q0*q1+q2*q3), 1-2*(q1^2+q2^2));
theta_e = asin(max(-1,min(1, 2*(q0*q2-q3*q1))));
psi_e = atan2(2*(q0*q3+q1*q2), 1-2*(q2^2+q3^2));
sensors = sensor_models_hifi(state_hifi, p, t);
t = t + p.dt; k = k + 1;
time_log(k) = t;
true_log(k,:) = [phi_e, theta_e, psi_e, state_hifi(11)];
true_vel_log(k,:) = state_hifi(4:6)';
true_alt_log(k) = altitude;
true_pos(k,:) = true_pos(k-1,:) + p.dt*0.5*(true_vel_log(k,:) + true_vel_log(k-1,:));
mag_log(k,:) = sensors.mag_body';
accel_log(k,:) = sensors.accel_body';
accel_true_log(k,:) = sensors.accel_body_true';
accel_NP_log(k,:) = sensors.accel_NP';
if sensors.gps_available
gps_vel_log(k,:) = sensors.gps_vel';
end
end
time_log = time_log(1:k);
true_log = true_log(1:k,:);
true_vel_log = true_vel_log(1:k,:);
true_alt_log = true_alt_log(1:k);
true_pos = true_pos(1:k,:);
mag_log = mag_log(1:k,:);
accel_log = accel_log(1:k,:);
accel_true_log = accel_true_log(1:k,:);
accel_NP_log = accel_NP_log(1:k,:);
gps_vel_log = gps_vel_log(1:k,:);
N2 = k;
fprintf('Truth model: %d steps, %.2f s\n', N2, t);

% ================================================
% STAGE C: MAGNETOMETER DE-BIASING
% NOTE (Phase 2 limitation): still uses true_vel_log for R_field_i
% geometry. Not yet closed-loop.
% ================================================
by_est_log = zeros(N2,1);
bz_est_log = zeros(N2,1);
x_debias = [0;0];
P_debias = diag([500^2, 500^2]);
Q_debias = diag([0.0005, 0.0005]);
mag_y = mag_log(:,2);
mag_z = mag_log(:,3);
for i = 1:N2
vel_i = true_vel_log(i,:)';
psi_v = atan2(vel_i(2), vel_i(1));
theta_v = atan2(-vel_i(3), sqrt(vel_i(1)^2+vel_i(2)^2));
R_i2f_i = dcm(psi_v, theta_v);
B_f_i = R_i2f_i * p.B_earth;
R_field_i = norm(B_f_i(2:3));
R_noise_i = (2*R_field_i*p.sigma_mag)^2;
[x_debias, P_debias] = mag_debias_filter_step(x_debias, P_debias, ...
mag_y(i), mag_z(i), R_field_i, Q_debias, R_noise_i);
by_est_log(i) = x_debias(1);
bz_est_log(i) = x_debias(2);
end
mag_y_clean = mag_y - by_est_log;
mag_z_clean = mag_z - bz_est_log;
fprintf('Stage C: bias = [%.1f, %.1f] nT (true [%.0f, %.0f])\n', ...
by_est_log(end), bz_est_log(end), p.bias_mag(2), p.bias_mag(3));

% ================================================
% STAGE D: ROLL EXTRACTION
% NOTE (Phase 2 limitation): build_S_frame's chi_i still uses
% true_vel_log. Not yet closed-loop.
% ================================================
phi_est_stageD = zeros(N2,1);
for i = 1:N2
vel_i = true_vel_log(i,:)';
[~, chi_i] = build_S_frame(vel_i, p.B_earth);
phi_raw = atan2(mag_y_clean(i), mag_z_clean(i)) - chi_i;
phi_est_stageD(i) = atan2(sin(phi_raw), cos(phi_raw));
end
phi_err_stageD = rad2deg(wrapToPi(true_log(:,1) - phi_est_stageD));
fprintf('Stage D: RMS phi error = %.4f deg\n', rms(phi_err_stageD));

% ================================================
% STAGE 2: pb measurement (zero-crossing) -- independent of velocity
% ================================================
idx_cal = find(mag_y(1:end-1) < 0 & mag_y(2:end) >= 0, 1, 'first');
phi_ref_offset = true_log(idx_cal,1);
[~, pb_est_stage2] = roll_spin_extractor(mag_y, time_log, phi_ref_offset);

% ================================================
% STAGE E / LAYER 6A / TRACK1 / TRACK2 / LAYER 6C -- STATE INIT
% ================================================

% --- Stage E state ---
x_E = [0; 0; 0; 0; p.p0; 0];
P_E = diag([deg2rad(2)^2, deg2rad(2)^2, deg2rad(1)^2, deg2rad(1)^2, ...
(2*pi*20)^2, 300^2]);
Q_E = diag([deg2rad(0.05)^2, deg2rad(0.05)^2, ...
deg2rad(0.001)^2, deg2rad(0.001)^2, ...
(2*pi*0.5)^2, 0.001^2]);
R_pb = (2*pi*1.2)^2;

theta_est_E = zeros(N2,1);
psi_est_E = zeros(N2,1);
pb_est_E = zeros(N2,1);
axial_bias_est_E = zeros(N2,1);
q_est_E = zeros(N2,1);
r_est_E = zeros(N2,1);
dtheta_S_log = zeros(N2,1);
dpsi_S_log = zeros(N2,1);

theta_est_E(1) = theta0;
psi_est_E(1) = psi0;
pb_est_E(1) = p.p0;

P_log_E = zeros(6,6,N2);
P_log_E(:,:,1) = P_E;

% --- Layer 6A (NP estimator) state ---
sigma_b_init = 0.5 * p.g;
sigma_r_init = 0.005;
x_np_rad = [0; 0; 0.002; 0.0015];
P_np_rad = diag([sigma_b_init^2, sigma_b_init^2, sigma_r_init^2, sigma_r_init^2]);
Q_np_rad = diag([(0.001*p.g)^2, (0.001*p.g)^2, (1e-6)^2, (1e-6)^2]);
x_np_ax = 0;
P_np_ax = (5*p.g)^2;
Q_np_ax = (0.001*p.g)^2;
R_np_ax = (p.sigma_accel*p.g)^2;

NP_est_log = zeros(N2, 3);
accel_clean_log = zeros(N2, 3);
r_offset_log = zeros(N2, 2);

% --- Track1 (aero/dynamic) state ---
pos_d_state = [0; 0; 0]; vel_d_state = [u0; v0; w0];
% --- Track2 (inertial) state ---
pos_a_state = [0; 0; 0]; vel_a_state = [u0; v0; w0];

P_initial = diag([1.0^2, 1.0^2, 1.0^2, 1.0^2, 1.0^2, 1.0^2]);
P_d_state = P_initial + 1e-3*eye(6);
P_a_state = P_initial + 2e-3*eye(6);
P_da_state = P_initial * 0.98;

sigma_a_model = 0.1; sigma_pos_model = 0.01;
Q_d = diag([sigma_pos_model^2, sigma_pos_model^2, sigma_pos_model^2, ...
(sigma_a_model*p.dt)^2, (sigma_a_model*p.dt)^2, (sigma_a_model*p.dt)^2]);
sigma_a_inertial = 0.5; sigma_pos_inertial = 0.01;
Q_a = diag([sigma_pos_inertial^2, sigma_pos_inertial^2, sigma_pos_inertial^2, ...
(sigma_a_inertial*p.dt)^2, (sigma_a_inertial*p.dt)^2, (sigma_a_inertial*p.dt)^2]);
Q_da = zeros(6,6);

P_DA_EIG_CAP = 1e4;

R_gps = diag([repmat(5.0^2, 1, 3), repmat(0.4^2, 1, 3)]);
H_gps = eye(6);
next_gps_time = 3.0;
gps_interval = 4.0;

pos_c_log = zeros(N2, 3); vel_c_log = zeros(N2, 3);
L_norm_log = zeros(N2, 1);
denom_mineig_log = zeros(N2, 1);
P_da_norm_log = zeros(N2,1);
P_da_cap_triggered_log = false(N2,1);
gps_active_log = false(N2,1);

% --- Standalone Track1 (aero) and Track2 (INS) history logs ---
% These record each track's OWN state after every step (including at
% GPS resets, when both are anchored to x_c_corrected -- same as the
% combiner's internal states). This lets us see how each track performs
% on its own, not just the cascaded combiner output.
pos_d_log = zeros(N2, 3); vel_d_log = zeros(N2, 3);
pos_a_log = zeros(N2, 3); vel_a_log = zeros(N2, 3);
pos_d_log(1,:) = pos_d_state'; vel_d_log(1,:) = vel_d_state';
pos_a_log(1,:) = pos_a_state'; vel_a_log(1,:) = vel_a_state';

% --- Bootstrap: t=0 fused state = honest gun-derived initial velocity ---
pos_c_log(1,:) = pos_d_state';
vel_c_log(1,:) = vel_d_state'; % = [u0, v0, w0]

% ================================================
% UNIFIED CAUSAL TIME LOOP
% [Stage E] -> [Layer 6A] -> [Track1 & Track2] -> [Layer 6C] -> [GPS]
% ================================================
for i = 2:N2

[rho_i,~] = atmosphere(max(0, true_alt_log(i)));

% --- CAUSAL VELOCITY REFERENCE (ground-truth leak fix) ---
% Stage E's flight-geometry angles now come from the PREVIOUS
% timestep's FUSED velocity estimate (unit delay z^-1), not from
% true_vel_log or gps_vel_log.
vel_est = vel_c_log(i-1,:)';

V_i = norm(vel_est);
if V_i < 1; V_i = 1; end
qbar_i = 0.5*rho_i*V_i^2;

pb_meas_i = pb_est_stage2(i);
if isnan(pb_meas_i); pb_meas_i = x_E(5); end

if time_log(i) < 3.0
R_axial_i = (10 * p.sigma_mag)^2;
else
R_axial_i = (p.sigma_mag)^2;
end

% ---------------- STAGE E STEP ----------------
[x_E, P_E] = attitude_ekf_step_StageE(x_E, P_E, mag_log(i,:)', ...
phi_est_stageD(i), vel_est, p, qbar_i, V_i, p.dt, Q_E, R_axial_i, ...
pb_meas_i, R_pb);
P_log_E(:,:,i) = P_E;

psi_v_i = atan2(vel_est(2), vel_est(1));
theta_v_i = atan2(-vel_est(3), sqrt(vel_est(1)^2+vel_est(2)^2));
theta_est_E(i) = theta_v_i + x_E(1);
psi_est_E(i) = psi_v_i + x_E(2);
pb_est_E(i) = x_E(5);
axial_bias_est_E(i) = x_E(6);
q_est_E(i) = x_E(3);
r_est_E(i) = x_E(4);
dtheta_S_log(i) = x_E(1);
dpsi_S_log(i) = x_E(2);

% ---------------- LAYER 6A STEP ----------------
phi_roll_i = phi_est_stageD(i);
theta_est_i = theta_est_E(i);
psi_est_i = psi_est_E(i);
Ry_i = [cos(theta_est_i), 0, -sin(theta_est_i); 0, 1, 0; sin(theta_est_i), 0, cos(theta_est_i)];
Rz_i = [cos(psi_est_i), sin(psi_est_i), 0; -sin(psi_est_i), cos(psi_est_i), 0; 0, 0, 1];
Rx_i = [1, 0, 0; 0, cos(phi_roll_i), -sin(phi_roll_i); 0, sin(phi_roll_i), cos(phi_roll_i)];
R_i2b_est = Rx_i * Ry_i * Rz_i;
vel_b_est = R_i2b_est * vel_est;
u_b = vel_b_est(1); v_b = vel_b_est(2); w_b = vel_b_est(3);
alpha_i = atan2(w_b, u_b);
beta_i = atan2(v_b, u_b);
alpha_total_i = sqrt(alpha_i^2 + beta_i^2);

f_x_pred = -(qbar_i * p.A * (p.CD0 + p.CX_a2*alpha_total_i^2)) / p.m;
f_y_pred = (qbar_i * p.A * (-p.CNa*beta_i + p.CYpa*pb_est_E(i)*p.d/(2*V_i))) / p.m;
f_z_pred = (qbar_i * p.A * (-p.CNa*alpha_i)) / p.m;
f_aero_pred = [f_x_pred; f_y_pred; f_z_pred];

sigma_ball_j = p.balloting_j * exp(-time_log(i)/p.balloting_tau);
sigma_ball_k = p.balloting_k * exp(-time_log(i)/p.balloting_tau);
R_np_rad_i = diag([(p.sigma_accel*p.g + sigma_ball_j)^2, ...
(p.sigma_accel*p.g + sigma_ball_k)^2]);

[x_np_rad, P_np_rad, x_np_ax, P_np_ax, accel_clean_i] = ...
accel_NP_estimator_step(x_np_rad, P_np_rad, x_np_ax, P_np_ax, ...
accel_log(i,:)', f_aero_pred, pb_est_E(i), ...
Q_np_rad, R_np_rad_i, Q_np_ax, R_np_ax);

pb2_i = pb_est_E(i)^2;
NP_j_est = x_np_rad(1) + x_np_rad(3)*pb2_i;
NP_k_est = x_np_rad(2) + x_np_rad(4)*pb2_i;
NP_est_log(i,:) = [x_np_ax; NP_j_est; NP_k_est];
accel_clean_log(i,:) = accel_clean_i;
r_offset_log(i,:) = [x_np_rad(3), x_np_rad(4)];

% ---------------- TRACK1 (AERO/DYNAMIC) PROPAGATION ----------------
[pos_d, vel_d, Phi_d, P_d] = track1_dynamic_step(...
pos_d_state, vel_d_state, P_d_state, ...
dtheta_S_log(i), dpsi_S_log(i), ...
pb_est_E(i), q_est_E(i), r_est_E(i), ...
p.B_earth, rho_i, p.dt, p, Q_d);
x_d = [pos_d; vel_d];

% ---------------- TRACK2 (INS) PROPAGATION ----------------
[pos_a, vel_a, Phi_a, P_a] = track2_inertial_step(...
pos_a_state, vel_a_state, P_a_state, ...
accel_clean_log(i-1,:)', accel_clean_log(i,:)', ...
phi_est_stageD(i-1), theta_est_E(i-1), psi_est_E(i-1), ...
phi_est_stageD(i), theta_est_E(i), psi_est_E(i), ...
p.dt, p, Q_a);
x_a = [pos_a; vel_a];

% ---------------- LAYER 6C: CASCADED FUSION COMBINER ----------------
P_da = Phi_d * P_da_state * Phi_a' + Q_da;

% Hard eigenvalue cap on P_da every step (structural safety net).
[V_da, D_da] = eig(P_da);
d_da = diag(D_da);
if any(d_da > P_DA_EIG_CAP)
d_da = min(d_da, P_DA_EIG_CAP);
P_da = V_da * diag(d_da) * V_da';
P_da = (P_da + P_da')/2;
P_da_cap_triggered_log(i) = true;
end
P_da_norm_log(i) = norm(P_da, 'fro');

cross_da = P_da;
cross_ad = P_da';
numerator = P_d - cross_da;
denominator = P_a + P_d - cross_ad - cross_da;
reg = diag(diag(denominator) * 1e-8);
reg(reg < 1e-10) = 1e-10;
L = numerator / (denominator + reg);

x_c_minus = x_d + L * (x_a - x_d);
P_c_minus = P_d - L * (P_d - cross_ad);
P_c_minus = (P_c_minus + P_c_minus') / 2;

L_norm_log(i) = norm(L, 'fro');
denom_mineig_log(i) = min(eig(denominator));

% ---------------- GPS CORRECTION ANCHOR LAYER ----------------
t = time_log(i);
is_gps_available = (t <= 40.0) && (t >= next_gps_time);
gps_active_log(i) = is_gps_available;

if is_gps_available
next_gps_time = t + gps_interval;
z_gps = [true_pos(i, :)'; true_vel_log(i, :)'] + [5.0*randn(3,1); 0.4*randn(3,1)];
y_residual = z_gps - H_gps * x_c_minus;
S_residual = H_gps * P_c_minus * H_gps' + R_gps;
K_gps = P_c_minus * H_gps' / S_residual;
x_c_corrected = x_c_minus + K_gps * y_residual;
I_mat = eye(6);
P_c_corrected = (I_mat - K_gps*H_gps) * P_c_minus * (I_mat - K_gps*H_gps)' + K_gps*R_gps*K_gps';
P_c_corrected = (P_c_corrected + P_c_corrected') / 2;

% Independent-track states anchored to corrected fused state.
pos_d_state = x_c_corrected(1:3); vel_d_state = x_c_corrected(4:6);
pos_a_state = x_c_corrected(1:3); vel_a_state = x_c_corrected(4:6);

% Distinct perturbations so the denominator never restarts
% exactly singular right after a reset.
eps_d = 1e-3;
eps_a = 2e-3;
P_d_state = P_c_corrected + eps_d*eye(6);
P_a_state = P_c_corrected + eps_a*eye(6);
P_da_state = P_c_corrected * 0.98;

pos_c_log(i,:) = x_c_corrected(1:3)';
vel_c_log(i,:) = x_c_corrected(4:6)';
else
% Tracks propagate independently -- NOT overwritten from x_c.
pos_d_state = pos_d; vel_d_state = vel_d; P_d_state = P_d;
pos_a_state = pos_a; vel_a_state = vel_a; P_a_state = P_a;
P_da_state = P_da;

pos_c_log(i,:) = x_c_minus(1:3)';
vel_c_log(i,:) = x_c_minus(4:6)';
end

% Log standalone track states AFTER GPS resolution above, so a reset
% at this step is reflected in both the combiner's internal state
% and these standalone histories (matches what pos_d_state/pos_a_state
% actually carry into the NEXT iteration's propagation).
pos_d_log(i,:) = pos_d_state'; vel_d_log(i,:) = vel_d_state';
pos_a_log(i,:) = pos_a_state'; vel_a_log(i,:) = vel_a_state';
end

% ================================================
% RESULTS
% ================================================
theta_err_E = rad2deg(wrapToPi(true_log(:,2) - theta_est_E));
psi_err_E = rad2deg(wrapToPi(true_log(:,3) - psi_est_E));
pb_err_E = (true_log(:,4) - pb_est_E) / (2*pi);

fprintf('\n=== Stage E (causal, closed-loop) ===\n');
fprintf('RMS theta = %.3f deg, RMS psi = %.3f deg\n', rms(theta_err_E(2:end)), rms(psi_err_E(2:end)));
fprintf('RMS pb = %.3f Hz\n', rms(pb_err_E(2:end)));
fprintf('Final axial bias = %.1f nT (true: %.0f)\n', axial_bias_est_E(end), p.bias_mag(1));

NP_err = accel_NP_log - NP_est_log;
idx_settled = time_log > 3*p.balloting_tau;
fprintf('\n=== Layer 6A: NP Estimator ===\n');
fprintf(' Axial (X): %.4f m/s^2\n', rms(NP_err(idx_settled,1)));
fprintf(' Radial (Y): %.4f m/s^2\n', rms(NP_err(idx_settled,2)));
fprintf(' Radial (Z): %.4f m/s^2\n', rms(NP_err(idx_settled,3)));

vel_err_c = true_vel_log - vel_c_log;
pos_err_c = true_pos - pos_c_log;
fprintf('\n=========================================================\n');
fprintf('FINAL FUSION RESULTS (causal, closed-loop Stage E)\n');
fprintf('=========================================================\n');
fprintf('CASCADED COMBINED Velocity RMS: [%.2f, %.2f, %.2f] m/s\n', ...
rms(vel_err_c(2:end,1)), rms(vel_err_c(2:end,2)), rms(vel_err_c(2:end,3)));
fprintf('CASCADED COMBINED Final position error: [%.1f, %.1f, %.1f] m\n\n', ...
pos_err_c(end,1), pos_err_c(end,2), pos_err_c(end,3));
fprintf('L norm: min=%.4f max=%.4f\n', min(L_norm_log(2:end)), max(L_norm_log(2:end)));
fprintf('Denominator min-eigenvalue: min=%.4e max=%.4e\n', ...
min(denom_mineig_log(2:end)), max(denom_mineig_log(2:end)));
fprintf('P_da cap triggered on %.1f%% of steps\n\n', 100*mean(P_da_cap_triggered_log(2:end)));

[worst_L, worst_idx] = max(L_norm_log(2:end));
fprintf('Worst L norm (%.2f) occurs at t=%.2fs (GPS active at that time: %d)\n\n', ...
worst_L, time_log(worst_idx+1), gps_active_log(worst_idx+1));

% --- Individual track performance (independent of combiner) ---
vel_err_d = true_vel_log - vel_d_log;
pos_err_d = true_pos - pos_d_log;
vel_err_a = true_vel_log - vel_a_log;
pos_err_a = true_pos - pos_a_log;

fprintf('=========================================================\n');
fprintf('STANDALONE TRACK PERFORMANCE (vs cascaded combiner above)\n');
fprintf('=========================================================\n');
fprintf('TRACK1 (aero/dynamic) Velocity RMS: [%.2f, %.2f, %.2f] m/s\n', ...
rms(vel_err_d(2:end,1)), rms(vel_err_d(2:end,2)), rms(vel_err_d(2:end,3)));
fprintf('TRACK1 Final position error: [%.1f, %.1f, %.1f] m\n\n', ...
pos_err_d(end,1), pos_err_d(end,2), pos_err_d(end,3));
fprintf('TRACK2 (INS) Velocity RMS: [%.2f, %.2f, %.2f] m/s\n', ...
rms(vel_err_a(2:end,1)), rms(vel_err_a(2:end,2)), rms(vel_err_a(2:end,3)));
fprintf('TRACK2 Final position error: [%.1f, %.1f, %.1f] m\n\n', ...
pos_err_a(end,1), pos_err_a(end,2), pos_err_a(end,3));

% --- Time-resolved diagnostic plots ---
figure('Name','Time-Resolved Fusion Diagnostics (Causal Stage E)','Position',[100 100 1200 900]);
subplot(4,1,1);
plot(time_log, L_norm_log); grid on; title('||L||_F over time'); ylabel('norm');
hold on; xline(40,'r--','GPS cutoff');
subplot(4,1,2);
plot(time_log, denom_mineig_log); grid on; title('Denominator min eigenvalue over time');
set(gca,'YScale','log'); ylabel('eig'); hold on; xline(40,'r--');
subplot(4,1,3);
plot(time_log, P_da_norm_log); grid on; title('||P_{da}||_F over time (capped)');
ylabel('norm'); hold on; xline(40,'r--');
subplot(4,1,4);
plot(time_log, gps_active_log,'g-'); grid on; title('GPS active flag');
xlabel('t [s]'); ylim([-0.1 1.1]);

% --- Trajectory comparison plots: True vs Combined vs Track1 vs Track2 ---
figure('Name','Trajectory: True vs Fused vs Independent Tracks','Position',[150 150 1200 800]);
labels = {'X','Y','Z'};
for ax = 1:3
subplot(3,1,ax);
plot(time_log, true_pos(:,ax), 'k-', 'LineWidth', 1.4); hold on;
plot(time_log, pos_c_log(:,ax), 'b-', 'LineWidth', 1.1);
plot(time_log, pos_d_log(:,ax), 'r--', 'LineWidth', 0.9);
plot(time_log, pos_a_log(:,ax), 'g:', 'LineWidth', 1.1);
grid on; ylabel(sprintf('Pos %s [m]', labels{ax}));
if ax == 1; legend('True','Combined (6C)','Track1 (aero)','Track2 (INS)','Location','best'); end
end
xlabel('t [s]');
sgtitle('Position: True vs Cascaded-Fused vs Independent Tracks');

figure('Name','Velocity: True vs Fused vs Independent Tracks','Position',[200 200 1200 800]);
for ax = 1:3
subplot(3,1,ax);
plot(time_log, true_vel_log(:,ax), 'k-', 'LineWidth', 1.4); hold on;
plot(time_log, vel_c_log(:,ax), 'b-', 'LineWidth', 1.1);
plot(time_log, vel_d_log(:,ax), 'r--', 'LineWidth', 0.9);
plot(time_log, vel_a_log(:,ax), 'g:', 'LineWidth', 1.1);
grid on; ylabel(sprintf('Vel %s [m/s]', labels{ax}));
if ax == 1; legend('True','Combined (6C)','Track1 (aero)','Track2 (INS)','Location','best'); end
end
xlabel('t [s]');
sgtitle('Velocity: True vs Cascaded-Fused vs Independent Tracks');

% --- Position error comparison (all three vs truth, per axis) ---
figure('Name','Position Error: Combined vs Track1 vs Track2','Position',[250 150 1200 800]);
for ax = 1:3
subplot(3,1,ax);
plot(time_log, pos_err_c(:,ax), 'b-', 'LineWidth', 1.1); hold on;
plot(time_log, pos_err_d(:,ax), 'r--', 'LineWidth', 0.9);
plot(time_log, pos_err_a(:,ax), 'g:', 'LineWidth', 1.1);
grid on; ylabel(sprintf('Pos Err %s [m]', labels{ax}));
if ax == 1; legend('Combined (6C)','Track1 (aero)','Track2 (INS)','Location','best'); end
end
xlabel('t [s]');
sgtitle('Position Error vs Truth: Combined vs Independent Tracks');

fprintf('\nNOTE: P_log_E (6x6xN2 Stage E covariance history) is available in the workspace.\n');
fprintf('NOTE: Stage C/D geometry (R_field_i, chi_i) still consumes true_vel_log --\n');
fprintf('flagged Phase 2 limitation, intentionally out of scope for this merge.\n');
