% =========================================================================
% LAYER 7 v4 + Layer 6A: Full pipeline -- 8-STATE STAGE E (gyro-augmented)
% Truth -> Stage C -> Stage D -> Stage E (8-state, +gyro) -> NP Estimator
% Patent US 9,702,674 B2, FIG.7 + FIG.9
%
% CHANGES FROM PREVIOUS (6-state, axial-mag-only) VERSION:
% 1. Stage E state grows 6 -> 8: adds gyro rate biases [b_q; b_r].
% 2. A third measurement update is added each step: direct gyro rate
% readings [q_body; r_body], via attitude_ekf_step_StageE_v8.m.
% 3. gyro_meas is synthesized in the truth-model loop from eom_hifi's
% TRUE body rates q_b, r_b (state_hifi(12), state_hifi(13)) +
% assumed constant bias + noise -- NOT from the fixed-plane
% dtheta_S_dot/dpsi_S_dot, per the body-vs-S-frame discussion.
% 4. R_gyro is spin-adaptive (Phase 2: gyro_R_scaling), inflating at
% high pb so the filter structurally distrusts the gyro during
% high-spin/saturation-risk conditions and leans on it more as
% spin decays -- this is NOT yet validated against real saturation
% behavior, just a documented modeling assumption (see Phase 2).
%
% STILL UNCHANGED / OUT OF SCOPE (same limitations as before):
% - Stage C, Stage D still use true_vel_log for internal geometry.
% - Radial magnetometer Stage E update remains REVERTED (circular
% dependency with Stage D, see prior history) -- not reintroduced.
% =========================================================================
clc; clear all; close all;
p = projectile_params();

% --- Gyro model parameters (ASSUMPTIONS -- no datasheet, see Phase 2) ---
p.sigma_gyro = deg2rad(0.05); % rad/s, gyro noise floor (1-sigma) at zero spin
p.sigma_gyro_spin = deg2rad(0.5); % rad/s, additional 1-sigma noise AT pb_ref (=p.p0)
p.gyro_R_exponent = 2; % n in R_gyro(pb) = R_floor + R_scale*(pb/pb_ref)^n
p.bias_gyro_q = deg2rad(0.3); % assumed constant true gyro bias, q-channel [rad/s]
p.bias_gyro_r = deg2rad(-0.2); % assumed constant true gyro bias, r-channel [rad/s]
p.sigma_gyro_meas = deg2rad(0.05); % 1-sigma measurement noise added when SYNTHESIZING
% gyro_meas in the truth loop below (kept separate
% from sigma_gyro/sigma_gyro_spin, which define the
% FILTER's assumed R -- these can be tuned to differ
% if you want to test filter mismatch robustness)

% ================================================
% RUN HIFI TRUTH MODEL
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
true_qr_log = zeros(N,2); % NEW: true body-frame [q_b, r_b], for gyro synthesis + validation
gyro_meas_log = zeros(N,2); % NEW: synthesized raw gyro measurement [q_meas, r_meas]

k = 1; t = 0;
true_log(1,:) = [phi0, theta0, psi0, p.p0];
true_vel_log(1,:) = [u0, v0, w0];
true_alt_log(1) = 0;
true_qr_log(1,:) = [0, 0];
gyro_meas_log(1,:) = [p.bias_gyro_q, p.bias_gyro_r];

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
mag_log(k,:) = sensors.mag_body';
accel_log(k,:) = sensors.accel_body';
accel_true_log(k,:)= sensors.accel_body_true';
accel_NP_log(k,:) = sensors.accel_NP';
if sensors.gps_available
gps_vel_log(k,:) = sensors.gps_vel';
end

% --- NEW: true body-frame q_b, r_b (state_hifi indices 12, 13) ---
q_b_true = state_hifi(12);
r_b_true = state_hifi(13);
true_qr_log(k,:) = [q_b_true, r_b_true];

% --- NEW: synthesize raw gyro measurement = true body rate +
% assumed constant bias + noise. This is a MODELING ASSUMPTION
% (no real gyro hardware/datasheet) -- flagged same as Phase 2. ---
gyro_meas_log(k,:) = [q_b_true + p.bias_gyro_q + p.sigma_gyro_meas*randn(), ...
r_b_true + p.bias_gyro_r + p.sigma_gyro_meas*randn()];
end
time_log = time_log(1:k);
true_log = true_log(1:k,:);
true_vel_log = true_vel_log(1:k,:);
true_alt_log = true_alt_log(1:k);
mag_log = mag_log(1:k,:);
accel_log = accel_log(1:k,:);
accel_true_log = accel_true_log(1:k,:);
accel_NP_log = accel_NP_log(1:k,:);
gps_vel_log = gps_vel_log(1:k,:);
true_qr_log = true_qr_log(1:k,:);
gyro_meas_log = gyro_meas_log(1:k,:);
N2 = k;
fprintf('Truth model: %d steps, %.2f s\n', N2, t);

% ================================================
% STAGE C: MAGNETOMETER DE-BIASING (unchanged)
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
% STAGE D: ROLL EXTRACTION (unchanged)
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
% STAGE 2: pb measurement (zero-crossing) -- unchanged
% ================================================
idx_cal = find(mag_y(1:end-1) < 0 & mag_y(2:end) >= 0, 1, 'first');
phi_ref_offset = true_log(idx_cal,1);
[~, pb_est_stage2] = roll_spin_extractor(mag_y, time_log, phi_ref_offset);
pb_err_stage2 = (true_log(:,4) - pb_est_stage2) / (2*pi);
fprintf('Stage 2 pb: RMS = %.3f Hz\n', rms(pb_err_stage2(~isnan(pb_est_stage2))));

% ================================================
% STAGE E: MAIN ATTITUDE EKF -- 8-STATE (gyro-augmented)
% ================================================
% State: [dtheta_S; dpsi_S; q; r; pb; axial_bias; b_q; b_r]
x_E = [0; 0; 0; 0; p.p0; 0; 0; 0];

% P_E: old 6x6 block unchanged, new b_q/b_r initial variance ASSUMED
% (no datasheet) -- tune against validation, same caveat as Phase 2.
P_E = blkdiag(...
diag([deg2rad(2)^2, deg2rad(2)^2, deg2rad(1)^2, deg2rad(1)^2, ...
(2*pi*20)^2, 300^2]), ...
diag([deg2rad(1)^2, deg2rad(1)^2]));

% Q_E: old 6x6 block unchanged, new b_q/b_r process noise ASSUMED
% (slowly-varying bias, small random walk).
Q_E = blkdiag(...
diag([deg2rad(0.05)^2, deg2rad(0.05)^2, ...
deg2rad(0.001)^2, deg2rad(0.001)^2, ...
(2*pi*0.5)^2, 0.001^2]), ...
diag([deg2rad(0.001)^2, deg2rad(0.001)^2]));

R_pb = (2*pi*1.2)^2;

theta_est_E = zeros(N2,1);
psi_est_E = zeros(N2,1);
pb_est_E = zeros(N2,1);
axial_bias_est_E = zeros(N2,1);
q_est_E = zeros(N2,1);
r_est_E = zeros(N2,1);
b_q_est_E = zeros(N2,1); % NEW
b_r_est_E = zeros(N2,1); % NEW
dtheta_S_log = zeros(N2,1);
dpsi_S_log = zeros(N2,1);

theta_est_E(1) = theta0;
psi_est_E(1) = psi0;
pb_est_E(1) = p.p0;

P_log_E = zeros(8,8,N2); % NEW: now 8x8
P_log_E(:,:,1) = P_E;

for i = 2:N2
if ~any(isnan(gps_vel_log(i,:)))
vel_est = gps_vel_log(i,:)';
else
vel_est = true_vel_log(i,:)';
end
[rho_i,~] = atmosphere(max(0, true_alt_log(i)));
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

% --- NEW: raw gyro measurement for this step ---
gyro_meas_i = gyro_meas_log(i,:)';

% --- 8-state EKF step (Stage E v8, includes gyro update) ---
[x_E, P_E] = attitude_ekf_step_StageE(x_E, P_E, mag_log(i,:)', ...
phi_est_stageD(i), vel_est, p, qbar_i, V_i, p.dt, Q_E, R_axial_i, ...
pb_meas_i, R_pb, gyro_meas_i);

P_log_E(:,:,i) = P_E;
psi_v_i = atan2(vel_est(2), vel_est(1));
theta_v_i = atan2(-vel_est(3), sqrt(vel_est(1)^2+vel_est(2)^2));
theta_est_E(i) = theta_v_i + x_E(1);
psi_est_E(i) = psi_v_i + x_E(2);
pb_est_E(i) = x_E(5);
axial_bias_est_E(i) = x_E(6);
q_est_E(i) = x_E(3);
r_est_E(i) = x_E(4);
b_q_est_E(i) = x_E(7); % NEW
b_r_est_E(i) = x_E(8); % NEW
dtheta_S_log(i) = x_E(1);
dpsi_S_log(i) = x_E(2);
end
theta_err_E = rad2deg(wrapToPi(true_log(:,2) - theta_est_E));
psi_err_E = rad2deg(wrapToPi(true_log(:,3) - psi_est_E));
pb_err_E = (true_log(:,4) - pb_est_E) / (2*pi);

fprintf('\n=== Stage E (8-state, gyro-augmented) ===\n');
fprintf('RMS theta = %.3f deg, RMS psi = %.3f deg\n', rms(theta_err_E), rms(psi_err_E));
fprintf('RMS pb = %.3f Hz\n', rms(pb_err_E));
fprintf('Final axial bias = %.1f nT (true: %.0f)\n', axial_bias_est_E(end), p.bias_mag(1));

% --- NEW: gyro bias estimation diagnostics ---
fprintf('Final b_q estimate = %.4f deg/s (true assumed bias: %.4f deg/s)\n', ...
rad2deg(b_q_est_E(end)), rad2deg(p.bias_gyro_q));
fprintf('Final b_r estimate = %.4f deg/s (true assumed bias: %.4f deg/s)\n', ...
rad2deg(b_r_est_E(end)), rad2deg(p.bias_gyro_r));

% --- NEW: compare estimated q,r (S-frame) against true q,r (body-frame,
% rotated through ESTIMATED phi_roll for a fair like-for-like check) ---
q_r_true_Sframe = zeros(N2,2);
for i = 1:N2
phi_i = phi_est_stageD(i);
Rmat = [cos(phi_i), sin(phi_i); -sin(phi_i), cos(phi_i)]; % inverse of body-rotation
q_r_true_Sframe(i,:) = (Rmat * true_qr_log(i,:)')';
end
q_err_E = rad2deg(q_est_E - q_r_true_Sframe(:,1));
r_err_E = rad2deg(r_est_E - q_r_true_Sframe(:,2));
fprintf('RMS q (S-frame) error = %.3f deg/s, RMS r (S-frame) error = %.3f deg/s\n', ...
rms(q_err_E(2:end)), rms(r_err_E(2:end)));
fprintf('(Compare against pre-gyro-update baseline to see if the gyro channel actually helps.)\n');

% ================================================
% LAYER 6A: ACCELEROMETER NP ESTIMATOR (unchanged, consumes theta_est_E/
% psi_est_E/pb_est_E as before -- now from the 8-state filter)
% ================================================
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
for i = 2:N2
if ~any(isnan(gps_vel_log(i,:)))
vel_est = gps_vel_log(i,:)';
else
vel_est = true_vel_log(i,:)';
end
V_i = norm(vel_est);
if V_i < 1; V_i = 1; end
[rho_i,~] = atmosphere(max(0, true_alt_log(i)));
qbar_i = 0.5*rho_i*V_i^2;
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
end

% ================================================
% LAYER 6C (BASIC): velocity and position from clean accel (unchanged)
% ================================================
vel_est_accel = zeros(N2, 3);
pos_est_accel = zeros(N2, 3);
vel_est_accel(1,:) = [u0, v0, w0];
pos_est_accel(1,:) = [0, 0, 0];
g_inertial = [0; 0; 9.81];
for i = 2:N2
phi_i = phi_est_stageD(i);
theta_i = theta_est_E(i);
psi_i = psi_est_E(i);
Ry_i = [cos(theta_i),0,-sin(theta_i);0,1,0;sin(theta_i),0,cos(theta_i)];
Rz_i = [cos(psi_i),sin(psi_i),0;-sin(psi_i),cos(psi_i),0;0,0,1];
Rx_i = [1,0,0;0,cos(phi_i),-sin(phi_i);0,sin(phi_i),cos(phi_i)];
R_b2i = (Rx_i*Ry_i*Rz_i)';
f_inertial = R_b2i * accel_clean_log(i,:)';
a_inertial = f_inertial + g_inertial;
vel_est_accel(i,:) = vel_est_accel(i-1,:) + ...
0.5*p.dt*(a_inertial + (R_b2i*accel_clean_log(i-1,:)' + g_inertial))';
pos_est_accel(i,:) = pos_est_accel(i-1,:) + ...
0.5*p.dt*(vel_est_accel(i,:) + vel_est_accel(i-1,:));
end
vel_err = true_vel_log - vel_est_accel;
pos_err = [cumsum(true_vel_log(:,1))*p.dt, ...
cumsum(true_vel_log(:,2))*p.dt, ...
cumsum(true_vel_log(:,3))*p.dt] - pos_est_accel;
fprintf('\n=== Inertial integration ===\n');
fprintf('Velocity RMS error: [%.2f, %.2f, %.2f] m/s\n', ...
rms(vel_err(:,1)), rms(vel_err(:,2)), rms(vel_err(:,3)));
fprintf('Final position error: [%.1f, %.1f, %.1f] m\n', ...
pos_err(end,1), pos_err(end,2), pos_err(end,3));
NP_err = accel_NP_log - NP_est_log;
fprintf('\n=== Layer 6A: NP Estimator (physical 4-state) ===\n');
fprintf('NP estimation RMS error (after t=15s, balloting gone):\n');
idx_settled = time_log > 3*p.balloting_tau;
fprintf(' Axial (X): %.4f m/s^2 (true NP = %.3f m/s^2)\n', ...
rms(NP_err(idx_settled,1)), rms(accel_NP_log(idx_settled,1)));
fprintf(' Radial (Y): %.4f m/s^2 (true NP RMS = %.1f m/s^2)\n', ...
rms(NP_err(idx_settled,2)), rms(accel_NP_log(idx_settled,2)));
fprintf(' Radial (Z): %.4f m/s^2 (true NP RMS = %.1f m/s^2)\n', ...
rms(NP_err(idx_settled,3)), rms(accel_NP_log(idx_settled,3)));
fprintf('\nNOTE: P_log_E (8x8xN2 Stage E covariance history) is now available in the workspace.\n');
fprintf('NOTE: Run this WITH and WITHOUT the gyro update (e.g. by comparing against\n');
fprintf('the previous 6-state script''s dpsi_S/psi_est_E RMS) to see whether the\n');
fprintf('gyro channel genuinely improves the known yaw-observability weakness\n');
fprintf('(finding C, prior handoff doc) or not -- that comparison is the actual\n');
fprintf('point of this whole exercise and has NOT been run yet.\n');
