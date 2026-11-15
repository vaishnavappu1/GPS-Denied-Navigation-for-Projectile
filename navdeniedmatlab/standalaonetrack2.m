% STANDALONE TEST: Track 2 inertial predictor vs truth
% Run AFTER main_layer7_v3.m (reuses workspace: true_vel_log, time_log,
% accel_clean_log, phi_est_stageD, theta_est_E, psi_est_E, p, N2)
%
% Verifies track2_inertial_step.m (refactored to add Phi_a/P_a) gives
% the SAME standalone results as the original inline Track 2 integration
% (velocity RMS [40.88, 9.06, 46.03] m/s, position error
% [-3183.2, -366.1, 3771.6] m, from the last confirmed main script run).
% If these numbers match, the refactor introduced no regression and
% it's safe to proceed to building the combiner.

fprintf('=== TRACK 2 STANDALONE TEST (covariance-propagation version) ===\n\n');

elev_rad = deg2rad(p.elev);
az_rad   = deg2rad(p.az);
u0 = p.V0*cos(elev_rad)*cos(az_rad);
v0 = p.V0*cos(elev_rad)*sin(az_rad);
w0 = -p.V0*sin(elev_rad);

pos_a = [0; 0; 0];
vel_a = [u0; v0; w0];

% Same initial covariance as Track 1 standalone test, for consistency
sigma_pos0 = 1.0;
sigma_vel0 = 1.0;
P_combined_a = diag([sigma_pos0^2, sigma_pos0^2, sigma_pos0^2, ...
    sigma_vel0^2, sigma_vel0^2, sigma_vel0^2]);

% Process noise: inertial track drifts more than dynamic track,
% reflecting residual NP error + sensor noise after integration.
% sigma_accel_err ~ 0.5 m/s^2 effective (NP residual ~0.4 m/s^2 radial,
% 0.05 m/s^2 axial, per Layer 6A validated results)
sigma_a_inertial = 0.5;  % m/s^2
Q_a = diag([0, 0, 0, ...
    (sigma_a_inertial*p.dt)^2, ...
    (sigma_a_inertial*p.dt)^2, ...
    (sigma_a_inertial*p.dt)^2]);

pos_a_log = zeros(N2, 3);
vel_a_log = zeros(N2, 3);
pos_a_log(1,:) = pos_a';
vel_a_log(1,:) = vel_a';

for i = 2:N2
    [pos_a, vel_a, Phi_a, P_combined_a] = track2_inertial_step(...
        pos_a, vel_a, P_combined_a, ...
        accel_clean_log(i-1,:)', accel_clean_log(i,:)', ...
        phi_est_stageD(i-1), theta_est_E(i-1), psi_est_E(i-1), ...
        phi_est_stageD(i),   theta_est_E(i),   psi_est_E(i), ...
        p.dt, p, Q_a);

    pos_a_log(i,:) = pos_a';
    vel_a_log(i,:) = vel_a';
end

% True position from integrating true velocity
true_pos = zeros(N2,3);
for i = 2:N2
    true_pos(i,:) = true_pos(i-1,:) + ...
        0.5*p.dt*(true_vel_log(i,:) + true_vel_log(i-1,:));
end

vel_err_a = true_vel_log - vel_a_log;
pos_err_a = true_pos - pos_a_log;

fprintf('Track 2 (inertial, covariance-propagation version):\n');
fprintf('  Velocity RMS error: [%.2f %.2f %.2f] m/s\n', ...
    rms(vel_err_a(:,1)), rms(vel_err_a(:,2)), rms(vel_err_a(:,3)));
fprintf('  Final position error: [%.1f %.1f %.1f] m\n', ...
    pos_err_a(end,1), pos_err_a(end,2), pos_err_a(end,3));
fprintf('  Final velocity error: [%.2f %.2f %.2f] m/s\n', ...
    vel_err_a(end,1), vel_err_a(end,2), vel_err_a(end,3));

fprintf('\nFor comparison, original (main script) Track 2 result was:\n');
fprintf('  Velocity RMS error: [40.88, 9.06, 46.03] m/s\n');
fprintf('  Final position error: [-3183.2, -366.1, 3771.6] m\n');
fprintf('\nIf the numbers above closely match these, the refactor is\n');
fprintf('verified and it is safe to proceed to building the combiner.\n');




track2_final_pos   = pos_a_log(end, 1:3); 
fprintf('Track 2:   [%8.2f, %8.2f, %8.2f]\n', track2_final_pos(1), track2_final_pos(2), track2_final_pos(3));

