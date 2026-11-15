% =========================================================================
% TRACK 1 STANDALONE TEST -- PRODUCTION REALISTIC (using Stage E's REAL
% estimates, not ground-truth override)
%
% This is the honest end-to-end test: Track 1 driven by whatever Stage E
% actually estimated (dtheta_S_log, dpsi_S_log, pb_est_E, q_est_E, r_est_E),
% exactly as the real deployed system would have available. All prior
% standalone tests tonight used ground-truth angles to isolate bugs in
% Track1's own math -- this is the first test of the REAL, full pipeline.
%
% Expect this to perform WORSE than the truth-angle tests, especially in
% Z near apogee, since Stage E's own dtheta_S estimate was shown to
% degrade near apogee (correlation -0.42 vs true alpha) -- this run
% combines that estimation error with whatever residual Track1 error
% remains, rather than isolating them.
% =========================================================================
fprintf('=== TRACK 1 STANDALONE TEST (PRODUCTION REALISTIC) ===\n\n');

elev_rad = deg2rad(p.elev);
az_rad   = deg2rad(p.az);
u0 = p.V0*cos(elev_rad)*cos(az_rad);
v0 = p.V0*cos(elev_rad)*sin(az_rad);
w0 = -p.V0*sin(elev_rad);

pos_prev = [0; 0; 0];
vel_prev = [u0; v0; w0];

pos_d_log = zeros(N2, 3);
vel_d_log = zeros(N2, 3);
pos_d_log(1,:) = pos_prev';
vel_d_log(1,:) = vel_prev';
P_prev_combined = eye(6);
Q_d             = zeros(6,6);

for i = 2:N2
    [rho_i,~] = atmosphere(max(0, true_alt_log(i)));

    % --- USE STAGE E'S REAL ESTIMATES (production-realistic) ---
    dtheta_S_real = dtheta_S_log(i-1);
    dpsi_S_real   = dpsi_S_log(i-1);
    pb_real       = pb_est_E(i-1);
    q_real        = q_est_E(i-1);
    r_real        = r_est_E(i-1);

    [pos_d, vel_d, Phi_d, P_d] = track1_dynamic_step(...
        pos_prev, vel_prev, P_prev_combined, ...
        dtheta_S_real, dpsi_S_real, ...
        pb_real, q_real, r_real, ...
        p.B_earth, ...
        rho_i, p.dt, p, Q_d);

    pos_d_log(i,:) = pos_d';
    vel_d_log(i,:) = vel_d';

    pos_prev = pos_d;
    vel_prev = vel_d;
    P_prev_combined = P_d;
end

true_pos = zeros(N2,3);
for i = 2:N2
    true_pos(i,:) = true_pos(i-1,:) + ...
        0.5*p.dt*(true_vel_log(i,:) + true_vel_log(i-1,:));
end

vel_err_d = true_vel_log - vel_d_log;
pos_err_d = true_pos - pos_d_log;

fprintf('Track 1 (PRODUCTION REALISTIC, using Stage E real estimates):\n');
fprintf('  Velocity RMS error: [%.2f %.2f %.2f] m/s\n', ...
    rms(vel_err_d(:,1)), rms(vel_err_d(:,2)), rms(vel_err_d(:,3)));
fprintf('  Final position error: [%.1f %.1f %.1f] m\n', ...
    pos_err_d(end,1), pos_err_d(end,2), pos_err_d(end,3));
fprintf('  Final velocity error: [%.2f %.2f %.2f] m/s\n\n', ...
    vel_err_d(end,1), vel_err_d(end,2), vel_err_d(end,3));

fprintf('For comparison, ground-truth-angle baseline was: [18.39, 10.87, 35.25] m/s\n\n');

% ================================================
% PLOTS
% ================================================
figure('Name','Track1 Production-Realistic: Velocity Tracking','Position',[50 50 1400 800]);

subplot(3,2,1);
plot(time_log, true_vel_log(:,1), 'b-', time_log, vel_d_log(:,1), 'r--', 'LineWidth', 1);
legend('true u','est u'); grid on; title('X-velocity (u)'); xlabel('t [s]'); ylabel('m/s');

subplot(3,2,2);
plot(time_log, vel_err_d(:,1), 'k-'); grid on; title('X-velocity error'); xlabel('t [s]'); ylabel('m/s');

subplot(3,2,3);
plot(time_log, true_vel_log(:,2), 'b-', time_log, vel_d_log(:,2), 'r--', 'LineWidth', 1);
legend('true v','est v'); grid on; title('Y-velocity (v)'); xlabel('t [s]'); ylabel('m/s');

subplot(3,2,4);
plot(time_log, vel_err_d(:,2), 'k-'); grid on; title('Y-velocity error'); xlabel('t [s]'); ylabel('m/s');

subplot(3,2,5);
plot(time_log, true_vel_log(:,3), 'b-', time_log, vel_d_log(:,3), 'r--', 'LineWidth', 1);
legend('true w','est w'); grid on; title('Z-velocity (w)'); xlabel('t [s]'); ylabel('m/s');

subplot(3,2,6);
plot(time_log, vel_err_d(:,3), 'k-'); grid on; title('Z-velocity error'); xlabel('t [s]'); ylabel('m/s');

figure('Name','Track1 Production-Realistic: Position Tracking','Position',[100 100 1400 500]);
subplot(1,3,1);
plot(time_log, true_pos(:,1), 'b-', time_log, pos_d_log(:,1), 'r--'); grid on;
legend('true X','est X'); title('X position'); xlabel('t [s]'); ylabel('m');

subplot(1,3,2);
plot(time_log, true_pos(:,2), 'b-', time_log, pos_d_log(:,2), 'r--'); grid on;
legend('true Y','est Y'); title('Y position'); xlabel('t [s]'); ylabel('m');

subplot(1,3,3);
plot(time_log, true_pos(:,3), 'b-', time_log, pos_d_log(:,3), 'r--'); grid on;
legend('true Z','est Z'); title('Z position'); xlabel('t [s]'); ylabel('m');

% Overlay velocity magnitude to visually check apogee correlation
figure('Name','Z-error vs Velocity (Apogee Check)','Position',[150 150 900 500]);
subplot(2,1,1);
plot(time_log, sqrt(sum(true_vel_log.^2,2)), 'k-'); grid on;
title('True Velocity Magnitude'); ylabel('m/s');
subplot(2,1,2);
plot(time_log, abs(vel_err_d(:,3)), 'm-'); grid on;
title('|Z-velocity error| (production realistic)'); xlabel('t [s]'); ylabel('m/s');