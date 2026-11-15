% =========================================================================
% DIAGNOSTIC: b_r divergence investigation
% Run this AFTER main_layer7_6A_8state_gyro.m has populated the workspace
% (needs: b_r_est_E, b_q_est_E, r_est_E, q_est_E, phi_est_stageD,
% true_log, time_log, dpsi_S_log, dtheta_S_log, p, true_qr_log)
% =========================================================================
fprintf('=== b_r DIVERGENCE DIAGNOSTIC ===\n\n');

phi_err_stageD_i = rad2deg(wrapToPi(true_log(:,1) - phi_est_stageD));

% --- HYPOTHESIS 1: r / b_r non-identifiability (compensating drift) ---
fprintf('--- H1: r_est_E vs b_r_est_E correlation ---\n');
valid_idx = 2:length(b_r_est_E);
cc = corrcoef(r_est_E(valid_idx), b_r_est_E(valid_idx));
corr_r_br = cc(1,2);
fprintf('corr(r_est_E, b_r_est_E) = %.3f (strong negative -> supports H1, non-identifiability)\n\n', corr_r_br);

% --- HYPOTHESIS 2: Stage D roll error correlates with b_r drift rate ---
fprintf('--- H2: |phi_roll error| vs |d(b_r)/dt| correlation ---\n');
db_r_dt = diff(b_r_est_E) ./ diff(time_log);
db_r_dt = [0; db_r_dt];
cc = corrcoef(abs(phi_err_stageD_i(valid_idx)), abs(db_r_dt(valid_idx)));
corr_phi_dbr = cc(1,2);
fprintf('corr(|phi_roll err|, |d(b_r)/dt|) = %.3f (strong positive -> supports H2, roll-leakage)\n\n', corr_phi_dbr);

% --- HYPOTHESIS 3: b_r divergence timing vs known dpsi_S weakness ---
% (finding C: dpsi_S estimate severely underperforms -- check whether
% b_r's divergence timeline matches dpsi_S estimation degrading)
fprintf('--- H3: dpsi_S_log vs b_r_est_E correlation ---\n');
cc = corrcoef(dpsi_S_log(valid_idx), b_r_est_E(valid_idx));
corr_dpsi_br = cc(1,2);
fprintf('corr(dpsi_S_log, b_r_est_E) = %.3f (strong correlation -> supports H3, yaw-weakness leakage)\n\n', corr_dpsi_br);

% --- Monotonic drift vs oscillation check (distinguishes divergence type) ---
fprintf('--- Drift character: monotonic (integrator-like) vs oscillatory (noisy but bounded) ---\n');
sign_changes = sum(diff(sign(db_r_dt(valid_idx))) ~= 0);
fprintf('Number of sign changes in d(b_r)/dt: %d (out of %d steps)\n', sign_changes, length(valid_idx));
fprintf(' Few sign changes -> monotonic drift (integrator-like, e.g. H1 or H2 accumulating)\n');
fprintf(' Many sign changes -> oscillatory/noisy but not necessarily diverging (less concerning)\n\n');

% --- Time-resolved plots ---
figure('Name','b_r Divergence Diagnostic','Position',[100 100 1300 900]);

subplot(3,2,1);
plot(time_log, rad2deg(b_q_est_E), 'b-'); hold on;
plot(time_log, rad2deg(b_r_est_E), 'r-');
yline(rad2deg(p.bias_gyro_q), 'b--');
yline(rad2deg(p.bias_gyro_r), 'r--');
grid on; legend('b_q est','b_r est','b_q true','b_r true','Location','best');
title('Gyro bias estimates vs assumed truth'); ylabel('deg/s'); xlabel('t [s]');

subplot(3,2,2);
plot(time_log, rad2deg(q_est_E), 'b-'); hold on;
plot(time_log, rad2deg(r_est_E), 'r-');
grid on; legend('q est (S-frame)','r est (S-frame)','Location','best');
title('S-frame rate estimates'); ylabel('deg/s'); xlabel('t [s]');

subplot(3,2,3);
plot(time_log, phi_err_stageD_i, 'k-');
grid on; title('Stage D roll error over time'); ylabel('deg'); xlabel('t [s]');

subplot(3,2,4);
plot(time_log, rad2deg(dtheta_S_log), 'b-'); hold on;
plot(time_log, rad2deg(dpsi_S_log), 'r-');
grid on; legend('dtheta_S','dpsi_S','Location','best');
title('Attitude deviation states'); ylabel('deg'); xlabel('t [s]');

subplot(3,2,5);
scatter(r_est_E(valid_idx), b_r_est_E(valid_idx), 5, time_log(valid_idx), 'filled');
colorbar; grid on;
xlabel('r_est_E [rad/s]'); ylabel('b_r_est_E [rad/s]');
title(sprintf('r vs b_r phase plot (color=time), corr=%.2f', corr_r_br));

subplot(3,2,6);
scatter(abs(phi_err_stageD_i(valid_idx)), abs(rad2deg(db_r_dt(valid_idx))), 5, time_log(valid_idx), 'filled');
colorbar; grid on;
xlabel('|phi_{roll} error| [deg]'); ylabel('|d(b_r)/dt| [deg/s per s]');
title(sprintf('Roll error vs b_r drift rate (color=time), corr=%.2f', corr_phi_dbr));

fprintf('=== INTERPRETATION GUIDE ===\n');
fprintf('- If H1 correlation is strongly negative (e.g. < -0.6): r/b_r are trading off\n');
fprintf(' against each other -- a genuine non-identifiability problem needing either\n');
fprintf(' better Q tuning (constrain r more via process model) or accepting b_r is\n');
fprintf(' not separably recoverable from this measurement alone.\n');
fprintf('- If H2 correlation is strongly positive: Stage D''s ~2.2 deg roll error is\n');
fprintf(' leaking into the gyro rotation and getting absorbed into b_r -- would need\n');
fprintf(' either a better Stage D or accepting this as a coupled-error source.\n');
fprintf('- If H3 correlation is strong: b_r is absorbing the KNOWN dpsi_S weakness\n');
fprintf(' (finding C) rather than being a clean, independent bias estimate.\n');
fprintf('- Check which hypothesis has the HIGHEST |correlation| -- that''s your\n');
fprintf(' primary suspect. They are not mutually exclusive; more than one may\n');
fprintf(' contribute.\n');