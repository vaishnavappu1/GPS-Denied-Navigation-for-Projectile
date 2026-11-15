% =========================================================================
% LAYER 6C: POSITION/VELOCITY ESTIMATOR -- CORRECTED (DECOUPLED)
% Track 1 (dynamic) + Track 2 (inertial) + Combiner (patent eq.41-43)
% Run AFTER main_layer7_v3.m (reuses workspace variables)
%
% FIX: Track1 and Track2 now propagate from THEIR OWN independent state
% history, not re-seeded from the combined estimate x_c every step. This
% matches the patent's explicit requirement (col. 19): "The integrate
% accelerometer model 650... and the Integrate dynamic model 640... both
% produce INDEPENDENT predictions of the position and velocity states."
% The previous version re-seeded pos_prev/vel_prev = x_c(1:3)/x_c(4:6)
% every step, which is the SAME feedback-contamination bug identified
% at the very start of this debugging session (destructive feedback loop
% where the combined state pollutes each track's own propagation).
% =========================================================================
fprintf('=== LAYER 6C: FULL COMBINED PIPELINE (CORRECTED, DECOUPLED) ===\n\n');

% ================================================
% INITIALIZATION
% ================================================
elev_rad = deg2rad(p.elev);
az_rad   = deg2rad(p.az);

u0 = p.V0*cos(elev_rad)*cos(az_rad);
v0 = p.V0*cos(elev_rad)*sin(az_rad);
w0 = -p.V0*sin(elev_rad);

% --- Combined output state (REPORTING ONLY, never fed back) ---
pos_c = [0; 0; 0];
vel_c = [u0; v0; w0];

% --- NEW: independent state histories for each track ---
pos_d_state = [0; 0; 0];
vel_d_state = [u0; v0; w0];
pos_a_state = [0; 0; 0];
vel_a_state = [u0; v0; w0];

sigma_pos0 = 1.0;
sigma_vel0 = 1.0;
P_combined = diag([sigma_pos0^2, sigma_pos0^2, sigma_pos0^2, ...
                    sigma_vel0^2, sigma_vel0^2, sigma_vel0^2]);

% Each track also needs its OWN covariance now that they're decoupled
P_d_state = P_combined;
P_a_state = P_combined;

% ================================================
% PROCESS NOISE
% ================================================
sigma_a_model    = 0.1;
sigma_pos_model  = 0.01;
Q_d = diag([sigma_pos_model^2, sigma_pos_model^2, sigma_pos_model^2, ...
            (sigma_a_model*p.dt)^2, (sigma_a_model*p.dt)^2, (sigma_a_model*p.dt)^2]);

sigma_a_inertial   = 0.5;
sigma_pos_inertial = 0.01;
Q_a = diag([sigma_pos_inertial^2, sigma_pos_inertial^2, sigma_pos_inertial^2, ...
            (sigma_a_inertial*p.dt)^2, (sigma_a_inertial*p.dt)^2, (sigma_a_inertial*p.dt)^2]);

% ================================================
% PREALLOCATED LOGS
% ================================================
pos_c_log = zeros(N2, 3);
vel_c_log = zeros(N2, 3);
pos_d_log = zeros(N2, 3);   % now the ACTUAL track feeding the combiner (not a duplicate diagnostic)
pos_a_log = zeros(N2, 3);
vel_d_log = zeros(N2, 3);
vel_a_log = zeros(N2, 3);

L_max = zeros(N2,1);
L_min = zeros(N2,1);
L_norm = zeros(N2,1);
P_eig_min = zeros(N2,1);
P_eig_max = zeros(N2,1);

pos_c_log(1,:) = pos_c';
vel_c_log(1,:) = vel_c';
pos_d_log(1,:) = pos_d_state';
pos_a_log(1,:) = pos_a_state';
vel_d_log(1,:) = vel_d_state';
vel_a_log(1,:) = vel_a_state';

print_interval = 5.0 / p.dt;

% ================================================
% MAIN LOOP
% ================================================
for i = 2:N2
    [rho_i,~] = atmosphere(max(0, true_alt_log(i)));

    % --- TRACK 1: propagates from ITS OWN previous state ---
    [pos_d, vel_d, Phi_d, P_d] = track1_dynamic_step(...
        pos_d_state, vel_d_state, P_d_state, ...
        dtheta_S_log(i), dpsi_S_log(i), ...
        pb_est_E(i), q_est_E(i), r_est_E(i), ...
        p.B_earth, rho_i, p.dt, p, Q_d);
    x_d = [pos_d; vel_d];

    % --- TRACK 2: propagates from ITS OWN previous state ---
    [pos_a, vel_a, Phi_a, P_a] = track2_inertial_step(...
        pos_a_state, vel_a_state, P_a_state, ...
        accel_clean_log(i-1,:)', accel_clean_log(i,:)', ...
        phi_est_stageD(i-1), theta_est_E(i-1), psi_est_E(i-1), ...
        phi_est_stageD(i),   theta_est_E(i),   psi_est_E(i), ...
        p.dt, p, Q_a);
    x_a = [pos_a; vel_a];

    if mod(i, print_interval) == 0
        true_accel   = (true_vel_log(i,:) - true_vel_log(i-1,:)) / p.dt;
        accel_d_calc = (vel_d - vel_d_state) / p.dt;
        accel_a_calc = (vel_a - vel_a_state) / p.dt;
        fprintf('--- Snapshot at t=%.1fs (Step %d) ---\n', i*p.dt, i);
        fprintf('  Truth Accel:    [%8.2f, %8.2f, %8.2f]\n', true_accel(1), true_accel(2), true_accel(3));
        fprintf('  Track 1 (Aero): [%8.2f, %8.2f, %8.2f]\n', accel_d_calc(1), accel_d_calc(2), accel_d_calc(3));
        fprintf('  Track 2 (Iner): [%8.2f, %8.2f, %8.2f]\n', accel_a_calc(1), accel_a_calc(2), accel_a_calc(3));
        fprintf('  Track1 Pos:     [%8.1f, %8.1f, %8.1f]\n', pos_d(1), pos_d(2), pos_d(3));
        fprintf('  Track2 Pos:     [%8.1f, %8.1f, %8.1f]\n', pos_a(1), pos_a(2), pos_a(3));
    end

    % --- COMBINER (blends for REPORTING ONLY -- does not feed back) ---
    [x_c, P_combined, L] = layer6C_combiner_step(...
        x_d, P_d, Phi_d, x_a, P_a, Phi_a, P_combined);

    L_max(i)  = max(L(:));
    L_min(i)  = min(L(:));
    L_norm(i) = norm(L, 'fro');

    eig_P = eig(P_combined);
    P_eig_min(i) = min(eig_P);
    P_eig_max(i) = max(eig_P);

    % --- UPDATE each track's OWN state history for next iteration ---
    pos_d_state = pos_d;
    vel_d_state = vel_d;
    P_d_state   = P_d;

    pos_a_state = pos_a;
    vel_a_state = vel_a;
    P_a_state   = P_a;

    % --- Logging ---
    pos_c_log(i,:) = x_c(1:3)';
    vel_c_log(i,:) = x_c(4:6)';
    pos_d_log(i,:) = pos_d';
    pos_a_log(i,:) = pos_a';
    vel_d_log(i,:) = vel_d';
    vel_a_log(i,:) = vel_a';
end

% ================================================
% TRUE POSITION REFERENCE
% ================================================
true_pos = zeros(N2,3);
for i = 2:N2
    true_pos(i,:) = true_pos(i-1,:) + ...
        0.5*p.dt*(true_vel_log(i,:) + true_vel_log(i-1,:));
end

% ================================================
% RESULTS
% ================================================
vel_err_c = true_vel_log - vel_c_log;
pos_err_c = true_pos - pos_c_log;
pos_err_d = true_pos - pos_d_log;
pos_err_a = true_pos - pos_a_log;
vel_err_d = true_vel_log - vel_d_log;
vel_err_a = true_vel_log - vel_a_log;

fprintf('\nCOMBINED ESTIMATE (Layer 6C, patent eq.41-43, DECOUPLED tracks):\n');
fprintf('  Velocity RMS error: [%.2f %.2f %.2f] m/s\n', ...
    rms(vel_err_c(:,1)), rms(vel_err_c(:,2)), rms(vel_err_c(:,3)));
fprintf('  Final position error: [%.1f %.1f %.1f] m\n', ...
    pos_err_c(end,1), pos_err_c(end,2), pos_err_c(end,3));

fprintf('\nTrack 1 (INDEPENDENT, decoupled) Velocity RMS: [%.2f %.2f %.2f] m/s\n', ...
    rms(vel_err_d(:,1)), rms(vel_err_d(:,2)), rms(vel_err_d(:,3)));
fprintf('Track 1 Final position error: [%.1f %.1f %.1f] m\n', ...
    pos_err_d(end,1), pos_err_d(end,2), pos_err_d(end,3));

fprintf('\nTrack 2 (INDEPENDENT, decoupled) Velocity RMS: [%.2f %.2f %.2f] m/s\n', ...
    rms(vel_err_a(:,1)), rms(vel_err_a(:,2)), rms(vel_err_a(:,3)));
fprintf('Track 2 Final position error: [%.1f %.1f %.1f] m\n', ...
    pos_err_a(end,1), pos_err_a(end,2), pos_err_a(end,3));

fprintf('\n--- Combiner gain diagnostics ---\n');
fprintf('L norm: min=%.4f max=%.4f\n', min(L_norm(2:end)), max(L_norm(2:end)));
fprintf('P_combined eigenvalues: min=%.4e max=%.4e\n', min(P_eig_min(2:end)), max(P_eig_max(2:end)));

fprintf('\n--- WHAT TO CHECK ---\n');
fprintf('Compare Track1/Track2 RMS here against your earlier standalone tests --\n');
fprintf('they should now be genuinely independent (no shared-seed contamination).\n');
fprintf('Then check if the COMBINED estimate RMS is better than EITHER track alone\n');
fprintf('on each axis -- that is the real test of whether fusion is working.\n');