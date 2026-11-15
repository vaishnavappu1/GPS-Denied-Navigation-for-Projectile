% =========================================================================
% LAYER 6C: CASCADED FUSION WITH TIME-RESOLVED DIAGNOSTICS + P_da CAP
% Adds: (1) full time-history logging of L norm / denominator min-eig, to
% see EXACTLY when instability occurs, and (2) a hard eigenvalue cap on
% P_da every step -- a structural fix that doesn't rely on GPS resets to
% bail out an otherwise-unbounded cross-covariance.
% =========================================================================
fprintf('=== CASCADED FUSION WITH DIAGNOSTICS + P_da EIGENVALUE CAP ===\n\n');

elev_rad = deg2rad(p.elev); az_rad = deg2rad(p.az);
u0 = p.V0*cos(elev_rad)*cos(az_rad); v0 = p.V0*cos(elev_rad)*sin(az_rad); w0 = -p.V0*sin(elev_rad);

pos_d_state = [0; 0; 0]; vel_d_state = [u0; v0; w0];
pos_a_state = [0; 0; 0]; vel_a_state = [u0; v0; w0];

P_initial = diag([1.0^2, 1.0^2, 1.0^2, 1.0^2, 1.0^2, 1.0^2]);
% --- FIX: avoid identical P_d/P_a/P_da at t=0 (same near-singular
% denominator issue as GPS resets -- this was the actual source of the
% worst L-norm spike, 413.6, occurring at t=0.00s even after the GPS
% reset fix was applied, since initialization is a separate code path). ---
P_d_state  = P_initial + 1e-3*eye(6);
P_a_state  = P_initial + 2e-3*eye(6);
P_da_state = P_initial * 0.98;

sigma_a_model = 0.1; sigma_pos_model = 0.01;
Q_d = diag([sigma_pos_model^2, sigma_pos_model^2, sigma_pos_model^2, ...
            (sigma_a_model*p.dt)^2, (sigma_a_model*p.dt)^2, (sigma_a_model*p.dt)^2]);
sigma_a_inertial = 0.5; sigma_pos_inertial = 0.01;
Q_a = diag([sigma_pos_inertial^2, sigma_pos_inertial^2, sigma_pos_inertial^2, ...
            (sigma_a_inertial*p.dt)^2, (sigma_a_inertial*p.dt)^2, (sigma_a_inertial*p.dt)^2]);
Q_da = zeros(6,6);

% --- NEW: P_da eigenvalue cap ---
% Set relative to typical P_d/P_a scale so it can never dominate the
% denominator regardless of how long P_da propagates unconstrained.
P_DA_EIG_CAP = 1e4;  % tune based on observed P_d/P_a scales; start conservative

% --- NEW: continuous per-step soft deflation parameters (GPS-denied tail fix) ---
% Equivalent-decay calibration: previous periodic version used factor=0.85
% every 5.0s. To get similar (or slightly stronger, since it's smoother)
% cumulative decay per-step: rate^(steps_per_interval) = 0.85
% steps_per_interval = 5.0/p.dt = 5000 (at dt=0.001)
% rate = 0.85^(1/5000) ~ 0.99996749 -- start here, then tune stronger/weaker
DEFLATE_RATE_PER_STEP = 0.85^(1/(5.0/p.dt));

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

pos_c_log(1,:) = pos_d_state'; vel_c_log(1,:) = vel_d_state';

for i = 2:N2
    t = i * p.dt;
    [rho_i,~] = atmosphere(max(0, true_alt_log(i)));

    [pos_d, vel_d, Phi_d, P_d] = track1_dynamic_step(...
        pos_d_state, vel_d_state, P_d_state, ...
        dtheta_S_log(i), dpsi_S_log(i), ...
        pb_est_E(i), q_est_E(i), r_est_E(i), ...
        p.B_earth, rho_i, p.dt, p, Q_d);
    x_d = [pos_d; vel_d];

    [pos_a, vel_a, Phi_a, P_a] = track2_inertial_step(...
        pos_a_state, vel_a_state, P_a_state, ...
        accel_clean_log(i-1,:)', accel_clean_log(i,:)', ...
        phi_est_stageD(i-1), theta_est_E(i-1), psi_est_E(i-1), ...
        phi_est_stageD(i),   theta_est_E(i),   psi_est_E(i), ...
        p.dt, p, Q_a);
    x_a = [pos_a; vel_a];

    P_da = Phi_d * P_da_state * Phi_a' + Q_da;

    % --- Deflation REMOVED after analysis: the growth of P_da (and
    % therefore L) in the GPS-denied tail was found to be smooth,
    % monotonic, and physically legitimate -- it reflects genuine growing
    % correlation between Track1/Track2 errors (both derive from shared
    % Stage E attitude estimates, so their errors correlate more over
    % time without correction). Artificially deflating P_da discards real
    % information and measurably hurt Z-axis accuracy (20.38->20.60) for
    % only a cosmetic reduction in max L norm (31.3->26.5), which was not
    % actually indicating a numerical problem (denominator stayed healthy,
    % min-eigenvalue ~3.1e-3, nowhere near singular). Reverted.

    % --- hard eigenvalue cap on P_da every step (safety net, kept as
    % insurance against genuine runaway growth in other scenarios) ---
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

        pos_d_state = x_c_corrected(1:3); vel_d_state = x_c_corrected(4:6);
        pos_a_state = x_c_corrected(1:3); vel_a_state = x_c_corrected(4:6);

        % --- FIX: avoid setting P_d, P_a, P_da to IDENTICAL matrices.
        % Doing so makes (P_a + P_d - P_ad - P_da) near-exactly singular
        % right after every reset (confirmed via diagnostic: denominator
        % min-eigenvalue crashes to ~1e-7 at every single GPS fix).
        % Give each track a small, DISTINCT perturbation so the denominator
        % never starts from an exactly-singular configuration, and make
        % P_da_state slightly smaller than P_d/P_a (representing imperfect,
        % not-yet-fully-correlated tracks immediately post-reset). ---
        eps_d = 1e-3;
        eps_a = 2e-3;  % deliberately different from eps_d
        P_d_state  = P_c_corrected + eps_d*eye(6);
        P_a_state  = P_c_corrected + eps_a*eye(6);
        P_da_state = P_c_corrected * 0.98;  % slightly less than full correlation

        pos_c_log(i,:) = x_c_corrected(1:3)';
        vel_c_log(i,:) = x_c_corrected(4:6)';
    else
        pos_d_state = pos_d; vel_d_state = vel_d; P_d_state = P_d;
        pos_a_state = pos_a; vel_a_state = vel_a; P_a_state = P_a;
        P_da_state  = P_da;

        pos_c_log(i,:) = x_c_minus(1:3)';
        vel_c_log(i,:) = x_c_minus(4:6)';
    end
end

% ================================================
% RESULTS
% ================================================
vel_err_c = true_vel_log - vel_c_log;
pos_err_c = true_pos - pos_c_log;

fprintf('=========================================================\n');
fprintf('FINAL FUSION RESULTS (WITH P_da CAP)\n');
fprintf('=========================================================\n');
fprintf('CASCADED COMBINED Velocity RMS: [%.2f, %.2f, %.2f] m/s\n', ...
    rms(vel_err_c(:,1)), rms(vel_err_c(:,2)), rms(vel_err_c(:,3)));
fprintf('CASCADED COMBINED Final position error: [%.1f, %.1f, %.1f] m\n\n', ...
    pos_err_c(end,1), pos_err_c(end,2), pos_err_c(end,3));

fprintf('L norm: min=%.4f max=%.4f\n', min(L_norm_log(2:end)), max(L_norm_log(2:end)));
fprintf('Denominator min-eigenvalue: min=%.4e max=%.4e\n', ...
    min(denom_mineig_log(2:end)), max(denom_mineig_log(2:end)));
fprintf('P_da cap triggered on %.1f%% of steps\n\n', 100*mean(P_da_cap_triggered_log(2:end)));

% --- Find WHEN the worst L values occur ---
[worst_L, worst_idx] = max(L_norm_log(2:end));
fprintf('Worst L norm (%.2f) occurs at t=%.2fs (GPS active at that time: %d)\n\n', ...
    worst_L, time_log(worst_idx+1), gps_active_log(worst_idx+1));

% --- Time-resolved plots ---
figure('Name','Time-Resolved Fusion Diagnostics','Position',[100 100 1200 900]);
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

fprintf('--- INTERPRETATION ---\n');
fprintf('Check WHERE the worst L norm occurs (printed above). If it''s after\n');
fprintf('t=40s (GPS cutoff): confirms the instability resumes once resets stop,\n');
fprintf('exactly as predicted -- the cap should now be preventing this. If L\n');
fprintf('norm stays bounded (<~5-10) throughout, especially post-40s, the cap\n');
fprintf('is working as a real structural fix, not just delaying the problem.\n');