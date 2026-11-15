% =========================================================================
% PHASE 1b: CORRECTED GYRO MEASUREMENT MODEL (body-frame rotation fix)
%
% SUPERSEDES Phase 1's gyro_h_StageE / gyro_jacobian_numeric.
%
% ROOT ISSUE FOUND (design review, not yet a runtime bug since Phase 1's
% version was never wired into the EKF): the original gyro_h_StageE
% assumed state (q, r, fixed-plane/S-frame rates) and measurement
% (gyro, body-frame rates) were the same quantity. They are not. A real
% gyro is body-mounted and reads body-frame rates; q/r in this EKF's
% state vector are S-frame rates (dtheta_S_dot, dpsi_S_dot).
%
% FIX: S-frame and body frame share the same x-axis (roll axis) and
% differ only by a rotation about it through phi_roll (Stage D,
% external input, same role it already plays in mag_axial_h_StageE and
% Layer 6A/6C's R_i2b construction). So:
%
% [q_body] [ cos(phi_roll) -sin(phi_roll) ] [q_S]
% [r_body] = [ sin(phi_roll) cos(phi_roll) ] [r_S]
%
% h_gyro = [q_body_pred + b_q; r_body_pred + b_r]
%
% This is the SAME Rx(phi_roll) rotation used everywhere else in this
% project's pipeline -- not a new assumption, just correctly applied
% here for the first time.
%
% CONSEQUENCE FOR VALIDATION: Phase 1's Test 2 (unity sensitivity,
% H(1,3)=H(2,4)=1 everywhere) and Test 3 (zero coupling between q/r)
% no longer hold in general -- H(1,3)=cos(phi_roll), H(1,4)=-sin(phi_roll),
% H(2,3)=sin(phi_roll), H(2,4)=cos(phi_roll). These ARE genuinely coupled
% now, by design, because a single body-frame gyro axis reads a mixture
% of both S-frame rates whenever phi_roll != 0. The tests below are
% rewritten to check the CORRECT structure instead.
% =========================================================================
fprintf('=== PHASE 1b VALIDATION: CORRECTED GYRO MEASUREMENT MODEL ===\n\n');

% --- TEST 1: zero roll -> reduces to old Phase 1 behavior exactly ---
x_test1 = [0;0; deg2rad(3); deg2rad(-2); 100; 200; 0.01; -0.005];
h_zero_roll = gyro_h_StageE_v2(x_test1, 0);
q_S = x_test1(3); r_S = x_test1(4); b_q = x_test1(7); b_r = x_test1(8);
expected_zero_roll = [q_S + b_q; r_S + b_r];
fprintf('TEST 1 (phi_roll=0): h_gyro = [%.6f; %.6f] (expect [%.6f; %.6f], same as unrotated case)\n', ...
h_zero_roll(1), h_zero_roll(2), expected_zero_roll(1), expected_zero_roll(2));
if all(abs(h_zero_roll - expected_zero_roll) < 1e-10)
fprintf(' PASS\n\n');
else
fprintf(' *** FAIL -- rotation should collapse to identity at phi_roll=0, do not proceed ***\n\n');
end


% --- TEST 2: 90 deg roll -> q/r should swap (with a sign flip) exactly ---
x_test2 = [0;0; deg2rad(5); deg2rad(-3); 100; 200; 0; 0];
h_90 = gyro_h_StageE_v2(x_test2, pi/2);
q_S = x_test2(3); r_S = x_test2(4);
expected_90 = [-r_S; q_S]; % cos(90)=0, sin(90)=1
fprintf('TEST 2 (phi_roll=90deg): h_gyro = [%.6f; %.6f] (expect [%.6f; %.6f], q/r swapped w/ sign)\n', ...
h_90(1), h_90(2), expected_90(1), expected_90(2));
if all(abs(h_90 - expected_90) < 1e-10)
fprintf(' PASS\n\n');
else
fprintf(' *** FAIL -- 90deg rotation not behaving as a clean swap, check rotation matrix ***\n\n');
end

% --- TEST 3: rotation preserves vector magnitude (energy/norm check) ---
fprintf('--- TEST 3: rotation preserves ||[q;r]|| magnitude across roll angles ---\n');
phi_sweep = deg2rad(0:15:345);
x_test3 = [0;0; deg2rad(4); deg2rad(-6); 100; 200; 0; 0];
q_S = x_test3(3); r_S = x_test3(4);
norm_true = sqrt(q_S^2 + r_S^2);
max_norm_err = 0;
for phi_r = phi_sweep
h_r = gyro_h_StageE_v2(x_test3, phi_r);
norm_h = sqrt(h_r(1)^2 + h_r(2)^2);
max_norm_err = max(max_norm_err, abs(norm_h - norm_true));
end
fprintf('Max |norm(h_gyro) - norm([q_S;r_S])| across full roll sweep: %.3e (expect ~0)\n', max_norm_err);
if max_norm_err < 1e-10
fprintf(' PASS -- rotation is orthonormal, magnitude preserved as expected\n\n');
else
fprintf(' *** FAIL -- rotation matrix is not orthonormal, check for typo ***\n\n');
end

% --- TEST 4: Jacobian structure -- q/r columns now a rotation block,
% NOT a clean diagonal; other columns (dtheta_S, dpsi_S, pb, axial_bias)
% still exactly zero; bias columns still exactly 1 (on their own row) ---
fprintf('--- TEST 4: Jacobian structure check ---\n');
phi_test = deg2rad(37); % arbitrary non-special angle
x_test4 = [deg2rad(1.5); deg2rad(-0.8); deg2rad(5); deg2rad(-3); 180; 150; 0.02; -0.015];
H4 = gyro_jacobian_numeric_v2(x_test4, phi_test);
fprintf('Full H matrix (2x8) at phi_roll=%.1f deg:\n', rad2deg(phi_test));
disp(H4);

expected_H_qr = [cos(phi_test), -sin(phi_test); sin(phi_test), cos(phi_test)];
actual_H_qr = H4(:, 3:4);
fprintf('Expected q/r block:\n'); disp(expected_H_qr);
fprintf('Actual q/r block (numeric):\n'); disp(actual_H_qr);
qr_block_err = max(abs(expected_H_qr - actual_H_qr), [], 'all');
fprintf('Max q/r block error: %.3e (expect ~0)\n', qr_block_err);

independent_cols = [1,2,5,6];
max_leak = max(abs(H4(:,independent_cols)), [], 'all');
fprintf('Max |H| among columns that SHOULD still be zero (dtheta_S,dpsi_S,pb,axial_bias): %.2e\n', max_leak);

bias_diag_err = max(abs([H4(1,7); H4(2,8)] - 1));
bias_offdiag_leak = max(abs([H4(1,8); H4(2,7)]));
fprintf('Bias self-sensitivity error (expect ~0, i.e. H(1,7)=H(2,8)=1): %.2e\n', bias_diag_err);
fprintf('Bias cross-leak (expect ~0, i.e. H(1,8)=H(2,7)=0): %.2e\n', bias_offdiag_leak);

if qr_block_err < 1e-4 && max_leak < 1e-6 && bias_diag_err < 1e-4 && bias_offdiag_leak < 1e-6
fprintf(' PASS -- Jacobian has the correct rotation-block structure\n\n');
else
fprintf(' *** FAIL -- Jacobian structure does not match the intended rotation model ***\n\n');
end

fprintf('\n=== PHASE 1b VALIDATION COMPLETE ===\n');
fprintf('If all four tests PASSED, gyro_h_StageE_v2/gyro_jacobian_numeric_v2\n');
fprintf('should REPLACE gyro_h_StageE8/gyro_jacobian_numeric8 inside\n');
fprintf('attitude_ekf_step_StageE_v8.m -- the gyro update call there must be\n');
fprintf('changed to pass phi_roll through to both h and H (it does not\n');
fprintf('currently, since that file used the un-rotated Phase 1 model).\n');

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function h_gyro = gyro_h_StageE_v2(x, phi_roll)
% CORRECTED gyro measurement model. Rotates the state's S-frame (fixed
% plane) q, r into body frame via phi_roll (Stage D, external input),
% matching the same Rx(phi_roll) convention used in mag_axial_h_StageE
% and Layer 6A/6C's R_i2b construction.
q_S = x(3);
r_S = x(4);
b_q = x(7);
b_r = x(8);

q_body_pred = cos(phi_roll)*q_S - sin(phi_roll)*r_S;
r_body_pred = sin(phi_roll)*q_S + cos(phi_roll)*r_S;

h_gyro = [q_body_pred + b_q; r_body_pred + b_r];
end

function H = gyro_jacobian_numeric_v2(x, phi_roll)
% Numeric Jacobian of gyro_h_StageE_v2 w.r.t. all 8 states, at fixed
% (externally supplied) phi_roll -- same treatment as
% mag_axial_jacobian_numeric holding theta_v/psi_v fixed.
n = length(x);
H = zeros(2, n);
eps_ang = 1e-6;
eps_rate = 1e-6;
eps_bias_mag = 1e-3;
eps_bias_gyro = 1e-6;

h0 = gyro_h_StageE_v2(x, phi_roll);
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
h1 = gyro_h_StageE_v2(x + dx, phi_roll);
H(:,i) = (h1 - h0) / dx(i);
end
end

