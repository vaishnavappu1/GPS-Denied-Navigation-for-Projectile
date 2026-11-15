function [pos_new, vel_new, Phi_d, P_d] = track1_dynamic_step(...
    pos_prev, vel_prev, P_prev_combined, ...
    dtheta_S_est, dpsi_S_est, pb_est, q_est, r_est, ...
    B_earth, rho, dt, p, Q_d)
% TRACK 1: Flight-Dynamic Position/Velocity Predictor
% Patent US 9,702,674 B2, claims 7-9, FIG.6 block 640, FIG.7 block 750, FIG.8
%
% DIGITAL ANTI-ALIASING INTEGRATION GUARD REGISTERED THIS BLOCK.

phi_c = 0;  % unguided baseline

% ================================================
% BASE STEP: Module 750 Composition
% ================================================
[alpha_est, beta_est, R_g2f_base] = compute_alpha_beta_module750(...
    vel_prev, dtheta_S_est, dpsi_S_est, B_earth);

% Acceleration at current state using the unified frame
a_prev = compute_accel(vel_prev, alpha_est, beta_est, pb_est, phi_c, p, rho, R_g2f_base, dt);

% ================================================
% TRAPEZOIDAL INTEGRATION (patent eq.12)
vel_mid = vel_prev + dt * a_prev;
[~, ~, R_g2f_mid] = compute_alpha_beta_module750(vel_mid, dtheta_S_est, dpsi_S_est, B_earth);
a_mid = compute_accel(vel_mid, alpha_est, beta_est, pb_est, phi_c, p, rho, R_g2f_mid, dt);

vel_new = vel_prev + 0.5*dt*(a_prev + a_mid);
pos_new = pos_prev + 0.5*dt*(vel_prev + vel_new);

% ================================================
% ANALYTIC PHI_D (Jacobian evaluated with frozen R_g2f_base)
% ================================================
eps_v  = 1e-4;
J_vel  = zeros(3,3);
for j = 1:3
    dv       = zeros(3,1);
    dv(j)    = eps_v;
    vel_test = vel_prev + dv;
    
    % The perturbed velocity alters alpha/beta because vel_f = R_g2f * vel_test
    [alpha_test, beta_test, ~] = compute_alpha_beta_module750(...
        vel_test, dtheta_S_est, dpsi_S_est, B_earth);
        
    % CRITICAL: Keep R_g2f_base FROZEN during the perturbation loop to 
    % cleanly isolate the aerodynamic partial derivatives w.r.t velocity.
    a_pert     = compute_accel(vel_test, alpha_test, beta_test, pb_est, phi_c, p, rho, R_g2f_base, dt);
    J_vel(:,j) = (a_pert - a_prev) / eps_v;
end

Phi_d = [eye(3),    dt*eye(3)        ;
         zeros(3,3), eye(3) + dt*J_vel];

% ================================================
% COVARIANCE PROPAGATION (patent eq.15)
% ================================================
P_d     = Phi_d * P_prev_combined * Phi_d' + Q_d;
P_d     = (P_d + P_d') / 2;  % symmetry guard
end

% ================================================
function [alpha, beta, R_g2f] = compute_alpha_beta_module750(vel_i, dtheta_S, dpsi_S, B_earth)
% Patent FIG.7 block 750 composition: R_g2f = R_s2f * R_i2S
[R_i2S, ~] = build_S_frame(vel_i, B_earth);
dtheta_S = dtheta_S;   % <-- ADD THIS LINE
dpsi_S   = -dpsi_S;     % <-- ADD THIS LINE

Ry_s = [cos(dtheta_S), 0, -sin(dtheta_S); 0, 1, 0; sin(dtheta_S), 0, cos(dtheta_S)];
Rz_s = [cos(dpsi_S), sin(dpsi_S), 0; -sin(dpsi_S), cos(dpsi_S), 0; 0, 0, 1];
R_s2f = Ry_s * Rz_s;   % back to original
R_g2f = R_s2f * R_i2S;
vel_f = R_g2f * vel_i;
norm_v = norm(vel_f);
if norm_v < 1; norm_v = 1; end

alpha = asin(max(-1, min(1, vel_f(3)/norm_v)));
beta  = asin(max(-1, min(1, vel_f(2)/norm_v)));
end

function a_inertial = compute_accel(vel_i, alpha, beta, pb, phi_c, p, rho, R_g2f, dt)
V = norm(vel_i);
if V < 1; V = 1; end
qbar = 0.5 * rho * V^2;

% Use the explicitly passed aerodynamic angles to preserve analytical derivatives
alpha_total = sqrt(alpha^2 + beta^2);

% Resolve the baseline command angle directly from the verified input states
alpha_c = alpha * cos(phi_c) - beta * sin(phi_c);

% Aerodynamic forces (Resolved strictly from input tracking states)
X_B = qbar * p.A * (p.CD0 + p.CX_a2 * alpha_total^2);
Y_B = qbar * p.A * (-p.CNa * beta  + p.CYpa * pb * p.d / (2*V));
Z_B = qbar * p.A * (-p.CNa * alpha + p.CYpa * pb * p.d / (2*V));
X_C = qbar * p.A * p.CX0c;
Y_C = qbar * p.A * (p.CZc  * sin(phi_c) + p.CYac * alpha_c * cos(phi_c));
Z_C = qbar * p.A * (-p.CZc * cos(phi_c) + p.CYac * alpha_c * sin(phi_c));

F_f        = [-(X_B+X_C); (Y_B+Y_C); (Z_B+Z_C)];
F_i        = R_g2f' * F_f;
a_inertial = (F_i + [0; 0; p.m*p.g]) / p.m;
end