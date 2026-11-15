function dstate = eom(state, p, rho, phi_c)
% Equations of Motion — Patent US 9,702,674 B2
% Fixed-plane (-f) frame implementation
% This is exactly what the patent uses for its simulation
%
% KEY DESIGN DECISIONS (from patent text):
% 1. Fixed-plane frame tracks velocity vector direction
% 2. alpha and beta are near zero (axisymmetric assumption)
% 3. Spin dynamics integrated separately from trajectory
% 4. No full Euler angle body tracking needed
% --- Unpack state ---
vel_i = state(4:6);   % velocity in inertial frame [m/s]
phi   = state(7);     % roll angle [rad] — for sensor models only
theta = state(8);     % flight path elevation [rad]
psi   = state(9);     % flight path azimuth [rad]
pb    = state(10);    % spin rate [rad/s]
q     = state(11);    % pitch rate of fixed-plane frame [rad/s]
r     = state(12);    % yaw rate of fixed-plane frame [rad/s]
V    = norm(vel_i);
if V < 1; V = 1; end
qbar = 0.5 * rho * V^2;
% --- Fixed-plane DCM from velocity vector ---
% This is the patent's primary coordinate frame
% psi and theta track the velocity direction
psi_v   = atan2(vel_i(2), vel_i(1));
theta_v = atan2(-vel_i(3), ...
                sqrt(vel_i(1)^2 + vel_i(2)^2));
R_i2f = dcm(psi_v, theta_v);
vel_f = R_i2f * vel_i;  % = [V; ~0; ~0]
% --- Small angle of attack in fixed-plane frame ---
% These are near zero for stabilized shell
% Represent residual misalignment only
alpha = atan2(vel_f(3), vel_f(1));
beta  = atan2(vel_f(2), vel_f(1));
alpha_total = sqrt(alpha^2 + beta^2);
% --- Canard AOA (Patent Eq.11) ---
alpha_c = asin(max(-1, min(1, ...
    (cos(phi_c)*vel_f(3) - sin(phi_c)*vel_f(2)) / V)));
% --- Body forces in fixed-plane frame (Patent Eq.2) ---
X_B = qbar * p.A * (p.CD0 + p.CX_a2 * alpha_total^2);
Y_B = qbar * p.A * (-p.CNa * beta + ...
                     p.CYpa * pb * p.d / (2*V));
Z_B = qbar * p.A * (-p.CNa * alpha);
% --- Canard forces (Patent Eq.3) ---
X_C = qbar * p.A * p.CX0c;
Y_C = qbar * p.A * (p.CZc  * sin(phi_c) + ...
                     p.CYac * alpha_c * cos(phi_c));
Z_C = qbar * p.A * (-p.CZc * cos(phi_c) + ...
                     p.CYac * alpha_c * sin(phi_c));
% Total force — rotate to inertial
F_f = [-(X_B+X_C); (Y_B+Y_C); (Z_B+Z_C)];
F_i = R_i2f' * F_f;
accel_i = (F_i + [0;0;p.m*p.g]) / p.m;
% --- Moments in fixed-plane frame (Patent Eq.5) ---
L_B = qbar * p.A * p.d * ...
      (p.Clp * pb * p.d / (2*V));
M_B = qbar * p.A * p.d * ...
      (p.Cma  * alpha + ...
       p.Cmq  * q * p.d / (2*V) + ...
       p.Cnpa * pb * p.d / (2*V) * beta);
N_B = qbar * p.A * p.d * ...
      (-p.Cma  * beta + ...
       p.Cmq   * r * p.d / (2*V) + ...
       p.Cnpa  * pb * p.d / (2*V) * alpha);
M_C = qbar * p.A * p.d * ...
      ((p.CG-p.CPz)*p.CZc *cos(phi_c) - ...
       (p.CG-p.CPy)*p.CYac*alpha_c*sin(phi_c));
N_C = qbar * p.A * p.d * ...
      ((p.CG-p.CPz)*p.CZc *sin(phi_c) + ...
       (p.CG-p.CPy)*p.CYac*alpha_c*cos(phi_c));
% --- Spin decay (Patent Eq.4 line 1) ---
pb_dot = L_B / p.Ix;
% --- Fixed-plane angular rates (Patent Eq.4 lines 2-3) ---
% q and r are rates of the fixed-plane frame itself
% NOT body Euler rates
% Gyroscopic terms arise from spin-rate coupling
q_dot = pb * r * (p.Ix/p.Iy) + (M_B+M_C) / p.Iy;
r_dot = -pb * q * (p.Ix/p.Iy) + (N_B+N_C) / p.Iy;
% --- Roll accumulation ---
phi_dot = pb;
% --- Frame angle rates from acceleration ---
% theta and psi track velocity vector
% Their rates come from how acceleration rotates velocity
a_f = R_i2f * accel_i;
theta_dot = -a_f(3) / V;
cos_theta = cos(theta_v);
if abs(cos_theta) < 0.01
    cos_theta = 0.01;
end
psi_dot = a_f(2) / (V * cos_theta);
% --- NaN check ---
if any(isnan([accel_i; phi_dot; theta_dot; psi_dot; ...
              pb_dot; q_dot; r_dot]))
    accel_i = zeros(3,1);
    phi_dot = pb; theta_dot = 0; psi_dot = 0;
    pb_dot  = 0;  q_dot = 0;    r_dot = 0;
end
dstate = [vel_i; accel_i; ...
          phi_dot; theta_dot; psi_dot; ...
          pb_dot; q_dot; r_dot];
end