function dstate = eom_hifi(state, p, rho, phi_c)
% HIGH-FIDELITY 6-DOF Equations of Motion
% Spin-stabilized projectile — rigid body model
% Based on:
%   Baranowski (2013) Journal of Theoretical and Applied Mechanics
%   Gkritzapis et al. (2007) AIAA Paper 2007-6584
%   Fresconi et al. (2011) Journal of Aerospace Engineering
%
% State vector (13 elements):
%   state(1:3)  = position in inertial frame [m]
%   state(4:6)  = velocity in inertial frame [m/s]
%   state(7:10) = quaternion [q0,q1,q2,q3] (scalar first)
%   state(11:13)= body angular rates [p,q,r] [rad/s]
%                 p = spin (high), q = pitch rate, r = yaw rate
%
% Coordinate convention:
%   Inertial frame: x=North, y=East, z=Down (gun frame)
%   Body frame: x=spin axis (nose), y,z = lateral

% ================================================
% UNPACK STATE
% ================================================
pos_i = state(1:3);    % position [m]
vel_i = state(4:6);    % velocity inertial [m/s]
q0    = state(7);      % quaternion scalar part
q1    = state(8);
q2    = state(9);
q3    = state(10);
p_b   = state(11);     % spin rate [rad/s] — large ~1257
q_b   = state(12);     % pitch rate [rad/s] — small
r_b   = state(13);     % yaw rate [rad/s] — small

V = norm(vel_i);
if V < 1; V = 1; end
qbar = 0.5 * rho * V^2;

% ================================================
% QUATERNION TO ROTATION MATRIX
% R_i2b: inertial to body (maps inertial vectors to body frame)
% Convention: q = [q0, q1, q2, q3] scalar-first
% Reference: Baranowski (2013) Eq.4.17
% ================================================
R_i2b = 2 * [...
    0.5*(q0^2+q1^2-q2^2-q3^2), q1*q2+q0*q3,            q1*q3-q0*q2;
    q1*q2-q0*q3,                0.5*(q0^2-q1^2+q2^2-q3^2), q2*q3+q0*q1;
    q1*q3+q0*q2,                q2*q3-q0*q1,            0.5*(q0^2-q1^2-q2^2+q3^2)];

% Body to inertial = transpose
R_b2i = R_i2b';

% ================================================
% VELOCITY IN BODY FRAME
% This gives true aerodynamic angles
% ================================================
vel_b = R_i2b * vel_i;
u_b = vel_b(1);   % axial velocity component
v_b = vel_b(2);   % lateral y component
w_b = vel_b(3);   % lateral z component

% ================================================
% AERODYNAMIC ANGLES
% alpha = pitch-plane angle of attack
% beta  = yaw-plane sideslip
% alpha_t = total angle of attack
% Reference: Gkritzapis (2007), McCoy Chapter 4
% ================================================
alpha_t = asin(max(-1, min(1, sqrt(v_b^2 + w_b^2) / V)));
if alpha_t < 1e-6
    phi_aero = 0;
    
else
    phi_aero = atan2(w_b, v_b);   % aerodynamic roll angle
end

% Small-angle alpha and beta for moment equations
alpha = atan2(w_b, u_b);
beta  = atan2(v_b, u_b);
assignin('base','alpha_debug',rad2deg(alpha));
assignin('base','beta_debug',rad2deg(beta));
assignin('base','u_debug',u_b);
assignin('base','v_debug',v_b);
assignin('base','w_debug',w_b);

% ================================================
% CANARD AOA (Patent Eq.11)
% ================================================
alpha_c = asin(max(-1, min(1, ...
    (cos(phi_c)*w_b - sin(phi_c)*v_b) / V)));

% ================================================
% AERODYNAMIC FORCES IN BODY FRAME
% Reference: Fresconi et al. (2011), Patent Eq.2-3
% Forces expressed in fixed-plane (-f) frame then
% transformed — here we work directly in body frame
% ================================================

% Axial drag force
X_B = qbar * p.A * (p.CD0 + p.CX_a2 * alpha_t^2);

% Normal force (y-body)
Y_B = qbar * p.A * (-p.CNa * beta + ...
      p.CYpa * p_b * p.d / (2*V));

% Normal force (z-body)
Z_B = qbar * p.A * (-p.CNa * alpha + ...
      p.CYpa * p_b * p.d / (2*V));

% Canard forces
X_C = qbar * p.A * p.CX0c;
Y_C = qbar * p.A * (p.CZc  * sin(phi_c) + ...
      p.CYac * alpha_c * cos(phi_c));
Z_C = qbar * p.A * (-p.CZc * cos(phi_c) + ...
      p.CYac * alpha_c * sin(phi_c));


% Total aerodynamic force in body frame
F_aero_b = [-(X_B + X_C); (Y_B + Y_C); (Z_B + Z_C)];

% ================================================
% GRAVITY IN BODY FRAME
% Gravity is [0; 0; g] in inertial (z-down)
% Rotate to body frame
% ================================================
g_i = [0; 0; p.m * p.g];
F_grav_b = R_i2b * g_i;

% Total force in body frame
F_total_b = F_aero_b + F_grav_b;

% ================================================
% TRANSLATIONAL ACCELERATION
% Newton: F = m*a in body frame
% Then rotate to inertial for position integration
% Reference: Gkritzapis (2007)
% ================================================
% Acceleration in body frame (includes centripetal terms
% because body frame is rotating)
% a_inertial = R_b2i * (F/m) is cleaner for integration
accel_i = R_b2i * (F_total_b / p.m);

% ================================================
% AERODYNAMIC MOMENTS IN BODY FRAME
% Reference: Fresconi et al. (2011) Eq.5, Patent Eq.5
% ================================================

% Roll (spin) damping moment
L_B = qbar * p.A * p.d * ...
      (p.Clp * p_b * p.d / (2*V));

% Pitch moment
M_B = qbar * p.A * p.d * ...
      (p.Cma  * alpha + ...
       p.Cmq  * q_b * p.d / (2*V) + ...
       p.Cnpa * p_b * p.d / (2*V) * beta);

% Yaw moment
N_B = qbar * p.A * p.d * ...
      (-p.Cma  * beta + ...
       p.Cmq   * r_b * p.d / (2*V) + ...
       p.Cnpa  * p_b * p.d / (2*V) * alpha);

% Canard moments
M_C = qbar * p.A * p.d * ...
      ((p.CG - p.CPz) * p.CZc  * cos(phi_c) - ...
       (p.CG - p.CPy) * p.CYac * alpha_c * sin(phi_c));

N_C = qbar * p.A * p.d * ...
      ((p.CG - p.CPz) * p.CZc  * sin(phi_c) + ...
       (p.CG - p.CPy) * p.CYac * alpha_c * cos(phi_c));

% Total moment vector in body frame
L_total = L_B;
M_total = M_B + M_C;
N_total = N_B + N_C;

% ================================================
% ROTATIONAL DYNAMICS — Euler's equation
% I*omega_dot + omega x (I*omega) = M_total
% Reference: Baranowski (2013) Eq.4.7
% For axisymmetric projectile: Iy = Iz
% I = diag(Ix, Iy, Iy)
% ================================================
Ix = p.Ix;
Iy = p.Iy;

% omega x (I*omega) gyroscopic term
% omega = [p_b; q_b; r_b]
% I*omega = [Ix*p_b; Iy*q_b; Iy*r_b]
gyro_p = q_b * (Iy*r_b) - r_b * (Iy*q_b);   % = 0 for Iy=Iz
gyro_q = r_b * (Ix*p_b) - p_b * (Iy*r_b);
gyro_r = p_b * (Iy*q_b) - q_b * (Ix*p_b);

% Angular accelerations: I*omega_dot = M - gyro
p_dot = (L_total - gyro_p) / Ix;
q_dot = (M_total - gyro_q) / Iy;
r_dot = (N_total - gyro_r) / Iy;

% ================================================
% QUATERNION KINEMATICS
% Reference: Baranowski (2013) Eq.4.19
% dq/dt = 0.5 * Omega * q
% where Omega is skew-symmetric matrix of angular rates
% ================================================
q0_dot = 0.5 * (-q1*p_b - q2*q_b - q3*r_b);
q1_dot = 0.5 * ( q0*p_b - q3*q_b + q2*r_b);
q2_dot = 0.5 * ( q3*p_b + q0*q_b - q1*r_b);
q3_dot = 0.5 * (-q2*p_b + q1*q_b + q0*r_b);

% ================================================
% NaN SAFETY CHECK
% ================================================
dstate_check = [accel_i; q0_dot; q1_dot; q2_dot; q3_dot; ...
                p_dot; q_dot; r_dot];
if any(isnan(dstate_check)) || any(isinf(dstate_check))
    fprintf('HiFi NaN at V=%.1f alpha_t=%.2f\n', V, rad2deg(alpha_t));
    accel_i = zeros(3,1);
    q0_dot=0; q1_dot=0; q2_dot=0; q3_dot=0;
    p_dot=0;  q_dot=0;  r_dot=0;
end

% ================================================
% ASSEMBLE DERIVATIVE (13 elements)
% ================================================
dstate = [vel_i;                          % pos_dot = vel
          accel_i;                         % vel_dot = accel
          q0_dot; q1_dot; q2_dot; q3_dot; % quaternion kinematics
          p_dot; q_dot; r_dot];            % angular rate dynamics
end