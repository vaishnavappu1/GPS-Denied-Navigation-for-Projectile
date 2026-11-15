function sensors = sensor_models_hifi(state_hifi, p, t)
% SENSOR MODELS for HIGH-FIDELITY (quaternion) truth model
% state_hifi: 13-element [x,y,z,u,v,w(inertial),q0,q1,q2,q3,pb,q,r]

pos_i = state_hifi(1:3);
vel_i = state_hifi(4:6);
q0=state_hifi(7); q1=state_hifi(8); q2=state_hifi(9); q3=state_hifi(10);
pb = state_hifi(11); q = state_hifi(12); r = state_hifi(13);

V = norm(vel_i);
if V < 1; V = 1; end

% --- SAME quaternion->DCM as eom_hifi for consistency ---
R_i2b = 2 * [...
    0.5*(q0^2+q1^2-q2^2-q3^2), q1*q2+q0*q3,            q1*q3-q0*q2;
    q1*q2-q0*q3,                0.5*(q0^2-q1^2+q2^2-q3^2), q2*q3+q0*q1;
    q1*q3+q0*q2,                q2*q3-q0*q1,            0.5*(q0^2-q1^2-q2^2+q3^2)];

% True body-frame velocity and AOA
vel_b = R_i2b * vel_i;
u_b=vel_b(1); v_b=vel_b(2); w_b=vel_b(3);
alpha = atan2(w_b,u_b);
beta  = atan2(v_b,u_b);
alpha_total = sqrt(alpha^2+beta^2);

[rho,~] = atmosphere(max(-pos_i(3),0));
qbar = 0.5*rho*V^2;

% True specific force in body frame (forces computed directly in body)
f_x = -(qbar*p.A*(p.CD0+p.CX_a2*alpha_total^2))/p.m;
f_y =  (qbar*p.A*(-p.CNa*beta + p.CYpa*pb*p.d/(2*V)))/p.m;
f_z =  (qbar*p.A*(-p.CNa*alpha))/p.m;
f_body_true = [f_x; f_y; f_z];

% --- Accelerometer NP model (same as before) ---
phi_e = atan2(2*(q0*q1+q2*q3), 1-2*(q1^2+q2^2));  % for spin_angle
spin_angle = phi_e;

a_j = p.balloting_j*exp(-t/p.balloting_tau);
a_k = p.balloting_k*exp(-t/p.balloting_tau);
NP_sin_j = a_j*sin(spin_angle);
NP_sin_k = a_k*cos(spin_angle);

b_i=p.bias_accel(1)*p.g; b_j=p.bias_accel(2)*p.g; b_k=p.bias_accel(3)*p.g;
c_j=p.accel_offset_j*pb^2; c_k=p.accel_offset_k*pb^2; c_i=0.001*(q^2+r^2);

NP_i=b_i+c_i; NP_j=NP_sin_j+b_j+c_j; NP_k=NP_sin_k+b_k+c_k;
noise_a = p.sigma_accel*p.g*randn(3,1);

sensors.accel_body = [f_body_true(1)+NP_i+noise_a(1);
    f_body_true(2)+NP_j+noise_a(2);
    f_body_true(3)+NP_k+noise_a(3)];
sensors.accel_body_true = f_body_true;
sensors.accel_NP = [NP_i;NP_j;NP_k];

% --- Magnetometer (same R_i2b as eom_hifi) ---
B_body_true = R_i2b*p.B_earth;
noise_m = p.sigma_mag*randn(3,1);
sensors.mag_body = p.S_mag*B_body_true + p.bias_mag + noise_m;
sensors.mag_body_true = B_body_true;

% --- GPS ---
if t<=p.t_gps_end
    sensors.gps_pos=pos_i+p.sigma_gps_pos*randn(3,1);
    sensors.gps_vel=vel_i+p.sigma_gps_vel*randn(3,1);
    sensors.gps_available=true;
else
    sensors.gps_pos=NaN(3,1); sensors.gps_vel=NaN(3,1);
    sensors.gps_available=false;
end

sensors.phi_nose = 0.01*randn();
sensors.time=t; sensors.true_pos=pos_i; sensors.true_vel=vel_i;
sensors.true_pb=pb; sensors.true_alpha=alpha; sensors.true_beta=beta;
end