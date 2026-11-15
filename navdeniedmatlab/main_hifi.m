% =========================================
% HIGH-FIDELITY TRUTH MODEL
% Full 6-DOF quaternion-based simulation
% Generates sensor truth for navigation EKF
% Patent US 9,702,674 B2
% =========================================
clear; close all; clc;

p = projectile_params();

% --- Launch conditions ---
elev_rad = deg2rad(p.elev);
az_rad   = deg2rad(p.az);

u0 = p.V0 * cos(elev_rad) * cos(az_rad);
v0 = p.V0 * cos(elev_rad) * sin(az_rad);
w0 = -p.V0 * sin(elev_rad);

% --- Initial quaternion from Euler angles ---
% phi0=0, theta0=elev_rad, psi0=az_rad
% Convert Euler to quaternion (ZYX convention)
phi0   = 0;
theta0 = elev_rad;
psi0   = az_rad;

q0_init = cos(phi0/2)*cos(theta0/2)*cos(psi0/2) + ...
          sin(phi0/2)*sin(theta0/2)*sin(psi0/2);
q1_init = sin(phi0/2)*cos(theta0/2)*cos(psi0/2) - ...
          cos(phi0/2)*sin(theta0/2)*sin(psi0/2);
q2_init = cos(phi0/2)*sin(theta0/2)*cos(psi0/2) + ...
          sin(phi0/2)*cos(theta0/2)*sin(psi0/2);
q3_init = cos(phi0/2)*cos(theta0/2)*sin(psi0/2) - ...
          sin(phi0/2)*sin(theta0/2)*cos(psi0/2);

% --- Initial state (13 elements) ---
% [x,y,z, u,v,w, q0,q1,q2,q3, p,q,r]
state_hifi = [0; 0; 0; ...
              u0; v0; w0; ...
              q0_init; q1_init; q2_init; q3_init; ...
              p.p0; 0; 0];

phi_c = pi;


% --- Pre-allocate ---
N = round(p.t_max / p.dt) + 1;
alpha_log = zeros(N,1);%debug
beta_log = zeros(N,1);%debug
time_log   = zeros(N, 1);
state_log  = zeros(N, 13);
euler_log  = zeros(N, 3);   % [phi, theta, psi] for sensor models
mag_log    = zeros(N, 3);
accel_log  = zeros(N, 3);
gps_log    = zeros(N, 6);

time_log(1)    = 0;
state_log(1,:) = state_hifi';

% Extract initial Euler angles
euler_log(1,:) = [phi0, theta0, psi0];
k = 1;
t = 0;

% --- Main loop ---
while t < p.t_max

    altitude = -state_hifi(3);
    if t > 0.1 && altitude <= 0
        break;
    end

    [rho, ~] = atmosphere(max(altitude, 0));

    % RK4 integration of high-fidelity model
    k1 = eom_hifi(state_hifi, p, rho, phi_c);
    k2 = eom_hifi(state_hifi + 0.5*p.dt*k1, p, rho, phi_c);
    k3 = eom_hifi(state_hifi + 0.5*p.dt*k2, p, rho, phi_c);
    k4 = eom_hifi(state_hifi + p.dt*k3,     p, rho, phi_c);

    state_hifi = state_hifi + (p.dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    alpha_log(k) = alpha_debug;
    beta_log(k)  = beta_debug;
    % Re-normalize quaternion to prevent drift
    q_norm = norm(state_hifi(7:10));
    if q_norm > 0
        state_hifi(7:10) = state_hifi(7:10) / q_norm;
    end

    % Extract Euler angles from quaternion for sensor models
    q0 = state_hifi(7); q1 = state_hifi(8);
    q2 = state_hifi(9); q3 = state_hifi(10);

    phi_e   = atan2(2*(q0*q1 + q2*q3), ...
                    1 - 2*(q1^2 + q2^2));
    theta_e = asin(max(-1, min(1, 2*(q0*q2 - q3*q1))));
    psi_e   = atan2(2*(q0*q3 + q1*q2), ...
                    1 - 2*(q2^2 + q3^2));

    % Build 12-element state for sensor_models
    % sensor_models expects [x,y,z, u,v,w, phi,theta,psi, pb, q, r]
    state_for_sensors = [state_hifi(1:6); ...
                         phi_e; theta_e; psi_e; ...
                         state_hifi(11:13)];

    % Generate sensor measurements from hifi truth
    % Generate sensor measurements from hifi truth (use the hifi-specific model)
    sensors = sensor_models_hifi(state_hifi, p, t);

    t = t + p.dt;
    k = k + 1;

    time_log(k)    = t;
    state_log(k,:) = state_hifi';
    euler_log(k,:) = [phi_e, theta_e, psi_e];

    % --- Sensor logging (add these) ---
    mag_log(k,:)   = sensors.mag_body';
    accel_log(k,:) = sensors.accel_body';
    if sensors.gps_available
        gps_log(k,:) = [sensors.gps_pos; sensors.gps_vel]';
    else
        gps_log(k,:) = NaN(1,6);
    end
end

% Trim
time_log  = time_log(1:k);
state_log = state_log(1:k,:);
euler_log = euler_log(1:k,:);
mag_log   = mag_log(1:k,:);
accel_log = accel_log(1:k,:);
gps_log   = gps_log(1:k,:);
alpha_log = alpha_log(1:k);
beta_log  = beta_log(1:k);


fprintf('=== HiFi Truth Model Results ===\n');
fprintf('Flight time: %.2f s\n', t);
fprintf('Range:       %.2f km\n', state_log(end,1)/1000);
fprintf('Max alt:     %.2f km\n', max(-state_log(:,3))/1000);
fprintf('Final spin:  %.2f Hz\n', state_log(end,11)/(2*pi));
fprintf('East drift:  %.2f m\n',  state_log(end,2));
fprintf('Final phi:   %.2f deg\n', rad2deg(euler_log(end,1)));
fprintf('Final theta: %.2f deg\n', rad2deg(euler_log(end,2)));
fprintf('Final q:     %.4f rad/s\n', state_log(end,12));
fprintf('Final r:     %.4f rad/s\n', state_log(end,13));
fprintf('\n');
fprintf('Max alpha =%.2f deg\n',max(abs(alpha_log)));   
fprintf('Max beta =%.2f deg\n',max(abs(beta_log))); 
% --- Plots ---
figure('Name','HiFi Truth Model','Position',[100 100 1400 800]);

subplot(2,3,1);
plot(state_log(:,1)/1000, -state_log(:,3)/1000, 'b', 'LineWidth', 1.5);
xlabel('Range [km]'); ylabel('Altitude [km]');
title('Trajectory'); grid on;

V_all = sqrt(sum(state_log(:,4:6).^2, 2));
subplot(2,3,2);
plot(time_log, V_all, 'r', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Speed [m/s]');
title('Speed vs Time'); grid on;

subplot(2,3,3);
plot(time_log, state_log(:,11)/(2*pi), 'k', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Spin [Hz]');
title('Spin rate'); grid on;

subplot(2,3,4);
plot(time_log, rad2deg(euler_log(:,1)), 'g', 'LineWidth', 0.5);
xlabel('Time [s]'); ylabel('\phi [deg]');
title('Roll accumulation'); grid on;

subplot(2,3,5);
plot(time_log, rad2deg(euler_log(:,2)), 'm', 'LineWidth', 1.5);
hold on;
plot(time_log, rad2deg(euler_log(:,3)), 'c--', 'LineWidth', 1.5);
legend('\theta body pitch', '\psi body yaw');
xlabel('Time [s]'); ylabel('Angle [deg]');
title('Body Pitch and Yaw'); grid on;

subplot(2,3,6);
plot3(state_log(:,1)/1000, state_log(:,2)/1000, ...
      -state_log(:,3)/1000, 'b', 'LineWidth', 1.5);
xlabel('North [km]'); ylabel('East [km]');
zlabel('Alt [km]');
title('3D Path'); grid on;

% ================================================
% SENSOR PLOTS
% ================================================
figure('Name','HiFi Layer 5: Sensor Outputs','Position',[100 100 1400 900]);

subplot(3,3,1);
plot(time_log, mag_log(:,1), 'b', 'LineWidth', 0.5);
xlabel('Time [s]'); ylabel('B_x [nT]');
title('Magnetometer X (axial)'); grid on;

subplot(3,3,2);
plot(time_log, mag_log(:,2), 'r', 'LineWidth', 0.5);
xlabel('Time [s]'); ylabel('B_y [nT]');
title('Magnetometer Y (lateral)'); grid on;

subplot(3,3,3);
plot(time_log, mag_log(:,3), 'g', 'LineWidth', 0.5);
xlabel('Time [s]'); ylabel('B_z [nT]');
title('Magnetometer Z (lateral)'); grid on;

subplot(3,3,4);
plot(time_log, accel_log(:,1), 'b', 'LineWidth', 0.5);
xlabel('Time [s]'); ylabel('a_x [m/s²]');
title('Accelerometer X (axial)'); grid on;

subplot(3,3,5);
plot(time_log, accel_log(:,2), 'r', 'LineWidth', 0.5);
xlabel('Time [s]'); ylabel('a_y [m/s²]');
title('Accelerometer Y (radial)'); grid on;

subplot(3,3,6);
plot(time_log, accel_log(:,3), 'g', 'LineWidth', 0.5);
xlabel('Time [s]'); ylabel('a_z [m/s²]');
title('Accelerometer Z (radial)'); grid on;

subplot(3,3,7);
gps_valid = ~isnan(gps_log(:,1));
plot(time_log, gps_valid, 'k', 'LineWidth', 2);
xlabel('Time [s]'); ylabel('Available');
title('GPS Available (1=yes, 0=denied)'); grid on;
ylim([-0.1 1.1]);
xline(p.t_gps_end, '--r', 'GPS Denied');

subplot(3,3,8);
gps_err = sqrt(sum((gps_log(:,1:3) - state_log(:,1:3)).^2, 2));
plot(time_log(gps_valid), gps_err(gps_valid), 'm', 'LineWidth', 1);
xlabel('Time [s]'); ylabel('Error [m]');
title('GPS Position Error'); grid on;

% --- Nuisance parameter: spin is column 11 in hifi state ---
subplot(3,3,9);
centripetal = p.accel_offset_j * state_log(:,11).^2;
plot(time_log, centripetal, 'r', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Centripetal accel [m/s²]');
title('Nuisance Parameter (c×p_b²)'); grid on;

figure('Name','HiFi NP Component Analysis','Position',[100 100 1400 600]);

subplot(1,3,1);
plot(time_log(1:100), accel_log(1:100,2), 'r');
xlabel('Time [s]'); ylabel('[m/s²]');
title('Accel Y — first 0.01s (see sinusoid)'); grid on;

subplot(1,3,2);
plot(time_log, accel_log(:,2), 'r', 'LineWidth', 0.5);
xlabel('Time [s]'); ylabel('[m/s²]');
title('Accel Y — full flight'); grid on;

subplot(1,3,3);
cp_j = p.accel_offset_j * state_log(:,11).^2;
plot(time_log, cp_j, 'b', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('[m/s²]');
title('c×pb² centripetal decay'); grid on;