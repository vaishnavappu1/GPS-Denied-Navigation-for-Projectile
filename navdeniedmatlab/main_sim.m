
% =========================================
% LAYER 4B: Full 6-DOF with RK4
% State expanded to 12: adds q and r
% Patent US 9,702,674 B2
% =========================================
clear; close all; clc;

p = projectile_params();

elev_rad = deg2rad(p.elev);
az_rad   = deg2rad(p.az);

u0_i = p.V0 * cos(elev_rad) * cos(az_rad);
v0_i = p.V0 * cos(elev_rad) * sin(az_rad);
w0_i = -p.V0 * sin(elev_rad);

% State: [x,y,z, u,v,w, phi,theta,psi, pb, q, r]
state = [0; 0; 0; ...
         u0_i; v0_i; w0_i; ...
         0; elev_rad; az_rad; ...
         p.p0; 0; 0];




N = round(p.t_max / p.dt) + 1;
time_log  = zeros(N, 1);
state_log = zeros(N, 12);
time_log(1)    = 0;
state_log(1,:) = state';
k = 1;

t = 0;
mag_log = zeros(N, 3);
accel_log = zeros(N, 3);
gps_log = zeros(N, 6);
phi_c = 0; %20*pi/180;   % rad
while t < p.t_max

    pos_i    = state(1:3);
    altitude = -pos_i(3);
    if t > 0.1 && altitude <= 0
        break;
    end

    [rho, ~] = atmosphere(altitude);

    % RK4

    k1 = eom(state, p, rho, phi_c);
    k2 = eom(state + 0.5*p.dt*k1, p, rho, phi_c);
    k3 = eom(state + 0.5*p.dt*k2, p, rho, phi_c);
    k4 = eom(state + p.dt*k3,     p, rho, phi_c);

    state = state + (p.dt/6)*(k1 + 2*k2 + 2*k3 + k4);

    t = t + p.dt;
    k = k + 1;
    time_log(k)    = t;
    state_log(k,:) = state';
end

time_log  = time_log(1:k);
state_log = state_log(1:k,:);
alpha_log = zeros(N,1);
beta_log  = zeros(N,1);



fprintf('=== Layer 4B Results ===\n');
fprintf('Flight time: %.2f s\n', t);
fprintf('Range:       %.2f km\n', state_log(end,1)/1000);
fprintf('Max alt:     %.2f km\n', max(-state_log(:,3))/1000);
fprintf('Final spin: %.2f Hz\n', state_log(end,10)/(2*pi));
fprintf('East drift:  %.2f m\n',  state_log(end,2));
fprintf('Final q:     %.4f rad/s\n', state_log(end,11));
fprintf('Final r:     %.4f rad/s\n', state_log(end,12));

% Plots
figure('Name','Layer 4B: Full 6-DOF RK4', ...
       'Position',[100 100 1400 800]);

subplot(2,3,1);
plot(state_log(:,1)/1000, -state_log(:,3)/1000, ...
     'b', 'LineWidth', 1.5);
xlabel('Range [km]'); ylabel('Altitude [km]');
title('Trajectory'); grid on;

V_all = sqrt(sum(state_log(:,4:6).^2, 2));
subplot(2,3,2);
plot(time_log, V_all, 'r', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Speed [m/s]');
title('Speed vs Time'); grid on;

subplot(2,3,3);
plot(time_log, state_log(:,10) / (2*pi), 'r', 'LineWidth', 1.5);
ylabel('Spin [Hz]');
xlabel('Time [s]'); ylabel('Spin [Hz]');
title('Spin rate'); grid on;

subplot(2,3,4);
plot(time_log, rad2deg(state_log(:,7)), 'g', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('\phi [deg]');
title('Roll accumulation'); grid on;

subplot(2,3,5);
plot(time_log, rad2deg(state_log(:,8)), 'm', 'LineWidth', 1.5);
hold on;
plot(time_log, rad2deg(state_log(:,9)), 'c--', 'LineWidth', 1.5);
legend('\theta pitch', '\psi yaw');
xlabel('Time [s]'); ylabel('Angle [deg]');
title('Pitch and Yaw'); grid on;

subplot(2,3,6);
plot3(state_log(:,1)/1000, state_log(:,2)/1000, ...
      -state_log(:,3)/1000, 'b', 'LineWidth', 1.5);
xlabel('North [km]'); ylabel('East [km]');
zlabel('Alt [km]');
title('3D Path'); grid on;

