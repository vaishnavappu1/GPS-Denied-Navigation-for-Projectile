function p = projectile_params()
% All parameters for 155mm M795 Artillery Shell
% Patent US 9,702,674 B2

% MASS PROPERTIES
p.m = 43.2; % mass [kg]
p.D = 0.155; % diameter [m]
p.A = pi/4 * p.D^2; % reference area [m^2]
p.d = 0.155; % reference length [m]

% MOMENTS OF INERTIA (McCoy Modern Exterior Ballistics)
p.Ix = 0.1658; % axial [kg.m^2]
p.Iy = 1.724; % transverse [kg.m^2]

% AERODYNAMIC FORCE COEFFICIENTS
p.CD0 = 0.2; % zero-AOA drag
p.CX_a2 = 3.44; % AOA-squared drag
p.CNa = 2.48; % normal force slope [/rad]

% SPIN AND MAGNUS COEFFICIENTS
p.p0 = 50 * 2 * pi; % initial spin [rad/s]
p.Clp = -0.0192; % roll damping
p.CYpa = 0.0172; % Magnus force coefficient
p.Cnpa = 0.0126; % Magnus moment coefficient

% PITCH MOMENT COEFFICIENTS
p.Cma = -0.50; % pitch moment slope [/rad]
p.Cmq = -5.0; % pitch damping [/rad]

% CANARD (zeroed until guidance layer)
p.CX0c = 0.05;
p.CZc = 0.05;
p.CYac = 0.10;
p.CG = 0.55;
p.CPz = 0.10;
p.CPy = 0.10;

% PHYSICAL CONSTANTS
p.g = 9.81; % gravity [m/s^2]

% LAUNCH CONDITIONS
p.V0 = 827; % muzzle velocity [m/s]
p.elev = 45; % elevation [deg]
p.az = 0; % azimuth [deg]

% SIMULATION
p.dt = 0.001; % timestep [s]
p.t_max = 120; % max time [s]

% SENSOR PARAMETERS
p.B_earth = [28000; 900; -37000]; % Earth field [nT]
p.sigma_mag = 500; % mag noise [nT]
p.bias_mag = [200; -150; 300]; % mag bias [nT]
p.sigma_accel = 0.05; % accel noise [g]
p.bias_accel = [0.02; -0.015; 0.01]; % accel bias [g]
p.epsilon_c = 0.002; % sensor offset [m]
p.sigma_gps_pos = 5.0; % GPS pos noise [m]
p.sigma_gps_vel = 0.1; % GPS vel noise [m/s]
p.t_gps_end = 20.0; % GPS denied after [s]

% Accelerometer offset parameters
p.accel_offset_j = 0.002;    % [m] j-axis sensor offset
p.accel_offset_k = 0.0015;   % [m] k-axis sensor offset

% Balloting motion (muzzle exit disturbance)
p.balloting_j    = 50.0;     % [m/s^2] initial amplitude
p.balloting_k    = 40.0;     % [m/s^2]
p.balloting_tau  =  5.0;     % [s] decay time constant

% Magnetometer calibration matrix
p.S_mag = eye(3);            % ideal calibration
end
