function [rho, temp] = atmosphere(alt)
% Standard lapse-rate atmosphere
% Input: alt = altitude [m] (positive = up)
% Output: rho = density [kg/m^3]
% temp = temperature [K]

T0 = 288.15; % sea level temp [K]
rho0 = 1.225; % sea level density [kg/m^3]
L = 0.0065; % lapse rate [K/m]
R = 287.05; % gas constant [J/kg/K]
g = 9.81; % gravity [m/s^2]

alt = max(alt, 0); % clamp to sea level
temp = T0 - L * alt;
rho = rho0 * (temp / T0)^(g / (L * R));
end