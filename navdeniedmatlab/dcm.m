function R = dcm(psi, theta)
% Rotation matrix: inertial -> fixed-plane
% Rotates by azimuth psi then elevation theta
% NO roll — fixed-plane has no roll component
% Inputs in RADIANS

R = [ cos(theta)*cos(psi), cos(theta)*sin(psi), -sin(theta);
-sin(psi), cos(psi), 0;
sin(theta)*cos(psi), sin(theta)*sin(psi), cos(theta)];
end