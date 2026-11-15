function [R_i2S, chi] = build_S_frame(vel_est, B_earth)
% BUILD_S_FRAME — constructs the patent's steady-state "-S frame"
% (US 9,702,674 B2, col. 12, FIG. 7, block 740/750 region).

V = norm(vel_est);
if V < 1; V = 1; end

psi_v   = atan2(vel_est(2), vel_est(1));
theta_v = atan2(-vel_est(3), sqrt(vel_est(1)^2 + vel_est(2)^2));
R_i2f   = dcm(psi_v, theta_v);

B_f = R_i2f * B_earth;

% --- FIX 1: Correct the vector component ordering for roll extraction ---
% This aligns the cross-axis angle extraction with the actual fixed-plane reference horizon.
chi = atan2(B_f(2), B_f(3)); 

% --- FIX 2: Implement clean forward DCM rotation matching standard aerospace rules ---
% Avoid passing double-negatives into the sine entries. 
chi = atan2(B_f(2), B_f(3));
R_roll_chi = [1, 0,        0;
              0, cos(chi), -sin(chi);
              0, sin(chi),  cos(chi)];
R_i2S = R_roll_chi * R_i2f;


end