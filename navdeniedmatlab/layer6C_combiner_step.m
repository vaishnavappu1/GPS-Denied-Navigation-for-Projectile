function [x_c, P_c, L] = layer6C_combiner_step(...
    x_d, P_d, Phi_d, x_a, P_a, Phi_a, P_prev_combined)
% LAYER 6C: OPTIMAL COMBINATION OF TWO PREDICTION TRACKS
% Patent US 9,702,674 B2, FIG.6 block 660, FIG.11, equations 41-43
%
% NOTE: this function's math is unchanged and was not the source of the
% bug. The bug was in the CALLING script re-seeding both Track1 and
% Track2's propagation from x_c every step, violating the patent's
% explicit requirement (col. 19) that models 640 and 650 "produce
% INDEPENDENT predictions." See layer6C_main_corrected.m for the fix.
%
% This function still approximates P_da (cross-covariance) as
% Phi_d * P_prev_combined * Phi_a', which is only a rough proxy -- a
% fully rigorous implementation would propagate P_da as its own state
% across time (per the more complete architecture plan discussed).
% Flagged as a known simplification, not corrected in this pass.

cross_da = Phi_d * P_prev_combined * Phi_a';
cross_ad = Phi_a * P_prev_combined * Phi_d';

numerator = P_d - cross_da;
denominator = P_a + P_d - cross_ad - cross_da;

if any(isnan(denominator(:))) || any(isinf(denominator(:)))
    error('layer6C_combiner_step: NaN/Inf in denominator -- check parameters upstream.');
end

% Scale regularization safely across position/velocity partitions
reg_pos = 1e-9 * trace(denominator(1:3,1:3))/3 * eye(3);
reg_vel = 1e-9 * trace(denominator(4:6,4:6))/3 * eye(3);
reg = blkdiag(reg_pos, reg_vel);

L = numerator / (denominator + reg);

% State blend (Patent Eq. 41)
x_c = x_d + L * (x_a - x_d);

% Covariance blend
P_c = P_d - L * (P_d - cross_ad);
P_c = (P_c + P_c') / 2;

% Eigenvalue safety floor to prevent step-to-step collapse
[V_eig, D_eig] = eig(P_c);
d = diag(D_eig);
d(d < 1e-12) = 1e-12;
P_c = V_eig * diag(d) * V_eig';
P_c = (P_c + P_c') / 2;
end