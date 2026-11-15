function [phi_est, pb_est] = roll_spin_extractor(mag_y, time_log, phi_ref_offset)
% STAGE 2: Roll angle and spin rate extraction from magnetometer
% lateral oscillation, via zero-crossing detection.
%
% INPUTS:
%   mag_y          = lateral magnetometer channel (body Y), Nx1
%                    (this is mag_log(:,2) -- pick whichever lateral
%                    channel has the largest oscillation amplitude;
%                    check both Y and Z, use the stronger one or a
%                    combination)
%   time_log       = corresponding time vector, Nx1
%   phi_ref_offset = roll angle [rad] corresponding to an upward
%                    zero-crossing of mag_y (calibration constant --
%                    depends on B_earth direction + mounting; can be
%                    found by comparing to ONE truth sample if doing
%                    a validation test, or determined analytically)
%
% OUTPUTS:
%   phi_est = estimated roll angle [rad], wrapped to [-pi,pi], Nx1
%   pb_est  = estimated spin rate [rad/s], Nx1 (held constant between
%             zero-crossings, i.e. piecewise-constant estimate)
%
% METHOD: find every upward (negative-to-positive) zero-crossing of
% mag_y. Between crossing i and crossing i+1, the body completes
% exactly ONE full spin revolution (since mag_y is a sine-like signal
% at the spin frequency -- one full cycle = one full revolution,
% and a sine crosses zero twice per cycle: once up, once down, so
% consecutive UPWARD crossings are one full period apart).
%
% So: time between two consecutive UPWARD crossings = one full
% spin period T. pb = 2*pi/T over that interval.
% phi at the upward crossing = phi_ref_offset (by definition/calibration).
% Between crossings, phi is interpolated linearly using the
% locally-estimated constant pb.

N = length(mag_y);
phi_est = NaN(N,1);
pb_est  = NaN(N,1);

% --- Find upward zero-crossings ---
crossing_idx  = [];
crossing_time = [];
for i = 2:N
    if mag_y(i-1) < 0 && mag_y(i) >= 0
        % Linear interpolation for sub-sample crossing time
        frac    = -mag_y(i-1) / (mag_y(i) - mag_y(i-1));
        t_cross = time_log(i-1) + frac*(time_log(i)-time_log(i-1));
        crossing_idx(end+1,1)  = i; %#ok<AGROW>
        crossing_time(end+1,1) = t_cross; %#ok<AGROW>
    end
end

if length(crossing_idx) < 2
    warning('roll_spin_extractor: fewer than 2 zero-crossings found.');
    return;
end

% --- Compute pb between consecutive crossings ---
n_cross     = length(crossing_idx);
pb_segments = zeros(n_cross-1,1);
for c = 1:n_cross-1
    T_period = crossing_time(c+1) - crossing_time(c);
    if T_period > 0
        pb_segments(c) = 2*pi / T_period;
    else
        pb_segments(c) = NaN;
    end
end

% --- Fill phi_est, pb_est between each pair of crossings ---
for c = 1:n_cross-1
    idx_start = crossing_idx(c);
    idx_end   = crossing_idx(c+1);
    t_start   = crossing_time(c);
    pb_now    = pb_segments(c);

    for i = idx_start:idx_end
        dt_since_cross = time_log(i) - t_start;
        phi_raw        = phi_ref_offset + pb_now * dt_since_cross;
        phi_est(i)     = atan2(sin(phi_raw), cos(phi_raw)); % wrap
        pb_est(i)      = pb_now;
    end
end

% Fill before first crossing and after last crossing using
% nearest available pb (edge handling)
first_valid = crossing_idx(1);
last_valid  = crossing_idx(end);
if first_valid > 1
    phi_est(1:first_valid-1) = phi_est(first_valid);
    pb_est(1:first_valid-1)  = pb_est(first_valid);
end
if last_valid < N
    phi_est(last_valid+1:end) = phi_est(last_valid);
    pb_est(last_valid+1:end)  = pb_est(last_valid);
end

end