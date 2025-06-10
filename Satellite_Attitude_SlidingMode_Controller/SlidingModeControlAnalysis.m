% =========================================================================
% Sliding Mode Attitude Control Performance Analysis
%
% Description:
%   Compares sliding mode controllers (switching vs. saturation) by simulating
%   attitude response and computing quaternion and angular rate errors relative 
%   to reference Earth-pointing trajectories.
%
% Outputs:
%   - SlidingModeErrors.mat : struct with error metrics for both control modes
%   - Comparison plots saved as PNG files
%
% Author: Erin McNeil
% Last updated: 2024-11-16
% =========================================================================

clear; close all; clc;

%% --- Load Reference Trajectory Data ---
data = load('EarthPointingData.mat', 'qc', 'wc');
qc_ref = data.qc;    % 4×N reference quaternion trajectory
wc_ref = data.wc;    % 3×N reference angular velocity trajectory
N = size(qc_ref, 2);
T = 0:(N-1);         % Time vector [s]

%% --- Simulation Setup ---
y0 = [zeros(3,1); [0; 0; 0; 1]];  % Initial state: [angular velocity; quaternion]

control_modes = {'Switching', 'Saturation'};
mode_codes = [1, 2];

results = struct('mode', {}, 'theta_err', {}, 'omega_err', {}, ...
                 'theta_err_total', {}, 'omega_err_total', {});

%% --- Simulate and Compute Errors for Each Mode ---
for idx = 1:length(mode_codes)
    cmode = mode_codes(idx);
    mode_name = control_modes{idx};

    % Simulate attitude with custom ODE integrator
    ysim = ode4n(@eulerseqns, T, y0, cmode);  % Output size: (N)×7

    % Preallocate error arrays
    theta_err = zeros(1, N);
    omega_err = zeros(1, N);

    for k = 1:N
        % Extract simulated quaternion and angular velocity
        q_sim = ysim(k, 4:7).';
        w_sim = ysim(k, 1:3).';

        % Compute quaternion error angle [deg]
        Xi = [ qc_ref(4,k), -qc_ref(3,k),  qc_ref(2,k);
               qc_ref(3,k),  qc_ref(4,k), -qc_ref(1,k);
              -qc_ref(2,k),  qc_ref(1,k),  qc_ref(4,k);
              -qc_ref(1,k), -qc_ref(2,k), -qc_ref(3,k)];
        dq = Xi' * q_sim;
        theta_err(k) = 2 * asin(norm(dq)) * (180/pi);

        % Compute angular rate error [deg/s]
        omega_err(k) = norm(w_sim - wc_ref(:,k)) * (180/pi);
    end

    % Store results
    results(idx).mode = mode_name;
    results(idx).theta_err = theta_err;
    results(idx).omega_err = omega_err;
    results(idx).theta_err_total = sum(theta_err);
    results(idx).omega_err_total = sum(omega_err);
end

%% --- Plot Comparison Results ---

figure;
plot(T/3600, results(1).theta_err, 'LineWidth', 1.2); hold on;
plot(T/3600, results(2).theta_err, 'LineWidth', 1.2); hold off;
legend(control_modes, 'Location', 'best');
xlabel('Time [hours]');
ylabel('Quaternion Error [deg]');
title('Quaternion Tracking Error: Switching vs. Saturation');
saveas(gcf, 'QuaternionErrorComparison.png');

figure;
plot(T/3600, results(1).omega_err, 'LineWidth', 1.2); hold on;
plot(T/3600, results(2).omega_err, 'LineWidth', 1.2); hold off;
legend(control_modes, 'Location', 'best');
xlabel('Time [hours]');
ylabel('Angular Rate Error [deg/s]');
title('Angular Rate Tracking Error: Switching vs. Saturation');
saveas(gcf, 'AngularRateErrorComparison.png');

%% --- Save Error Metrics ---
save('SlidingModeErrors.mat', 'results', '-v7.3');
