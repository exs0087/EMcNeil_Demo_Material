% =========================================================================
% Universal-Variable Lambert (Gauss) Solver
%
% Description:
%   Demonstrates the universal-variable formulation for Lambert's problem.
%   Iteratively solves for the universal anomaly z using Newton-Raphson with
%   adaptive damping, computes Lagrange coefficients (f, g, g_dot), initial
%   and final velocity vectors, and extracts orbital elements from the initial state.
%
% Usage:
%   Run the script to execute a suite of test cases. Results and diagnostics
%   are displayed in the Command Window.
%
% Author: Erin McNeil
% Last updated: 2025-06-10
% =========================================================================

clear; clc;

%% Define constants
mu = 1;          % Gravitational parameter [DU^3/TU^2]
alpha = 0.05;    % Damping factor for Newton-Raphson iteration
tol = 1e-6;      % Convergence tolerance
max_iter = 1000; % Maximum iterations allowed

%% Define test cases as a cell array
% Each row: { r1 (DU), r2 (DU), T (TU), transfer type ('short' or 'long') }
cases = { ...
    [0.5, 0.6, 0.7],  [0.0, -1.0, 0.0], 20,    'long';   % Case (a)
    [0.3, 0.7, 0.4],  [0.6, -1.4, 0.8], 5,     'short';  % Case (b)
    [0.5, 0.6, 0.7],  [0.0, 1.0, 0.0],  1.2,   'long';   % Case (c)
    [-0.2, 0.6, 0.3], [0.4, 1.2, 0.6],  50,    'short';  % Case (d)
    [1.0, 0.0, 0.0],  [0.0, 1.0, 0.0],  0.0001,'short';  % Case (e)
    [-0.4, 0.6, -1.2],[0.2, -0.3, 0.6], 5,     'short'   % Case (f)
};

%% Loop over test cases
for k = 1:size(cases,1)
    % Extract case data
    r1 = cases{k,1}(:);
    r2 = cases{k,2}(:);
    T_target = cases{k,3};
    transfer_type = cases{k,4};  % 'short' or 'long'
    
    fprintf('==================== Case %d ====================\n', k);
    fprintf('r1 = [%.4f, %.4f, %.4f] DU\n', r1);
    fprintf('r2 = [%.4f, %.4f, %.4f] DU\n', r2);
    fprintf('Time of Flight (T) = %.6f TU, Transfer = %s-way\n', T_target, transfer_type);
    
    r1_mag = norm(r1);
    r2_mag = norm(r2);
    
    %% Compute the transfer angle, theta
    cos_theta = dot(r1, r2) / (r1_mag*r2_mag);
    % Clamp the value to [-1,1] 
    theta = acos( min(max(cos_theta, -1), 1) );
    if strcmpi(transfer_type, 'long')
        theta = 2*pi - theta;
    end
    fprintf('Transfer angle, theta = %.6f rad\n', theta);
    
    %% Compute parameter A
    A = sin(theta)*sqrt( (r1_mag*r2_mag)/(1 - cos(theta)) );
    if A == 0
        error('Parameter A is zero; no solution exists for this configuration.');
    end
    
    %% Choose initial guess for z
    if T_target > 10
        z = 2.5;  
    else
        z = 0;
    end

    %% Newton-Raphson iteration to solve T(z) = T_target
    iter = 0;
    converged = false;
    
    fprintf('\nIteration history:\n');
    fprintf('Iter\t z\t\t y\t\t x\t\t t(z)\t\t dt/dz\t\t z_next\n');
    
    while ~converged && iter < max_iter
        % Compute Stumpff functions for current z
        [C, S] = stumpff(z);
        
        % Compute y
        y = r1_mag + r2_mag + A*(z*S - 1)/sqrt(C);
        if y < 0
            error('Encountered negative y value (invalid solution).');
        end
        
        % Compute x and t(z)
        x = sqrt(y/C);
        t_z = (x^3 * S + A*sqrt(y)) / sqrt(mu);
        
        % Approximate dt/dz via finite differences
        dz = 1e-6;
        z_plus = z + dz;
        z_minus = z - dz;
        
        [C_plus, S_plus] = stumpff(z_plus);
        y_plus = r1_mag + r2_mag + A*(z_plus*S_plus - 1)/sqrt(C_plus);
        if y_plus < 0, y_plus = y; end  % safeguard
        x_plus = sqrt(y_plus/C_plus);
        t_plus = (x_plus^3 * S_plus + A*sqrt(y_plus)) / sqrt(mu);
        
        [C_minus, S_minus] = stumpff(z_minus);
        y_minus = r1_mag + r2_mag + A*(z_minus*S_minus - 1)/sqrt(C_minus);
        if y_minus < 0, y_minus = y; end
        x_minus = sqrt(y_minus/C_minus);
        t_minus = (x_minus^3 * S_minus + A*sqrt(y_minus)) / sqrt(mu);
        
        dt_dz = (t_plus - t_minus) / (2*dz);
        
        % Newton-Raphson update with damping factor alpha
        z_next = z + alpha * (T_target - t_z) / dt_dz;
        
        % Dynamically adjust the damping factor
        delta_z = abs(z_next - z);
        if delta_z > 1e-1
            alpha = 0.01;
        elseif delta_z > 1e-3
            alpha = 0.1;
        else
            alpha = 0.4;
        end
        
        fprintf('%d\t %.8f\t %.8f\t %.8f\t %.8f\t %.8f\t %.8f\n', iter, z, y, x, t_z, dt_dz, z_next);
        
        if abs(z_next - z) < tol
            converged = true;
        end
        
        z = z_next;
        iter = iter + 1;
    end
    
    if ~converged
        error('Newton-Raphson did not converge for case %d', k);
    end
    fprintf('Converged after %d iterations.\n', iter);
    
    %% Compute Lagrange coefficients:
    f = 1 - y / r1_mag;
    g = A * sqrt(y / mu);
    g_dot = 1 - y / r2_mag;
    
    %% Compute velocity vectors:
    v1 = (r2 - f*r1) / g;
    v2 = (g_dot*r2 - r1) / g;
    
    % Compute energy and angular momentum as a check
    energy = norm(v1)^2/2 - mu/r1_mag;
    h_vec = cross(r1, v1);
    h = norm(h_vec);
    fprintf('Diagnostic: Energy = %.6f DU^2/TU^2, |h| = %.6f DU^2/TU\n', energy, h);
    
    %% Compute orbital elements from (r1,v1)
    e_vec = (1/mu) * ((norm(v1)^2 - mu/r1_mag)*r1 - (dot(r1, v1))*v1);
    e = norm(e_vec);
    
    if abs(e - 1) < 1e-3
        fprintf('Case %d: Orbit is nearly parabolic (e=%.6f); some orbital elements may be undefined.', k, e);
        a = -mu/(2*energy);
        i = acos(h_vec(3)/h);
        RAAN = NaN;
        arg_perigee = NaN;
        true_anomaly = acos(dot(e_vec, r1)/(e*r1_mag));
    else
        if abs(energy) > 1e-10
            a = -mu / (2 * energy);
        else
            a = Inf;  % parabolic/hyperbolic borderline
        end
        
        i = acos(h_vec(3) / h);
        
        K = [0;0;1];
        n_vec = cross(K, h_vec);
        n = norm(n_vec);
        
        if n > 1e-8
            RAAN = acos(n_vec(1)/n);
            if n_vec(2) < 0
                RAAN = 2*pi - RAAN;
            end
        else
            RAAN = NaN;
        end
        
        if e > 1e-8 && n > 1e-8
            arg_perigee = acos(dot(n_vec, e_vec)/(n*e));
            if e_vec(3) < 0
                arg_perigee = 2*pi - arg_perigee;
            end
        else
            arg_perigee = NaN;
        end
        
        if e > 1e-8
            true_anomaly = acos(dot(e_vec, r1)/(e*r1_mag));
            if dot(r1,v1) < 0
                true_anomaly = 2*pi - true_anomaly;
            end
        else
            true_anomaly = NaN;
        end
    end
    
    %% Determine orbit type based on eccentricity for fun
    tol_e = 1e-3;  % tolerance for determining parabolic
    if abs(e - 1) < tol_e
        orbit_type = 'parabolic';
    elseif e < 1
        orbit_type = 'elliptical';
    else
        orbit_type = 'hyperbolic';
    end


    %% Display results
    fprintf('\nComputed Velocities:\n');
    fprintf('v1 = [%.8f, %.8f, %.8f] DU/TU\n', v1);
    fprintf('v2 = [%.8f, %.8f, %.8f] DU/TU\n', v2);
    
    fprintf('\nOrbital Elements from (r1,v1):\n');
    fprintf('Semi-major axis, a = %.8f DU\n', a);
    fprintf('Eccentricity, e = %.8f\n', e);
    fprintf('Inclination, i = %.8f rad\n', i);
    fprintf('RAAN = %.8f rad\n', RAAN);
    fprintf('Argument of perigee = %.8f rad\n', arg_perigee);
    fprintf('True anomaly = %.8f rad\n', true_anomaly);
    fprintf('Orbit type: %s\n', orbit_type);
    fprintf('------------------------------------------------------\n\n');
end

%% Local function: Stumpff functions
function [C, S] = stumpff(z)
    % Computes the Stumpff functions C(z) and S(z)
    tol_z = 1e-8;
    if z > tol_z
        sqrtz = sqrt(z);
        C = (1 - cos(sqrtz)) / z;
        S = (sqrtz - sin(sqrtz)) / (sqrtz^3);
    elseif z < -tol_z
        sqrt_minus_z = sqrt(-z);
        C = (1 - cosh(sqrt_minus_z)) / z;
        S = (sinh(sqrt_minus_z) - sqrt_minus_z) / (sqrt_minus_z^3);
    else
        % When z ~ 0 use the series expansion 
        C = 1/2;
        S = 1/6;
    end
end
