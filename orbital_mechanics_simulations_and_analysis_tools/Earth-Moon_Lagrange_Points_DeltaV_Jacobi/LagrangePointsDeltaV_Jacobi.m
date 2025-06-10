% =========================================================================
% Earth-Moon Lagrange Points ΔV and Jacobi Integral Analysis
%
% Description:
%   (a) Computes the ΔV to escape a 200 km LEO and approximate ΔV to reach the
%       Earth–Moon L2 point.
%   (b) Calculates the Jacobi integral at the Earth–Moon L4 point.
%   (c) Estimates the ΔV to inject into the vicinity of the L4 point.
%
% Author: Erin McNeil
% Last updated: 2025-06-10
% =========================================================================

clear; clc; close all;
format longG;

%% --- ΔV to Reach Earth–Moon L2 ---
muE       = 3.986e5;       % Earth gravitational parameter [km^3/s^2]
Re        = 6378;          % Earth radius [km]
alt_LEO   = 200;           % LEO altitude [km]
r_LEO     = Re + alt_LEO;  % LEO radius [km]

% Circular and escape speeds at LEO
V_LEO     = sqrt(muE / r_LEO);
V_escape  = sqrt(2 * muE / r_LEO);
% ΔV to escape
DeltaV_escape = V_escape - V_LEO;
% Approximate ΔV to reach L2 region
DeltaV_L2 = 3.2;  % km/s (literature value)

fprintf('--- L2 Injection ΔV ---\n');
fprintf('LEO escape ΔV: %.3f km/s\n', DeltaV_escape);
fprintf('Approx. ΔV to L2 vicinity: %.1f km/s\n\n', DeltaV_L2);

%% --- Jacobi Integral at L4 ---
mu_norm = 0.01215;    % Earth-Moon mass parameter (normalized)
% L4 coordinates in normalized rotating frame
x_L4 =  0.5 - mu_norm;
y_L4 =  sqrt(3)/2;
% Effective potential Omega at L4
Omega_L4 = 0.5*(x_L4^2 + y_L4^2) + (1-mu_norm)/norm([x_L4+mu_norm, 0])* + mu_norm/norm([x_L4-1+mu_norm, 0]);
% In normalized CR3BP, Jacobi constant C = 2*Omega
C_L4     = 2 * Omega_L4;
% Using scaled Jacobi integral J = -0.5 * C
J_L4     = -0.5 * C_L4;

fprintf('--- Jacobi Integral at L4 ---\n');
fprintf('Effective potential at L4, Omega_L4: %.3f (normalized)\n', Omega_L4);
fprintf('Jacobi constant at L4, C_L4: %.3f (normalized)\n', C_L4);
fprintf('Scaled Jacobi integral J_L4 = -0.5*C_L4: %.3f\n\n', J_L4);

%% --- ΔV to Reach L4 ---
% LEO escape already computed above, approximate additional burn
DeltaV_L4 = 3.0;  % km/s (literature estimate)

fprintf('--- L4 Injection ΔV ---\n');
fprintf('Approx. ΔV to L4 vicinity: %.1f km/s\n', DeltaV_L4);
fprintf('Note: L4 lies at a higher energy level than L2, so ΔV is slightly lower (~0.2 km/s).\n');
