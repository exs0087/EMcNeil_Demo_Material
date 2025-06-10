% =========================================================================
% Venus Flyby Mission Design Analysis
%
% Description:
%   (a) Computes an Earth–Venus transfer using Lambert's problem (universal variable).
%   (b) Determines Venus flyby hyperbolic trajectory parameters:
%       eccentricity, semimajor axis, turning angle (δ), asymptotic anomaly (ν∞).
%   (c) Evaluates post-flyby heliocentric orbit elements to assess whether the
%       gained energy is sufficient to reach Mars or Jupiter.
%
% Author: Erin McNeil
% Last updated: 2025-06-10
% =========================================================================

clear; clc; close all;

%% Constants & Conversion Factors
AU = 1.496e8;                     % km
muSun   = 1.32712440018e11;         % km^3/s^2 (Sun)
muEarth = 3.986e5;                  % km^3/s^2 (Earth)
muVenus = 3.24859e5;                % km^3/s^2 (Venus)

% Earth and Venus orbital radii
Re = 6378;                        % km, Earth radius
LEO_alt = 200;                    % km, LEO altitude
r_LEO = Re + LEO_alt;             % km

rEarth = 1 * AU;                  % km
rVenus = 0.723 * AU;              % km

vEarth = sqrt(muSun / rEarth);      % km/s (~29.78 km/s)
vVenus = sqrt(muSun / rVenus);      % km/s (~35.0 km/s)

% Venus flyby target periapsis (from Venus center)
r_p_flyby = 6500;                 % km

%% Nominal Transfer (Hohmann Baseline)
% For a Hohmann transfer from Earth (1 AU) to Venus (0.723 AU):
a_transfer = (rEarth + rVenus)/2;   
T_nom = pi * sqrt(a_transfer^3 / muSun);    % nominal time-of-flight [s]
T_nom_days = T_nom / (24*3600);               % ~146.6 days

% Nominal transfer speeds using the vis-viva equation:
v_transfer_dep = sqrt(muSun*(2/rEarth - 1/a_transfer));  % km/s at Earth
v_transfer_arr = sqrt(muSun*(2/rVenus  - 1/a_transfer));  % km/s at Venus

%% Earth Departure Injection: Hyperbolic Excess Speed Calculation
% Compute the Earth departure hyperbolic excess speed:
% v_inf_dep = |v_transfer_dep - vEarth|
v_inf_dep_nom = abs(v_transfer_dep - vEarth);  % km/s, expected ~2.5 km/s

% LEO circular speed:
v_LEO_circ = sqrt(muEarth / r_LEO);  % km/s
% Injection burn from LEO:
dV_injection_nom = sqrt(v_inf_dep_nom^2 + v_LEO_circ^2) - v_LEO_circ;  % km/s

%% Venus Arrival Phase Setup
% Assume that at the nominal departure date (July 27, 2026) Earth is at [1,0,0] AU.
% Choose a phase offset for Venus; for example, let Venus be 45° behind Earth.
phi0 = -45 * pi/180;  % in radians

% Venus mean motion (period ~224.7 days):
nVenus = 2*pi / 224.7;  % rad/day
% Nominal arrival phase for Venus:
phi_nom = phi0 + (T_nom_days * nVenus); 

%% Candidate Departure Date Variations
% Allow departure date adjustments of ±2 days from the nominal departure (July 27, 2026).
dep_offsets = -2:2;  % days
nCandidates = length(dep_offsets);

% Preallocate arrays:
v_inf_dep_all = repmat(v_inf_dep_nom, nCandidates, 1);  % constant for Hohmann solution
injection_burns = repmat(dV_injection_nom, nCandidates, 1);
gains = nan(nCandidates,1);
r2_candidates = zeros(3, nCandidates);
v_inf_Venus_all = nan(nCandidates,1);
delta_all = nan(nCandidates,1);

% Use Lambert's method for exact solution:
for i = 1:nCandidates
    offset = dep_offsets(i);
    % For fixed T_nom, arrival time is fixed; only Venus phase shifts.
    phi = phi_nom + (offset * (2*pi/224.7));  % shift in rad
    % Candidate Venus arrival position:
    r2 = rVenus * [cos(phi); sin(phi); 0];
    r2_candidates(:,i) = r2;
    
    % Earth departure position:
    r1 = [1; 0; 0] * AU;  % km
    % Solve Lambert's problem from r1 to r2 over time T_nom using universal variable formulation:
    try
        [v_dep_vec, v_arr_vec] = lambertSolver(r1, r2, T_nom, muSun);
    catch ME
        warning('Candidate %d (offset = %d days) did not converge: %s', i, offset, ME.message);
        continue;
    end
    
    % Compute Earth hyperbolic excess speed at departure:
    v_dep_mag = norm(v_dep_vec);
    v_inf_dep_all(i) = abs(v_dep_mag - vEarth);
    
    % Injection burn from LEO (for reference):
    injection_burns(i) = sqrt(v_inf_dep_all(i)^2 + v_LEO_circ^2) - v_LEO_circ;
    
    % Venus flyby calculations:
    v_Venus_vec = vVenus * [-sin(phi); cos(phi); 0];  % Venus circular velocity vector
    v_inf_Venus_vec = v_arr_vec - v_Venus_vec;         % hyperbolic excess vector at Venus
    v_inf_Venus_all(i) = norm(v_inf_Venus_vec);
    
    % Compute flyby turning angle δ:
    delta = 2 * asin(1 / (1 + (r_p_flyby * v_inf_Venus_all(i)^2) / muVenus));
    delta_all(i) = delta;
    
    % Simulate flyby: Rotate v_inf_Venus_vec by δ about z-axis:
    R = [cos(delta) -sin(delta) 0; sin(delta) cos(delta) 0; 0 0 1];
    v_inf_Venus_out_vec = R * v_inf_Venus_vec;
    
    % Compute post-flyby heliocentric velocity:
    v_after_vec = v_Venus_vec + v_inf_Venus_out_vec;
    gains(i) = norm(v_after_vec) - v_transfer_dep;
end

%% Candidate Selection
% required that the Earth departure hyperbolic excess speed is between 2.5 and 3.6 km/s.
tol_lower = 2.5;
tol_upper = 3.6;
valid_idx = find(v_inf_dep_all >= tol_lower & v_inf_dep_all <= tol_upper & gains > 2.5);
if isempty(valid_idx)
    error('No candidate meets the hyperbolic excess speed and gain requirements.');
end

[~, best_idx_rel] = max(gains(valid_idx));
idx = valid_idx(best_idx_rel);

optimal_offset = dep_offsets(idx);
optimal_phi = phi_nom + (optimal_offset * (2*pi/224.7));
optimal_r2 = r2_candidates(:,idx);
optimal_injection = injection_burns(idx);
optimal_v_inf_dep = v_inf_dep_all(idx);
optimal_gain = gains(idx);

fprintf('=== Problem 3 Part (a) Results ===\n');
fprintf('Optimal departure date: July %d, 2026 (offset %+d days from nominal July 27)\n', 27+optimal_offset, optimal_offset);
fprintf('Nominal flight time: %.2f days\n', T_nom_days);
fprintf('Earth departure hyperbolic excess speed: %.3f km/s\n', optimal_v_inf_dep);
fprintf('Earth departure injection burn from LEO: %.3f km/s\n', optimal_injection);
fprintf('Heliocentric velocity gain from flyby: %.3f km/s\n\n', optimal_gain);

%% Part (b): Venus Flyby Trajectory Elements
% Use optimal candidate's arrival phase:
v_Venus_opt = vVenus * [-sin(optimal_phi); cos(optimal_phi); 0];
v_inf_Venus_vec_opt = v_transfer_arr * [-sin(optimal_phi); cos(optimal_phi); 0] - v_Venus_opt;
v_inf_Venus_opt = norm(v_inf_Venus_vec_opt);
e_flyby = 1 + (r_p_flyby * v_inf_Venus_opt^2) / muVenus;
a_flyby = - muVenus / v_inf_Venus_opt^2;
delta_flyby = 2 * asin(1 / (1 + (r_p_flyby * v_inf_Venus_opt^2) / muVenus));
nu_inf = acos(-1 / e_flyby);

fprintf('=== Problem 3 Part (b): Venus Flyby Trajectory Elements ===\n');
fprintf('Hyperbolic excess speed at Venus: %.3f km/s\n', v_inf_Venus_opt);
fprintf('Flyby eccentricity, e = %.3f\n', e_flyby);
fprintf('Flyby semimajor axis, a = %.3f km\n', a_flyby);
fprintf('Turning angle, δ = %.3f deg\n', rad2deg(delta_flyby));
fprintf('Asymptotic true anomaly, ν∞ = %.3f deg\n\n', rad2deg(nu_inf));

%% Part (c): Post-Flyby Heliocentric Orbit Elements
% Compute post-flyby heliocentric velocity:
R_opt = [cos(delta_flyby) -sin(delta_flyby) 0; sin(delta_flyby) cos(delta_flyby) 0; 0 0 1];
v_inf_Venus_out_opt = R_opt * v_inf_Venus_vec_opt;
v_after_vec = v_Venus_opt + v_inf_Venus_out_opt;
v_after_mag = norm(v_after_vec);

% Use optimal_r2 as the position vector at Venus encounter:
r_arr = norm(optimal_r2);
energy = 0.5 * v_after_mag^2 - muSun/r_arr;
a_post = - muSun / (2 * energy);
h_vec = cross(optimal_r2, v_after_vec);
h_post = norm(h_vec);
e_post = sqrt(1 + 2 * energy * (h_post^2) / (muSun^2));

if e_post < 1
    r_aphelion = a_post * (1 + e_post);
else
    r_aphelion = Inf;
end

fprintf('=== Problem 3 Part (c): Post-Flyby Heliocentric Orbit Elements ===\n');
fprintf('Post-flyby heliocentric speed = %.3f km/s\n', v_after_mag);
fprintf('Specific orbital energy, ε = %.3f km^2/s^2\n', energy);
fprintf('Semimajor axis, a = %.3f km\n', a_post);
fprintf('Eccentricity, e = %.3f\n', e_post);
fprintf('Angular momentum, h = %.3f km^2/s\n', h_post);
if isfinite(r_aphelion)
    fprintf('Aphelion distance = %.3f AU\n', r_aphelion / AU);
else
    fprintf('The post-flyby orbit is hyperbolic (no aphelion).\n');
end

if isfinite(r_aphelion)
    if r_aphelion / AU >= 1.52 && r_aphelion / AU < 5.2
        reachable = 'Mars';
    elseif r_aphelion / AU >= 5.2
        reachable = 'Jupiter (or beyond)';
    else
        reachable = 'only near Earth orbit';
    end
else
    reachable = 'None (hyperbolic exit)';
end
fprintf('Based on the post-flyby orbit, the gained energy is sufficient to reach: %s\n', reachable);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Universal Variable Lambert Solver Function
function [v1, v2] = lambertSolver(r1, r2, T, mu)
    tol = 1e-6;
    max_iter = 500;
    r1_norm = norm(r1);
    r2_norm = norm(r2);
    
    dot_r1r2 = dot(r1, r2);
    theta = acos(min(max(dot_r1r2/(r1_norm*r2_norm), -1), 1));
    % For a short-way transfer; for a long-way, use: theta = 2*pi - theta;
    A = sin(theta) * sqrt((r1_norm * r2_norm) / (1 - cos(theta)));
    if A == 0
        error('Parameter A is zero; no solution exists.');
    end
    
    z = 0;
    ratio = 1;
    iter = 0;
    while (abs(ratio) > tol) && (iter < max_iter)
        [C, S] = stumpff(z);
        y = r1_norm + r2_norm + A*(z*S - 1) / sqrt(C);
        if y < 0
            error('Negative y encountered.');
        end
        x = sqrt(y/C);
        T_z = (x^3 * S + A * sqrt(y)) / sqrt(mu);
        dz = 1e-6;
        [C_plus, S_plus] = stumpff(z + dz);
        y_plus = r1_norm + r2_norm + A*((z+dz)*S_plus - 1) / sqrt(C_plus);
        x_plus = sqrt(y_plus/C_plus);
        T_plus = (x_plus^3 * S_plus + A * sqrt(y_plus)) / sqrt(mu);
        [C_minus, S_minus] = stumpff(z - dz);
        y_minus = r1_norm + r2_norm + A*((z-dz)*S_minus - 1) / sqrt(C_minus);
        x_minus = sqrt(y_minus/C_minus);
        T_minus = (x_minus^3 * S_minus + A * sqrt(y_minus)) / sqrt(mu);
        dTdz = (T_plus - T_minus) / (2 * dz);
        ratio = (T_z - T) / dTdz;
        z = z - ratio;
        iter = iter + 1;
    end
    if iter >= max_iter
        error('Lambert solver did not converge.');
    end
    
    f = 1 - y / r1_norm;
    g = A * sqrt(y / mu);
    g_dot = 1 - y / r2_norm;
    
    v1 = (r2 - f*r1)/g;
    v2 = (g_dot*r2 - r1)/g;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stumpff Functions
function [C, S] = stumpff(z)
    tol_z = 1e-8;
    if z > tol_z
        sqrtz = sqrt(z);
        C = (1 - cos(sqrtz))/z;
        S = (sqrtz - sin(sqrtz))/(sqrtz^3);
    elseif z < -tol_z
        sqrt_negz = sqrt(-z);
        C = (1 - cosh(sqrt_negz))/z;
        S = (sinh(sqrt_negz) - sqrt_negz)/(sqrt_negz^3);
    else
        C = 0.5;
        S = 1/6;
    end
end
