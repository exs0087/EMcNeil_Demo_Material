%% Earth-to-Mars Interplanetary Transfer
% This script computes the heliocentric transfer orbit for an interplanetary
% probe departing Earth en route to Mars. The transfer has a specified
% transfer angle (Δν = 170deg) and time-of-flight (TOF = 250 days). It is desired 
% that the probe arrive in Mars’ heliocentric orbit when Mars is at its periapsis.
%   (a) Compute the heliocentric departure (V1) and arrival (V2) velocities and 
%       the transfer orbit’s classical elements.
%   (b) Compute the Earth–Mars phase angle required at departure.
%   (c) Compute the synodic period for the geometric alignment.
%
% Author: Erin McNeil
% Last updated: 2025-06-10

clear; clc;

%% Constants and Unit Conversions
AU = 1.496e8;       % 1 AU in km
deg = pi/180;
mu_sun = 1.32712440018e11; % Sun's gravitational parameter, km^3/s^2

%% Given Transfer Parameters
TOF_days = 250;         % time-of-flight in days
TOF = TOF_days * 86400; % convert to seconds
dnu = 170 * deg;        % transfer angle in radians

%% Define Earth and Mars orbits
% Assume Earth's orbit is circular at 1 AU.
a_E = AU;              % km
e_E = 0;               % circular
% Earth’s position at departure on the ecliptic:
r1 = [a_E*cos(dnu); a_E*sin(dnu); 0];

% Mars orbit parameters
a_M_AU = 1.524;                    % semi-major axis in AU
a_M = a_M_AU * AU;                 % km
e_M = 0.0934;
i_M = 1.85 * deg;                  % inclination
% Assuming that Mars’ longitude of periapsis is 0 deg, so its periapsis lies along the x-axis.
rp_M = a_M*(1 - e_M);              % Mars periapsis distance in km
r2 = [rp_M; 0; 0];                 % set Mars at periapsis 

%% Solve Lambert's Problem 
[V1, V2, tof_out] = lambert_universal(r1, r2, TOF, mu_sun, 0);

% Check that the computed time-of-flight matches the input TOF
if abs(tof_out - TOF) > 1e-3
    warning('Computed TOF does not match input TOF.');
end

% Compute transfer orbit classical orbital elements from the departure state (r1, V1)
[coe, ~] = rv2coe(r1, V1, mu_sun);
a_transfer = coe(1);
e_transfer = coe(2);
i_transfer = coe(3)/deg;         % in degrees
RAAN_transfer = coe(4)/deg;        % in degrees
argp_transfer = coe(5)/deg;        % in degrees
nu_departure = coe(6)/deg;         % in degrees

fprintf('--- Transfer Orbit ---\n');
fprintf('Departure Velocity (V1): [%.6f, %.6f, %.6f] km/s\n', V1);
fprintf('Arrival   Velocity (V2): [%.6f, %.6f, %.6f] km/s\n', V2);
fprintf('Transfer Orbit Elements:\n');
fprintf('  Semi-major axis, a = %.6f km\n', a_transfer);
fprintf('  Eccentricity, e    = %.6f\n', e_transfer);
fprintf('  Inclination, i     = %.6f deg\n', i_transfer);
fprintf('  RAAN               = %.6f deg\n', RAAN_transfer);
fprintf('  Arg of Periapsis   = %.6f deg\n', argp_transfer);
fprintf('  True Anomaly (dep) = %.6f deg\n', nu_departure);

%% Compute Required Earth-Mars Phase Angle
T_M = 686.98 * 86400;   % Mars period in seconds
n_M = 2*pi / T_M;       % mean motion of Mars (rad/s)

phi_required = mod(360 - (n_M*TOF * (180/pi)), 360);
fprintf('\n--- Earth-Mars Phase Angle ---\n');
fprintf('Required Earth-Mars phase angle at departure: %.6f deg\n', phi_required);

%% Compute Synodic Period (Part c)
T_E = 365.25 * 86400;  % Earth period in seconds
n_E = 2*pi / T_E;      % mean motion of Earth (rad/s)
n_syn = abs(n_E - n_M);
T_syn = (2*pi / n_syn) / 86400; % in days
fprintf('\n--- Synodic Period ---\n');
fprintf('Synodic period: %.6f days\n', T_syn);


%% Local Function: lambert_universal
function [V1, V2, tof_out] = lambert_universal(r1, r2, tof, mu, dm)   
    r1_norm = norm(r1);
    r2_norm = norm(r2);
    cos_dnu = dot(r1, r2) / (r1_norm*r2_norm);
    cos_dnu = min(max(cos_dnu, -1), 1);
    dnu = acos(cos_dnu);
    if dm == 1
        dnu = 2*pi - dnu;
    end
    A = sin(dnu)*sqrt(r1_norm*r2_norm/(1 - cos_dnu));
    
    % Initialize iteration for universal variable z
    z = 0;
    tol_z = 1e-8;
    ratio = 1;
    count = 0;
    max_iter = 1000;
    while abs(ratio) > tol_z && count < max_iter
        [C, S] = stumpff(z);
        y = r1_norm + r2_norm + A*(z*S - 1)/sqrt(C);
        F = (y/C)^(3/2)*S + A*sqrt(y) - sqrt(mu)*tof;
        % Derivative dF/dz 
        dz = 1e-6;
        [C_p, S_p] = stumpff(z+dz);
        y_p = r1_norm + r2_norm + A*((z+dz)*S_p - 1)/sqrt(C_p);
        F_p = (y_p/C_p)^(3/2)*S_p + A*sqrt(y_p) - sqrt(mu)*tof;
        dFdz = (F_p - F)/dz;
        ratio = F / dFdz;
        z = z - ratio;
        count = count + 1;
    end
    if count >= max_iter
        error('Lambert solver did not converge');
    end
    [C, S] = stumpff(z);
    y = r1_norm + r2_norm + A*(z*S - 1)/sqrt(C);
    f = 1 - y/r1_norm;
    g = A * sqrt(y/mu);
    V1 = (r2 - f*r1)/g;
    g_dot = 1 - y/r2_norm;
    V2 = (g_dot*r2 - r1)/g;
    tof_out = ((y/C)^(3/2)*S + A*sqrt(y)) / sqrt(mu);
end

function [C, S] = stumpff(z)
    if z > 0
        sz = sqrt(z);
        C = (1 - cos(sz)) / z;
        S = (sz - sin(sz)) / (sz^3);
    elseif z < 0
        sz = sqrt(-z);
        C = (1 - cosh(sz)) / z;
        S = (sinh(sz) - sz) / (sz^3);
    else
        C = 1/2;
        S = 1/6;
    end
end

%% Local Function: rv2coe (State Vector to Orbital Elements)
function [coe, r_v] = rv2coe(r, v, mu)
    % Inputs:
    %   r, v: 3x1 state vectors in km and km/s
    %   mu: gravitational parameter in km^3/s^2
    % Outputs:
    %   coe = [a, e, i, RAAN, argp, nu]
    
    r_norm = norm(r);
    v_norm = norm(v);
    h_vec = cross(r, v);
    h_norm = norm(h_vec);
    i = acos(abs(h_vec(3)/h_norm));
    
    % Node vector
    K = [0; 0; 1];
    n_vec = cross(K, h_vec);
    n_norm = norm(n_vec);
    
    % Eccentricity vector
    e_vec = (1/mu)*((v_norm^2 - mu/r_norm)*r - dot(r,v)*v);
    e = norm(e_vec);
    
    a = 1 / (2/r_norm - v_norm^2/mu);
    
    if n_norm > 1e-8
        RAAN = acos(n_vec(1)/n_norm);
        if n_vec(2) < 0
            RAAN = 2*pi - RAAN;
        end
    else
        RAAN = 0;
    end
    
    if n_norm > 1e-8 && e > 1e-8
        argp = acos(dot(n_vec, e_vec)/(n_norm*e));
        if e_vec(3) < 0
            argp = 2*pi - argp;
        end
    else
        argp = 0;
    end
    
    if e > 1e-8
        nu = acos(dot(e_vec, r)/(e*r_norm));
        if dot(r, v) < 0
            nu = 2*pi - nu;
        end
    else
        nu = 0;
    end
    coe = [a, e, i, RAAN, argp, nu];
    r_v = [r; v];
end
