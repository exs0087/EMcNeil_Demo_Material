%% Mars Stay Time for a Human Mars Mission
% A human Mars mission is planned based on the outbound transfer in Problem 3.
% For the return transfer from Mars to Earth, we wish to depart Mars when the
% planetary alignment is such that the Lambert transfer (with a given TOF and
% transfer angle Δν) will succeed. We assume that Mars arrives 
% at periapsis (true anomaly = 0) and that the return transfer is arranged so that 
% the difference between Earth’s heliocentric true anomaly at arrival and Mars’ 
% heliocentric true anomaly at departure equals the desired transfer angle.
%   (a) TOF = 250 days and desired Δν = 170°
%   (b) TOF = 100 days and desired Δν = 80°
%
%
% Author: Erin McNeil
% Date: 03/27/2025

clear; clc;

%% Define fundamental constants
AU = 149597870.7;               % 1 AU in km
mu_sun = 1.32712440018e11;        % Sun's gravitational parameter, km^3/s^2

%% Define Earth orbital parameters (circular)
a_E = 1 * AU;                   % Earth's semi-major axis, km
T_E_days = 365.25;              % Earth's orbital period in days
n_E = 2*pi / (T_E_days*86400);    % Earth's mean motion in rad/s

%% Define Mars orbital parameters (heliocentric)
a_M = 1.524 * AU;               % Mars semi-major axis, km
e_M = 0.0934;                   % Mars eccentricity
i_M = deg2rad(1.85);            % Mars inclination (radians)

% Mars orbital period
T_M_sec = 2*pi*sqrt(a_M^3/mu_sun);
T_M_days = T_M_sec / 86400;     % in days
n_M = 2*pi / T_M_sec;           % Mars mean motion in rad/s

%% Define return transfer parameters for two scenarios

% Scenario 1:
TOF1_days = 250;              % TOF in days
TOF1 = TOF1_days * 86400;       % TOF in seconds
dnu1 = deg2rad(170);            % Desired transfer angle for return (rad)

% Scenario 2:
TOF2_days = 100;              % TOF in days
TOF2 = TOF2_days * 86400;       % TOF in seconds
dnu2 = deg2rad(80);             % Desired transfer angle for return (rad)

%% Compute Earth's heliocentric true anomaly at arrival for each scenario
nu_E_arr1 = mod(n_E * TOF1, 2*pi);  % in rad
nu_E_arr2 = mod(n_E * TOF2, 2*pi);  % in rad

fprintf('--- Return Transfer Phasing ---\n');
fprintf('Scenario 1 (TOF = %d days, Δν = %d°):\n', TOF1_days, 170);
fprintf('  Earth arrival true anomaly: %.2f deg\n', rad2deg(nu_E_arr1));
fprintf('Scenario 2 (TOF = %d days, Δν = %d°):\n', TOF2_days, 80);
fprintf('  Earth arrival true anomaly: %.2f deg\n\n', rad2deg(nu_E_arr2));

%% Determine required Mars departure true anomaly
nu_dep1 = mod(nu_E_arr1 - dnu1, 2*pi);
nu_dep2 = mod(nu_E_arr2 - dnu2, 2*pi);

fprintf('Required Mars departure true anomaly:\n');
fprintf('  Scenario 1: %.2f deg\n', rad2deg(nu_dep1));
fprintf('  Scenario 2: %.2f deg\n\n', rad2deg(nu_dep2));

%% Compute waiting time at Mars (time to go from periapsis, nu=0, to nu_dep)
E1 = 2 * atan( sqrt((1-e_M)/(1+e_M)) * tan(nu_dep1/2) );
if E1 < 0
    E1 = E1 + 2*pi;
end
t_wait1_sec = (E1 - e_M*sin(E1)) / n_M;   % seconds
t_wait1_days = t_wait1_sec / 86400;

E2 = 2 * atan( sqrt((1-e_M)/(1+e_M)) * tan(nu_dep2/2) );
if E2 < 0
    E2 = E2 + 2*pi;
end
t_wait2_sec = (E2 - e_M*sin(E2)) / n_M;   % seconds
t_wait2_days = t_wait2_sec / 86400;

fprintf('Calculated Mars stay (waiting) time:\n');
fprintf('  Scenario 1 (TOF = 250 days, Δν = 170°): %.2f days\n', t_wait1_days);
fprintf('  Scenario 2 (TOF = 100 days, Δν = 80°):  %.2f days\n', t_wait2_days);
