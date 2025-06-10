% =========================================================================
% Earth-Pointing Satellite Attitude Dynamics Simulation
%
% Description:
%   Simulates the desired quaternion attitude (qc) and angular velocity (wc)
%   & acceleration (wcd) profiles for a satellite in a highly elliptical orbit.
%   Outputs are saved as qc (4×N), wc (3×N), wcd (3×N) for
%   downstream scripts.
%
% Author: Erin McNeil
% Last updated: 2024-11-16
% =========================================================================

clear; close all; clc;
format long;

%% --- Inertia & Constants ---
Jx = 100; Jy = 90; Jz = 120;             % kg·m^2
J  = diag([Jx, Jy, Jz]);

mu = 3.986e5;                            % km^3/s^2
Re = 6378;                               % km

%% --- Orbital Elements (HEO Example) ---
a     = 26559;                          % km
e     = 0.704482;
i     = deg2rad(63.171);
Omega = deg2rad(206.346);
omega = deg2rad(281.646);
M0    = deg2rad(12.998);
t0    = 0;
n     = sqrt(mu / a^3);                 % rad/s

%% --- Time Setup ---
t_final = 64800;                        % seconds (18 h)
dt      = 1;                            % s
T       = 0:dt:t_final;                 
N       = numel(T);                    

%% --- Preallocate Outputs ---
qc  = zeros(4,  N);    % quaternion [x; y; z; w]
wc  = zeros(3,  N);    % angular velocity [rad/s]
wcd = zeros(3,  N);    % angular acceleration [rad/s^2]

%% --- Propagation & Attitude Loop ---
for k = 1:N
    t = T(k);
    M = M0 + n*(t - t0);            % mean anomaly
    
    % Solve Kepler's equation for E
    syms E_sym
    E = double(vpasolve(E_sym - e*sin(E_sym) == M, E_sym));
    
    % Orbital-plane position & velocity
    r = a*(1 - e*cos(E));
    x = a*(cos(E)-e);  y = a*sqrt(1-e^2)*sin(E);
    xdot = - (sqrt(mu/a)/(1-e*cos(E))) * sin(E);
    ydot =   (sqrt(mu/a)/(1-e*cos(E))) * sqrt(1-e^2) * cos(E);
    
    % Inertial transform
    R1 = [ cos(Omega) -sin(Omega) 0; sin(Omega) cos(Omega) 0; 0 0 1];
    R2 = [1 0 0; 0 cos(i) -sin(i); 0 sin(i) cos(i)];
    R3 = [ cos(omega) -sin(omega) 0; sin(omega) cos(omega) 0; 0 0 1];
    Q  = R1*R2*R3;
    
    r_comp  = Q*[x; y; 0];
    v_comp  = Q*[xdot; ydot; 0];
    rdotdot = [0;0;0];  % central gravity omitted in LVLH xdotdot
    
    % Build LVLH frame
    z_LVLH = -r_comp/norm(r_comp);
    y_LVLH = -cross(r_comp, v_comp)/norm(cross(r_comp, v_comp));
    x_LVLH = cross(y_LVLH, z_LVLH);
    A_IO    = [x_LVLH, y_LVLH, z_LVLH];
    
    % Quaternion [x; y; z; w]
    theta = acos((trace(A_IO)-1)/2);
    e_axis = (1/(2*sin(theta))) * [...
        A_IO(3,2)-A_IO(2,3);
        A_IO(1,3)-A_IO(3,1);
        A_IO(2,1)-A_IO(1,2)];
    q_cur = [ e_axis*sin(theta/2); cos(theta/2) ];
    qc(:,k) = q_cur;
    
    % Angular velocity in LVLH
    w_cur = [...
        0;
        -norm(cross(r_comp, v_comp))/norm(r_comp)^2;
        (norm(r_comp)*dot(y_LVLH, rdotdot)) / norm(cross(r_comp, v_comp)) ];
    wc(:,k) = w_cur;
    
    % Angular acceleration (no external torque)
    wcd(:,k) = J \ ( -cross(w_cur, J*w_cur) );
end

%% --- Save Outputs ---
save('EarthPointingData.mat', 'qc', 'wc', 'wcd', '-v7.3');

% --- Plotting Results ---
% Altitude vs Time
figure;
plot(T/3600, alt_km, 'LineWidth', 1.5);
xlabel('Time [hours]');
ylabel('Altitude [km]');
title('Satellite Altitude Over 18-hour Orbit');
saveas(gcf, 'altitude_vs_time.png');

% Quaternion Components vs Time
figure;
subplot(4,1,1); plot(T/3600, quat_desired(:,1)); ylabel('q_1');
subplot(4,1,2); plot(T/3600, quat_desired(:,2)); ylabel('q_2');
subplot(4,1,3); plot(T/3600, quat_desired(:,3)); ylabel('q_3');
subplot(4,1,4); plot(T/3600, quat_desired(:,4)); ylabel('q_4');
sgtitle('Desired Attitude Quaternion Trajectory');

% Angular Velocity vs Time
figure;
subplot(3,1,1); plot(T/3600, omega_desired(:,1)); ylabel('ω_1 [rad/s]');
subplot(3,1,2); plot(T/3600, omega_desired(:,2)); ylabel('ω_2 [rad/s]');
subplot(3,1,3); plot(T/3600, omega_desired(:,3)); ylabel('ω_3 [rad/s]');
sgtitle('Desired Angular Velocity in LVLH Frame');
