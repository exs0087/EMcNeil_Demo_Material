% =========================================================================
% Aerobraking Simulation for Mars Reconnaissance Orbiter (MRO)
%
% Description:
%   Simulates three atmospheric drag passes for the MRO using an exponential
%   Mars atmosphere model. Updates orbital elements (a, e, M) via Gauss' equations 
%   during drag-on phases and propagates in vacuum between passes, detecting entry 
%   and exit at a specified atmospheric threshold.
%
% Author: Erin McNeil
% Contact: your.email@domain.com
% Last updated: 2025-06-09
% =========================================================================

clear; clc;
format longG;

%% --- Mars and Spacecraft Parameters ---
mu_mars    = 4.282837e4;    % Mars gravitational parameter [km^3/s^2]
R_mars     = 3396.2;        % Mars equatorial radius [km]
A_ref      = 37.5e-6;       % Spacecraft reference area [km^2]
m_sc       = 2180;          % Spacecraft mass [kg]

%% --- Atmospheric Model ---
rho0_surface = 0.020 * 1e9;           % Surface density [kg/km^3]
H_scale      = 10.8 / 1e3;            % Scale height [km]
rho_thresh   = 1e-10 * 1e9;           % Density threshold [kg/km^3]

% Compute threshold radius for event detection
h_thresh  = -H_scale * log(rho_thresh / rho0_surface);
r_thresh  = R_mars + h_thresh;
fprintf('Atmospheric threshold at r = %.2f km\n', r_thresh);

%% --- Initial Orbit Elements at Entry ---
a0   = 25685.77;       % Semi-major axis [km]
e0   = 0.8619;         % Eccentricity
p0   = a0 * (1 - e0^2);

% Find true anomaly at entry (r = r_thresh) on inbound leg
cos_nu = (p0/r_thresh - 1) / e0;
nu_entry = -acos(max(min(cos_nu,1),-1));

% Convert to mean anomaly
E_entry = 2 * atan2(sqrt(1-e0)*sin(nu_entry/2), sqrt(1+e0)*cos(nu_entry/2));
M_entry = E_entry - e0*sin(E_entry);

Y_entry = [a0; e0; M_entry];  % [a; e; M] at drag-pass entry

%% --- Simulation Settings ---
dt       = 1;      % Time step [s]
T_max    = 1e5;    % Max time per phase [s]
num_pass = 3;      % Number of drag passes

% History storage
a_hist = zeros(1, num_pass);
e_hist = zeros(1, num_pass);

%% --- Main Aerobraking Loop ---
for pass = 1:num_pass
    % Drag-on phase: integrate until exit (r >= r_thresh)
    [~, ~, Y_exit] = integrate_until_event(Y_entry, dt, T_max, true,  r_thresh);
    a_new = Y_exit(1);
    e_new = Y_exit(2);
    fprintf('\nAfter drag pass %d: a = %.6f km, e = %.6f\n', pass, a_new, e_new);
    a_hist(pass) = a_new;
    e_hist(pass) = e_new;
    
    % Vacuum phase: integrate until next entry (r <= r_thresh)
    [~, ~, Y_entry] = integrate_until_event(Y_exit, dt, T_max, false, r_thresh);
end

%% --- Final Summary ---
fprintf('\nInitial a: %.6f km, Final a: %.6f km\n', a0, a_hist(end));
orbit_type = ternary(e_hist(end)<1, 'Elliptical', 'Hyperbolic or Parabolic');
fprintf('Final orbit type: %s\n', orbit_type);

%% --- Save Results ---
save('MRO_AerobrakingResults.mat', 'a_hist', 'e_hist', '-v7.3');

%% --- Nested Functions ---

function dY = orb_elem_deriv(~, Y, drag_on, r_thresh)
    % Computes [da; de; dM] using Gauss' equations with optional drag
    a = Y(1); e = Y(2); M = Y(3);
    persistent mu A m rho0 H R
    
    mu = 4.282837e4; A = 37.5e-6; m = 2180;
    rho0 = 0.020*1e9; H = 10.8/1e3; R = 3396.2;
    
    E = kepler_solver(M, e);
    nu = 2*atan2(sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2));
    
    p = a*(1-e^2);
    r = p/(1+e*cos(nu));
    h = sqrt(mu*p);
    
    Vr = (h/p)*e*sin(nu);
    Vt = (h/p)*(1+e*cos(nu));
    v  = sqrt(Vr^2 + Vt^2);
    
    if drag_on && (r < r_thresh)
        rho = rho0*exp(-(r-R)/H);
        Cd  = max(1.502 - 0.032*log(rho/1e9), 0.1);
        aD  = 0.5*rho*Cd*A/m*v^2;
        aR  = -aD*(Vr/v);
        aT  = -aD*(Vt/v);
    else
        aR = 0; aT = 0;
    end
    
    da_dt = (2*a^2/h)*(e*sin(nu)*aR + (p/r)*aT);
    de_dt = (sqrt(p/mu)/a)*(aR*sin(nu) + aT*(cos(nu) + (e+cos(nu))/(1+e*cos(nu))));
    dM_dt = sqrt(mu/a^3) - (2/a)*da_dt - (sqrt(1-e^2)/e)*de_dt*sin(nu);
    
    dY = [da_dt; de_dt; dM_dt];
end

function [tvals, Yvals, Y_event] = integrate_until_event(Y0, dt, T_max, drag_on, r_thresh)
    t = 0; Y = Y0(:); tvals = t; Yvals = Y.';
    r_prev = compute_r(Y);
    while t < T_max
        t_next = t + dt;
        Y_next = ode4(@(t,y) orb_elem_deriv(t,y,drag_on,r_thresh), [t t_next], Y);
        Y_next = Y_next(end,:).'; r_next = compute_r(Y_next);
        if (drag_on&&(r_prev<r_thresh&&r_next>=r_thresh)) || ...
           (~drag_on&&(r_prev>r_thresh&&r_next<=r_thresh))
            Y_event = interp_event(Y, Y_next, r_prev, r_next, r_thresh);
            return;
        end
        t = t_next; Y = Y_next; r_prev = r_next;
        tvals(end+1)=t; Yvals(end+1,:)=Y.';
    end
    Y_event = Y;
end

function r = compute_r(Y)
    E = kepler_solver(Y(3), Y(2));
    nu = 2*atan2(sqrt(1+Y(2))*sin(E/2), sqrt(1-Y(2))*cos(E/2));
    p = Y(1)*(1-Y(2)^2);
    r = p/(1+Y(2)*cos(nu));
end

function Y_event = interp_event(Y1, Y2, r1, r2, r_thresh)
    alpha = (r_thresh - r1)/(r2 - r1);
    Y_event = Y1 + alpha*(Y2-Y1);
end

function E = kepler_solver(M, e)
    M = mod(M,2*pi); E = M;
    for i=1:50
        dE = -(E - e*sin(E) - M)/(1 - e*cos(E));
        E = E + dE;
        if abs(dE)<1e-8, break; end
    end
end

function Y = ode4(odefun, tspan, y0)
    h = diff(tspan); Y = y0(:);
    for k=2:length(tspan)
        ti = tspan(k-1); yi = Y(:,end);
        k1 = odefun(ti,yi); k2 = odefun(ti+h/2, yi+h*k1/2);
        k3 = odefun(ti+h/2, yi+h*k2/2); k4 = odefun(tspan(k), yi+h*k3);
        Y = [Y, yi + h/6*(k1+2*k2+2*k3+k4)];
    end
    Y = Y.';
end

function out = ternary(cond, a, b)
    if cond, out = a; else out = b; end
end
