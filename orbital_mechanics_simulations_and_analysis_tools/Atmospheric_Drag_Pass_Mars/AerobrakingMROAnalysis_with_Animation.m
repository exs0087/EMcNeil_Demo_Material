%% AE 6353 - HW4 - Problem 2 - Aerobraking of MRO via Drag Passes
% This script simulates 3 drag passes of the Mars Reconnaissance Orbiter
% and animates the spacecraft trajectory in 3D.

clear; clc;

%% Declare global variables
global mu A m r_thresh rho0 H R_mars

%% Define constants (using km and seconds)
mu = 42828;           % Mars gravitational parameter, km^3/s^2
R_mars = 3396;        % Mars radius in km
A = 37.5e-6;          % Reference area in km^2
m = 2180;             % mass in kg

% Atmospheric model parameters:
rho0 = 0.020;         % kg/m^3 at surface 
rho0 = rho0 * 1e9;    % convert to kg/km^3
H = 10800/1000;       % scale height in km

% Define density threshold (in kg/km^3):
rho_thresh = 1e-10 * 1e9;

% Solve for altitude (h_thresh) and corresponding radius (r_thresh):
h_thresh = -H * log(rho_thresh / rho0);
r_thresh = R_mars + h_thresh;   % km threshold radius.
fprintf('Atmospheric density threshold corresponds to r = %.2f km\n', r_thresh);

%% Initial orbit parameters
a0 = 25685.77;      % km
e0 = 0.8619;
rp = 3546.19;       % km
p0 = a0*(1-e0^2);
cos_nu_entry = (p0/r_thresh - 1)/e0;
cos_nu_entry = min(max(cos_nu_entry, -1), 1);
nu_entry = -acos(cos_nu_entry);
E_entry = 2*atan2(sqrt(1-e0)*sin(nu_entry/2), sqrt(1+e0)*cos(nu_entry/2));
M_entry = E_entry - e0*sin(E_entry);
Y0 = [a0; e0; M_entry];

%% Integration parameters
dt = 1;         % time step in seconds
T_max = 1e5;    % maximum integration time [s]

%% Perform 3 drag passes
pass_count = 0;
Y_entry = Y0;
a_history = [];
e_history = [];
t_full = [];
r_full = [];
x_full = [];
y_full = [];
z_full = [];

while pass_count < 3
    [t_drag, Y_drag, Y_event] = integrate_until_event(Y_entry, dt, T_max, true, 'exit');
    Y_exit = Y_event;
    [a_new, e_new, ~, ~, ~, nu_new] = orbitalElements_from_elem(Y_exit, mu);
    pass_count = pass_count + 1;
    fprintf('\nAfter drag pass %d:\n', pass_count);
    fprintf('  Semi-major axis: %.6f km\n', a_new);
    fprintf('  Eccentricity: %.6f\n', e_new);
    a_history(end+1) = a_new;
    e_history(end+1) = e_new;

    [r_drag, nu_drag] = compute_r_from_Yarray(Y_drag);
    t_full = [t_full; t_drag + (pass_count-1)*T_max*2];
    r_full = [r_full; r_drag];
    x_drag = r_drag .* cos(nu_drag);
    y_drag = r_drag .* sin(nu_drag);
    z_drag = zeros(size(x_drag));
    x_full = [x_full; x_drag; NaN];
    y_full = [y_full; y_drag; NaN];
    z_full = [z_full; z_drag; NaN];

        % Propagate full vacuum orbit using Keplerian motion (1 full period)
    a = Y_exit(1);
    T_orbit = 2 * pi * sqrt(a^3 / mu);
    N_steps = ceil(T_orbit / dt);
    t_vac = linspace(0, T_orbit, N_steps)';
    Y_vac = zeros(N_steps, 3);
    Y_vac(1,:) = Y_exit';
    for j = 2:N_steps
        Y_temp = ode4(@(t,Y) orb_elem_deriv(t,Y,false), [0 dt], Y_vac(j-1,:)');
        Y_vac(j,:) = Y_temp(end,:);
    end
    Y_entry = Y_vac(end,:)';
    [~, ~, ~, ~, ~, nu_entry_new] = orbitalElements_from_elem(Y_entry, mu);
    fprintf('  (New entry state) True anomaly: %.6f rad\n', nu_entry_new);

    [r_vac, nu_vac] = compute_r_from_Yarray(Y_vac);
    t_full = [t_full; t_vac + (pass_count-1)*T_max*2 + t_drag(end)];
    r_full = [r_full; r_vac];
    x_vac = r_vac .* cos(nu_vac);
    y_vac = r_vac .* sin(nu_vac);
    z_vac = zeros(size(x_vac));
    x_full = [x_full; x_vac; NaN];
    y_full = [y_full; y_vac; NaN];
    z_full = [z_full; z_vac; NaN];
end

fprintf('\nAfter 3 drag passes:\n');
fprintf('  Initial semi-major axis: %.6f km\n', a0);
fprintf('  Final semi-major axis: %.6f km\n', a_history(end));
if abs(e_history(end) - 1) < 1e-3
    orbit_type = 'Parabolic';
elseif e_history(end) < 1
    orbit_type = 'Elliptical';
else
    orbit_type = 'Hyperbolic';
end
fprintf('Final orbit type: %s\n', orbit_type);

%% Animate the spacecraft trajectory in 3D
figure;
% Plot Mars as a translucent sphere
[XM, YM, ZM] = sphere(50);
surf(XM * R_mars, YM * R_mars, ZM * R_mars, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', [1 0.3 0.3]);
hold on;

% Break trajectory into individual segments (drag + vacuum arcs)
colors = lines(6); % Predefine up to 6 unique colors
segment_idx = find(isnan(x_full));
segment_starts = [1; segment_idx(1:end-1)+1];
segment_ends = segment_idx - 1;

% Plot each orbit segment in a unique color
orbit_labels = {"Drag Pass 1 Orbit", "Drag Pass 2 Orbit", "Drag Pass 3 Orbit"};
segment_colors = lines(length(segment_starts));

orbit_lines = gobjects(length(segment_starts), 1);
for i = 1:length(segment_starts)
    idx = segment_starts(i):segment_ends(i);
    label = "";
    if mod(i,2) == 0 && (i/2 <= length(orbit_labels))
        label = orbit_labels{i/2};
    end
    orbit_lines(i) = plot3(x_full(idx), y_full(idx), z_full(idx), '-', 'Color', segment_colors(i,:), 'LineWidth', 1.5, 'DisplayName', label);
end

% legend(orbit_lines(orbit_lines ~= 0), 'Location', 'best');
axis equal;
axis vis3d;
grid on;
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
title('3D Animation of MRO Aerobraking Trajectory');
view([120 30]);

trail_length = 300;
drag_pass_text = text(0, 0, 2*R_mars, 'Drag Pass #1', 'FontSize', 12, 'FontWeight', 'bold');
h_marker = plot3(x_full(1), y_full(1), z_full(1), 'ko', 'MarkerFaceColor', 'k');

for k = 1:100:length(x_full)
    idx = max(1, k-trail_length):k;
    if any(isnan(x_full(k)))
        continue;
    end
    set(h_marker, 'XData', x_full(k), 'YData', y_full(k), 'ZData', z_full(k));
    trail = plot3(x_full(idx), y_full(idx), z_full(idx), 'k-', 'LineWidth', 1);

    % Update orbit number label
    orbit_num = sum(isnan(x_full(1:k))) + 1;
    drag_pass_num = ceil(orbit_num / 2);
    set(drag_pass_text, 'String', sprintf('Drag Pass #%d', drag_pass_num));
    drawnow;
end

%% Function: Orbital Element Rates (Gauss equations)
function dYdt = orb_elem_deriv(~, Y, drag_on)
    global mu A m r_thresh rho0 H R_mars
    a = Y(1); e = Y(2); M = Y(3);
    E = kepler_solver(M, e);
    nu = 2 * atan2(sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2));
    p = a*(1-e^2);      
    r = p/(1+e*cos(nu));
    h = sqrt(mu*p);
    Vr = (h/p)*e*sin(nu);
    Vt = (h/p)*(1+e*cos(nu));
    v = sqrt(Vr^2+Vt^2);

    if drag_on && (r < r_thresh)
        rho = max(rho0 * exp(-(r-R_mars)/H), 1e-30);
        rho_m = rho / 1e9;
        Cd = max(1.502 - 0.032*log(rho_m), 0.1);
        a_D = 0.5 * rho * Cd * A/m * v^2;
        a_R = -a_D * (Vr/v);
        a_T = -a_D * (Vt/v);
    else
        a_R = 0;
        a_T = 0;
    end

    da_dt = (2*a^2/h) * ( e*sin(nu)*a_R + (p/r)*a_T );
    de_dt = (sqrt(p/mu)/a) * ( a_R*sin(nu) + a_T*(cos(nu) + (e+cos(nu))/(1+e*cos(nu)) ) );
    n = sqrt(mu/a^3);
    dM_dt = n - (2/(a))*da_dt - (sqrt(1-e^2)/e)*de_dt*sin(nu);

    dYdt = [da_dt; de_dt; dM_dt];
end

%% Function: Integrate Until an Event Occurs
function [t_vals, Y_vals, Y_event] = integrate_until_event(Y0, dt, T_max, drag_on, event_type)
    global r_thresh R_mars
    t = 0;
    Y = Y0;
    t_vals = t;
    Y_vals = Y';
    r_prev = compute_r(Y);
    event_occurred = false;

    while t < T_max && ~event_occurred
        t_next = t + dt;
        Y_next = ode4(@(t,Y) orb_elem_deriv(t,Y,drag_on), [t, t_next], Y);
        Y_next = Y_next(end,:)';
        r_next = compute_r(Y_next);

        if strcmpi(event_type, 'exit') && (r_prev < r_thresh) && (r_next >= r_thresh)
            event_occurred = true;
        elseif strcmpi(event_type, 'entry') && (r_prev > r_thresh) && (r_next <= r_thresh)
            event_occurred = true;
        end
        t = t_next;
        Y = Y_next;
        t_vals(end+1,1) = t;
        Y_vals(end+1,:) = Y';
        r_prev = r_next;
    end

    if size(Y_vals,1) < 2
        Y_event = Y_vals(end,:)';
        return;
    end
    r1 = compute_r(Y_vals(end-1,:)');
    r2 = compute_r(Y_vals(end,:)');
    if abs(r2 - r1) < 1e-12
        alpha = 0;
    else
        alpha = (r_thresh - r1) / (r2 - r1);
    end
    t1 = t_vals(end-1); t2 = t_vals(end);
    Y1 = Y_vals(end-1,:)'; Y2 = Y_vals(end,:)';
    Y_event = Y1 + alpha*(Y2 - Y1);
end

%% Function: Compute Radius from Orbital Elements
function r = compute_r(Y)
    global R_mars
    a = Y(1); e = Y(2); M = Y(3);
    E = kepler_solver(M, e);
    nu = 2*atan2(sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2));
    p = a*(1-e^2);
    r = p/(1+e*cos(nu));
end

%% Function: Radius and True Anomaly from Y array
function [r_vals, nu_vals] = compute_r_from_Yarray(Y_array)
    r_vals = zeros(size(Y_array,1),1);
    nu_vals = zeros(size(Y_array,1),1);
    for k = 1:size(Y_array,1)
        a = Y_array(k,1);
        e = Y_array(k,2);
        M = Y_array(k,3);
        E = kepler_solver(M, e);
        nu = 2*atan2(sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2));
        p = a*(1-e^2);
        r = p/(1+e*cos(nu));
        r_vals(k) = r;
        nu_vals(k) = nu;
    end
end

%% Function: Orbital Elements from State
function [a, e, i, RAAN, argp, nu] = orbitalElements_from_elem(Y, mu)
    a = Y(1);
    e = Y(2);
    E = kepler_solver(Y(3), e);
    nu = mod(2*atan2(sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2)), 2*pi);
    i = 0; RAAN = 0; argp = 0;
end

%% Function: Solve Kepler's Equation
function E = kepler_solver(M, e)
    M = mod(M, 2*pi);
    E = M;
    for iter = 1:100
        f = E - e*sin(E) - M;
        fp = 1 - e*cos(E);
        dE = -f/fp;
        E = E + dE;
        if abs(dE) < 1e-8
            break;
        end
    end
end

%% Fixed-Step 4th-Order Runge-Kutta Integrator
function Y = ode4(odefun, tspan, y0, varargin)
    h = diff(tspan);
    Y = y0(:);
    for i = 2:length(tspan)
        ti = tspan(i-1);
        yi = Y(:, end);
        k1 = odefun(ti, yi, varargin{:});
        k2 = odefun(ti+0.5*h, yi+0.5*h*k1, varargin{:});
        k3 = odefun(ti+0.5*h, yi+0.5*h*k2, varargin{:});
        k4 = odefun(tspan(i), yi+h*k3, varargin{:});
        Y_next = yi + h/6*(k1+2*k2+2*k3+k4);
        Y = [Y, Y_next];
    end
    Y = Y.';
end
