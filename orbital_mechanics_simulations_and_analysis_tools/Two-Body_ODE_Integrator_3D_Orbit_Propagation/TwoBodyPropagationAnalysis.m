% =========================================================================
% Two-Body Orbital Propagation and Tolerance Study
%
% Description:
%   Numerically integrates the two-body equations of motion for a spacecraft 
%   around Earth under default and custom integration tolerances. 
%   Compares trajectories and saves position/velocity data for each case.
%
% Author: Erin McNeil
% Last updated: 2025-06-09
% =========================================================================

clear; close all; clc;
format longG;

%% --- Earth and Initial Conditions ---
mu = 3.986e5;              % Earth's gravitational parameter [km^3/s^2]
Re = 6378.1;               % Earth radius [km]

r0 = [3858.213; -5798.143;  14.693];  % Initial position [km]
v0 = [-0.863;    -0.542;   7.497];    % Initial velocity [km/s]
y0 = [r0; v0];                         % Initial state [r; v]

t_final = 30000;           % Propagation duration [s]

%% --- Integration Cases ---
cases = { ...
    struct('name','default','dt',0.1, 'opts',[]), ...
    struct('name','tol1',   'dt',50,  'opts',odeset('RelTol',1e-5,'AbsTol',1e-7)), ...
    struct('name','tol2',   'dt',50,  'opts',odeset('RelTol',1e-10,'AbsTol',1e-12)) ...
};

results = struct();

for idx = 1:numel(cases)
    cfg = cases{idx};
    tspan = 0:cfg.dt:t_final;
    if isempty(cfg.opts)
        sol = ode45(@two_body_eom, tspan, y0);
    else
        sol = ode45(@two_body_eom, tspan, y0, cfg.opts);
    end

    % Extract time and state
    t = sol.x;                % 1×N
    y = sol.y;                % 6×N
    r = y(1:3,:);             % 3×N position [km]
    v = y(4:6,:);             % 3×N velocity [km/s]

    % Store
    results.(cfg.name).t = t;
    results.(cfg.name).r = r;
    results.(cfg.name).v = v;

    % Display final state
    fprintf('Case: %s | Final r = [%.3f, %.3f, %.3f] km | Final v = [%.3f, %.3f, %.3f] km/s\n', ...
        cfg.name, r(:,end), v(:,end));
end

%% --- Save Data ---
save('TwoBodyPropagationData.mat','results','-v7.3');

%% --- Plot Results ---
figure; hold on; axis equal; grid on;
% Plot Earth sphere
[xe, ye, ze] = ellipsoid(0,0,0,Re,Re,Re,20);
surf(xe, ye, ze, 'FaceColor',[0.5 0.7 1],'EdgeColor','none','FaceAlpha',0.3);

% Plot trajectories
colors = lines(numel(cases));
for idx = 1:numel(cases)
    cfg = cases{idx};
    r = results.(cfg.name).r;
    plot3(r(1,:), r(2,:), r(3,:), 'LineWidth',1.2, 'Color', colors(idx,:));
end

legend(['Earth', cases{1}.name, cases{2}.name, cases{3}.name], 'Location','best');
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-Body Orbital Trajectories under Different Tolerances');
saveas(gcf,'TwoBodyTrajectories.png');
