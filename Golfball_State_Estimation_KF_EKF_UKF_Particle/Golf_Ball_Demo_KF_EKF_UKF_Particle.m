%% Golf_Ball_Demo_KF_EKF_UKF_Particle.m
% Full demo comparing state-estimation filters on golf-ball trajectory data
% Methods: LUMVE, Standard KF, Extended KF, Unscented KF, Particle Filter
% Requires golfball_flight_data.xlsx with columns: time, range, azimuth, elevation

clear all; close all; clc;
randn('state', sum(clock));  % Seed RNG

%%--------------------------
% Read measurement data
%%--------------------------
Data      = readmatrix('golfball_flight_data.xlsx');
time      = Data(:,1);      % time [s]
range     = Data(:,2);      % range [m]
azimuth   = Data(:,3);      % azimuth [rad]
elevation = Data(:,4);      % elevation [rad]

%%--------------------------
% Physical & noise params
%%--------------------------
m   = 45.93e-3;          % mass [kg]
w   = [250; 0; 0];       % spin vector [rad/s]
g   = 9.81;              % gravity [m/s^2]
rho = 1.293;             % air density [kg/m^3]
dia = 42.67e-3;          % ball diameter [m]
A   = pi*(dia/2)^2;      % cross-sectional area [m^2]
Cd  = 0.2;               % drag coefficient
D   = 0.5*Cd*A*rho;      % drag constant
S   = 2.5e-5;            % Magnus constant

% Measurement noise std deviations
range_MeasNoise = 0.9144;    % [m]
az_MeasNoise    = 0.00174533; % [rad]
el_MeasNoise    = 0.00174533; % [rad]

% Covariances
R   = diag([range_MeasNoise^2, az_MeasNoise^2, el_MeasNoise^2]);
Rin = inv(R);
P0  = diag([1,1,1,0.1,0.1,0.1]); % initial covariance

%%====================
% 1) LUMVE Filter
%%====================
% initial states
xc  = [0.01; 0.02; 0.03; 0; 70.43; 15];   % initial guess
xc0 = xc;   % a priori
xcc = xc;   % corrected

% simulation parameters
dt = 0.1; tf = 7; dtPlot = dt; tPlot = -inf;
xhat_LUMVE = []; t_LUMVE = [];
stp = 0; iter = 0;

for t = 0:dt:tf
    iter = iter + 1;
    % save for plotting
    if t >= tPlot + dtPlot
        tPlot      = t + dtPlot - eps;
        xhat_LUMVE = [xhat_LUMVE, xcc];
        t_LUMVE    = [t_LUMVE, t];
    end
    % propagate
    xcdot = [
        xc(4);
        xc(5);
        xc(6);
        -D/m*xc(4)^2 + S/m*w(2)*xc(6) - S/m*w(3)*xc(5);
        -D/m*xc(5)^2 + S/m*w(3)*xc(4) - S/m*w(1)*xc(6);
        -D/m*xc(6)^2 + S/m*w(1)*xc(5) - S/m*w(2)*xc(4) - g
    ];
    xc = xc + xcdot * dt;
    % Jacobian H
    r    = norm(xc(1:3)); rho2 = norm(xc(1:2));
    H = [
        xc(1)/r, xc(2)/r, xc(3)/r, 0,0,0;
       -xc(2)/rho2^2, xc(1)/rho2^2, 0,0,0,0;
        xc(1)*xc(3)/(r^2*rho2), xc(2)*xc(3)/(r^2*rho2), -rho2/r^2, 0,0,0
    ];
    % measurement update timing
    if iter == dtPlot/dt && stp < numel(range)
        stp    = stp + 1;
        y_meas = [range(stp); azimuth(stp); elevation(stp)];
        iter   = 0;
    elseif t == 0
        stp    = 1;
        y_meas = [range(1); azimuth(1); elevation(1)];
    end
    % iterative correction
    xcc = xc; P = P0;
    for k = 1:10
        pr     = norm(xcc(1:3));
        pb     = atan2(xcc(2), xcc(1));
        pe     = atan2(xcc(3), norm(xcc(1:2)));
        y_pred = [pr; pb; pe] + [range_MeasNoise; az_MeasNoise; el_MeasNoise] .* randn(3,1);
        dy = y_meas - y_pred;
        M1 = H'*Rin*H + inv(P0);
        dx = M1 \ (H'*Rin*dy + inv(P0)*xc0);
        P  = inv(M1);
        xcc = xcc - dx;
    end
end

%%====================
% Truth via ODE45
%%====================
tspan = [0 5];
[tt_truth, X_truth] = ode45(@(t,X) odefcn_golfball_dynamics(t,X,D,S), tspan, xc0);
% align times
X_truth_L = zeros(numel(t_LUMVE), 6);
for i = 1:numel(t_LUMVE)
    [~, idx] = min(abs(tt_truth - t_LUMVE(i)));
    X_truth_L(i,:) = X_truth(idx,:);
end

%% Plot LUMVE
figure;
x1m = range.*cos(elevation).*cos(azimuth);
x2m = range.*cos(elevation).*sin(azimuth);
x3m = range.*sin(elevation);
plot3(x1m, x2m, x3m, 'mO'); hold on;
plot3(X_truth_L(:,1), X_truth_L(:,2), X_truth_L(:,3), 'm-');
plot3(xhat_LUMVE(1,:), xhat_LUMVE(2,:), xhat_LUMVE(3,:), 'b-');
grid on; view(45,45);
title('LUMVE Estimate'); xlabel('X'); ylabel('Y'); zlabel('Z');
legend('Meas','Truth','LUMVE');

% Reset variables
clearvars -except range azimuth elevation D S w m g rho dia A Cd range_MeasNoise az_MeasNoise el_MeasNoise
% recompute Cartesian measurement vectors
x1m = range.*cos(elevation).*cos(azimuth);
x2m = range.*cos(elevation).*sin(azimuth);
x3m = range.*sin(elevation);
% reinitialize initial state for consistency
xc0 = [0.01; 0.02; 0.03; 0; 70.43; 15];
% reinitialize initial state
xc0 = [0.01; 0.02; 0.03; 0; 70.43; 15];
Data  = readmatrix('golfball_flight_data.xlsx');
range = Data(:,2); azimuth = Data(:,3); elevation = Data(:,4);
R     = diag([range_MeasNoise^2, az_MeasNoise^2, el_MeasNoise^2]);
Q     = zeros(6);
P     = diag([1,1,1,0.1,0.1,0.1]);
xhat  = xc0;
% constant A, H0
a44 = -2*D/m * xhat(4); a45 = -S/m * w(3); a46 =  S/m * w(2);
a54 =  S/m * w(3); a55 = -2*D/m * xhat(5); a56 = -S/m * w(1);
a64 = -S/m * w(2); a65 =  S/m * w(1); a66 = -2*D/m * xhat(6);
A = [0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1; 0 0 0 a44 a45 a46;
     0 0 0 a54 a55 a56; 0 0 0 a64 a65 a66];
r0 = norm(xhat(1:3)); rho0 = norm(xhat(1:2));
H0 = [xhat(1)/r0, xhat(2)/r0, xhat(3)/r0, 0,0,0;
      -xhat(2)/rho0^2, xhat(1)/rho0^2,0,0,0,0;
       xhat(1)*xhat(3)/(r0^2*rho0), xhat(2)*xhat(3)/(r0^2*rho0), -rho0/r0^2,0,0,0];

dt=0.1; tf=7; dtPlot=dt; tPlot=-inf;
xhat_KF = []; trP_KF = []; t_KF = [];
stp=0; iter=0;
for t=0:dt:tf
    iter = iter + 1;
    if t >= tPlot+dtPlot
        tPlot   = t+dtPlot-eps;
        xhat_KF = [xhat_KF, xhat];
        trP_KF  = [trP_KF, trace(P)];
        t_KF    = [t_KF, t];
    end
    % Predict
    P    = A*P*A' + Q;
    cv   = cross(w, xhat(4:6));
    xdot = [xhat(4); xhat(5); xhat(6);
            -D/m*xhat(4)^2 + S/m*cv(1);
            -D/m*xhat(5)^2 + S/m*cv(2);
            -D/m*xhat(6)^2 + S/m*cv(3) - g];
    xhat = xhat + xdot*dt;
    % Update
    if iter==dtPlot/dt && stp<numel(range)
        stp     = stp+1;
        y_meas  = [range(stp); azimuth(stp); elevation(stp)];
        iter    = 0;
    elseif t==0
        stp    = 1;
        y_meas = [range(1); azimuth(1); elevation(1)];
    end
    K    = P*H0'/(H0*P*H0' + R);
    xhat = xhat + K*(y_meas - H0*xhat);
    P    = (eye(6) - K*H0)*P;
end

% compute Cartesian measurements for plotting
x1m = range.*cos(elevation).*cos(azimuth);
x2m = range.*cos(elevation).*sin(azimuth);
x3m = range.*sin(elevation);
% compute Cartesian measurements for plotting
x1m = range.*cos(elevation).*cos(azimuth);
x2m = range.*cos(elevation).*sin(azimuth);
x3m = range.*sin(elevation);
% compute Cartesian measurements for plotting
x1m = range.*cos(elevation).*cos(azimuth);
x2m = range.*cos(elevation).*sin(azimuth);
x3m = range.*sin(elevation);
% compute Cartesian measurements for plotting
x1m = range.*cos(elevation).*cos(azimuth);
x2m = range.*cos(elevation).*sin(azimuth);
x3m = range.*sin(elevation);
figure;
plot3(x1m,x2m,x3m,'rO'); hold on;
plot3(xhat_KF(1,:),xhat_KF(2,:),xhat_KF(3,:),'g-');
grid on; view(45,45);
title('Standard KF'); legend('Meas','KF');

%%============================
% 3) Extended Kalman Filter (EKF)
%%============================
clearvars -except range azimuth elevation D S w m g rho dia A Cd range_MeasNoise az_MeasNoise el_MeasNoise
% recompute Cartesian measurement vectors
x1m = range.*cos(elevation).*cos(azimuth);
x2m = range.*cos(elevation).*sin(azimuth);
x3m = range.*sin(elevation);
% reinitialize initial state for consistency
xc0 = [0.01; 0.02; 0.03; 0; 70.43; 15];
Data  = readmatrix('golfball_flight_data.xlsx');
range = Data(:,2); azimuth = Data(:,3); elevation = Data(:,4);
R = diag([range_MeasNoise^2, az_MeasNoise^2, el_MeasNoise^2]);
Q = zeros(6);
P = diag([1,1,1,0.1,0.1,0.1]);
xhat = xc0;

dt=0.1; tf=7; dtPlot=dt; tPlot=-inf;
xhat_EKF = []; trP_EKF = []; t_EKF = [];
stp=0; iter=0;
for t=0:dt:tf
    iter=iter+1;
    if t>=tPlot+dtPlot
        tPlot     = t+dtPlot-eps;
        xhat_EKF  = [xhat_EKF, xhat];
        trP_EKF   = [trP_EKF, trace(P)];
        t_EKF     = [t_EKF, t];
    end
    % Linearize A
    a44=-2*D/m*xhat(4); a45=-S/m*w(3); a46= S/m*w(2);
    a54= S/m*w(3);     a55=-2*D/m*xhat(5); a56=-S/m*w(1);
    a64=-S/m*w(2);     a65= S/m*w(1);     a66=-2*D/m*xhat(6);
    A=[0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1;0 0 0 a44 a45 a46;0 0 0 a54 a55 a56;0 0 0 a64 a65 a66];
    P=A*P*A'+Q;
    cv=cross(w,xhat(4:6));
    xdot=[xhat(4);xhat(5);xhat(6);-D/m*xhat(4)^2+S/m*cv(1);-D/m*xhat(5)^2+S/m*cv(2);-D/m*xhat(6)^2+S/m*cv(3)-g];
    xhat = xhat + xdot*dt;
    % Jacobian H
    r=norm(xhat(1:3)); rho2=norm(xhat(1:2));
    H=[xhat(1)/r,xhat(2)/r,xhat(3)/r,0,0,0; -xhat(2)/rho2^2,xhat(1)/rho2^2,0,0,0,0; xhat(1)*xhat(3)/(r^2*rho2),xhat(2)*xhat(3)/(r^2*rho2),-rho2/r^2,0,0,0];
    % get measurement
    if iter==dtPlot/dt && stp<numel(range)
        stp    = stp+1;
        y_meas = [range(stp); azimuth(stp); elevation(stp)];
        iter   = 0;
    elseif t==0
        stp    = 1;
        y_meas = [range(1); azimuth(1); elevation(1)];
    end
    K    = P*H'/(H*P*H'+R);
    pr   = norm(xhat(1:3)); pb=atan2(xhat(2),xhat(1)); pe=atan2(xhat(3),norm(xhat(1:2)));
    yhat = [pr;pb;pe];
    xhat = xhat + K*(y_meas - yhat);
    P    = (eye(6)-K*H)*P;
end
figure;
plot3(x1m,x2m,x3m,'cO'); hold on;
plot3(xhat_EKF(1,:),xhat_EKF(2,:),xhat_EKF(3,:),'c-');
grid on; view(45,45);
title('EKF'); legend('Meas','EKF');

%%============================
% 4) Unscented KF (UKF)
%%============================
clearvars -except range azimuth elevation D S w m g rho dia A Cd range_MeasNoise az_MeasNoise el_MeasNoise
% recompute Cartesian measurement vectors
x1m = range.*cos(elevation).*cos(azimuth);
x2m = range.*cos(elevation).*sin(azimuth);
x3m = range.*sin(elevation);
% reinitialize initial state for consistency
xc0 = [0.01; 0.02; 0.03; 0; 70.43; 15];
Data = readmatrix('golfball_flight_data.xlsx');
range = Data(:,2); azimuth=Data(:,3); elevation=Data(:,4);
R = diag([range_MeasNoise^2, az_MeasNoise^2, el_MeasNoise^2]);
Q = zeros(6);
P = diag([1,1,1,0.1,0.1,0.1]);
xhatUKF = xc0; Pukf = P;

% UKF weights and params
dt=0.1; tf=7; dtt=0.005;
n=6; T=dt;
W = 1/(2*n) * ones(2*n,1);
stp=0; iter=0;
xhat_UKF=[]; t_UKF=[];

for t=0:dt:tf
    iter=iter+1;
    if isempty(t_UKF) || t>=t_UKF(end)+dt
        xhat_UKF=[xhat_UKF, xhatUKF]; t_UKF=[t_UKF, t];
    end
    % generate sigma points
    [L,~]=chol(n*Pukf,'lower');
    sigmas = [xhatUKF+L, xhatUKF-L];
    % time update
    for i=1:2*n
        xi=sigmas(:,i);
        for tau=dtt:dtt:T
            cv=cross(w,xi(4:6));
            xdot=[xi(4);xi(5);xi(6); -D/m*xi(4)^2+S/m*cv(1); -D/m*xi(5)^2+S/m*cv(2); -D/m*xi(6)^2+S/m*cv(3)-g];
            xi=xi+xdot*dtt;
        end
        sigmas(:,i)=xi;
    end
    % predict mean+cov
    xhatUKF = sigmas*W;
    Pukf = Q;
    for i=1:2*n
        d=sigmas(:,i)-xhatUKF; Pukf=Pukf+W(i)*(d*d');
    end
    % measurement update
    % new measurement
    if iter==dt/dt && stp<numel(range)
        stp=stp+1; y_meas=[range(stp); azimuth(stp); elevation(stp)]; iter=0;
    elseif t==0
        stp=1; y_meas=[range(1); azimuth(1); elevation(1)];
    end
    [L,~]=chol(n*Pukf,'lower');
    sigmas=[xhatUKF+L, xhatUKF-L];
    zsig=zeros(3,2*n);
    for i=1:2*n
        xi=sigmas(:,i);
        pr=norm(xi(1:3)); pb=atan2(xi(2),xi(1)); pe=atan2(xi(3),norm(xi(1:2)));
        zsig(:,i)=[pr;pb;pe];
    end
    zhat=zsig*W;
    Py=R; Pxy=zeros(n,3);
    for i=1:2*n
        dz=zsig(:,i)-zhat; dx=sigmas(:,i)-xhatUKF;
        Py=Py+W(i)*(dz*dz'); Pxy=Pxy+W(i)*(dx*dz');
    end
    K=Pxy/Py;
    innov=y_meas-zhat;
    xhatUKF=xhatUKF+K*innov;
    Pukf=Pukf-K*Py*K';
end

figure;
plot3(x1m,x2m,x3m,'kO'); hold on;
plot3(xhat_UKF(1,:), xhat_UKF(2,:), xhat_UKF(3,:), 'r-');
grid on; view(45,45);
title('UKF'); legend('Meas','UKF');

%%============================
% 5) Particle Filter (PF)
%%============================
clearvars -except range azimuth elevation D S w m g rho dia A Cd range_MeasNoise az_MeasNoise el_MeasNoise
% recompute Cartesian measurement vectors
x1m = range.*cos(elevation).*cos(azimuth);
x2m = range.*cos(elevation).*sin(azimuth);
x3m = range.*sin(elevation);
% reinitialize initial state for consistency
xc0 = [0.01; 0.02; 0.03; 0; 70.43; 15];
Data = readmatrix('golfball_flight_data.xlsx');
range = Data(:,2); azimuth=Data(:,3); elevation=Data(:,4);
R = diag([range_MeasNoise^2, az_MeasNoise^2, el_MeasNoise^2]);

% initialize particles
Npart=1000;
P0=diag([1,1,1,0.1,0.1,0.1]);
particles = repmat(xc0,1,Npart) + chol(P0,'lower')*randn(6,Npart);
xhat_PF=xc0;

% sim parameters
dt=0.1; tf=7;
xhat_PF_hist=[]; t_PF=[];
stp=0; iter=0;
for t=0:dt:tf
    iter=iter+1;
    if isempty(t_PF) || t>=t_PF(end)+dt
        xhat_PF_hist=[xhat_PF_hist, xhat_PF]; t_PF=[t_PF, t];
    end
    % propagate particles
    for i=1:Npart
        xi=particles(:,i);
        cv=cross(w,xi(4:6));
        xdot=[xi(4);xi(5);xi(6); -D/m*xi(4)^2+S/m*cv(1); -D/m*xi(5)^2+S/m*cv(2); -D/m*xi(6)^2+S/m*cv(3)-g];
        particles(:,i)=xi+xdot*dt;
    end
    % measurement
    if iter==dt/dt && stp<numel(range)
        stp=stp+1; y_meas=[range(stp); azimuth(stp); elevation(stp)]; iter=0;
    elseif t==0
        stp=1; y_meas=[range(1); azimuth(1); elevation(1)];
    end
    % weight & resample
    wts=zeros(Npart,1);
    for i=1:Npart
        xi=particles(:,i);
        pr=norm(xi(1:3)); pb=atan2(xi(2),xi(1)); pe=atan2(xi(3),norm(xi(1:2)));
        innov=y_meas-[pr;pb;pe];
        wts(i)=exp(-0.5*(innov'*(R\innov)));
    end
    wts=wts/sum(wts);
    % resample via inverse transform sampling (no toolbox required)
cs = cumsum(wts);
u = rand(Npart,1);
idx = zeros(Npart,1);
for j = 1:Npart
    % select first index where cumulative sum exceeds random uniform
    loc = find(cs >= u(j), 1);
    if isempty(loc)
        loc = Npart;
    end
    idx(j) = loc;
end
particles = particles(:, idx);
    % roughening step: spread particles to prevent impoverishment
E = max(particles,[],2) - min(particles,[],2);
    for i=1:Npart
        particles(:,i)=particles(:,i) + 0.05 * E .* randn(6,1) / Npart^(1/6);
    end
    for i=1:Npart
        particles(:,i)=particles(:,i)+0.05*E.*randn(6,1)/Npart^(1/6);
    end
    xhat_PF=mean(particles,2);
end

figure;
plot3(x1m,x2m,x3m,'bO'); hold on;
plot3(xhat_PF_hist(1,:), xhat_PF_hist(2,:), xhat_PF_hist(3,:), 'k-');
grid on; view(45,45);
title('Particle Filter'); legend('Meas','PF');
