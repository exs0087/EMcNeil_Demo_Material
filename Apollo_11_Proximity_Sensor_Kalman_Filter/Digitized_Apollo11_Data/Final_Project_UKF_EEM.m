clear all; close all;clc;
%% AE6505 Final Project -- Apollo 11 Final Descent: Unscented Kalman Filter
% Use of an UKF to improve estimates of the true altitude from unreliable
% proximity sensor data
% Author: Erin McNeil
% Date: 4-25-2024

%% Choose Noisy Sensor Test Case to Run
method = 'Test 1'; %%%%SELECT TEST CASE%%%%

%% Load Digitized Data to Workspace
load('Time_A_trunc.mat');
load('Alt_trunc.mat');

load('Time_P_trunc.mat');
load('Pitch_trunc.mat');

load('Time_RA_trunc.mat');
load('RangeAccel_trunc.mat');


Altitude = Alt_trunc; %ft
RangeAccel = RangeAccel_trunc; %ft/sec^2
Pitch = Pitch_trunc*(pi/180); %rad

%% parameter inputs
g = 5.315; %gravity (ft/sec^2)

%% noise and covariance matrices
R = 1000; % Measurement noise covariance
P = eye(2)*1; % Initial state estimation covariance
Pukfs = P;
Q = eye(2)*1; %process noise

%% initial state  
x = [400; ... %altitude (ft)
    0.5;... %altitude rate (ft/sec)
    ]; % Initial state x0
xhatukfs = x;

y_meas = 400;

%% Iterations
dt = 0.2;
tf = 150;
N = length(x);
dtt = 0.005;
T = 0.2; % measurement time step

dtPlot = dt; % How often to plot results
tPlot = -inf;
tMeas = -inf;

% Initialize arrays for plotting at the end of the program
tArray = [];
xArray = [];
y_meaArr = [];
tMArray = [];
xhatukfsArray = [];

Ws = ones(2*N,1) / 2 / N; % standard UKF weights W = 1/2n

stp = 0;
iter = 0;
for t = 0 : dt : tf
    iter = iter + 1;
    
    if t >= tPlot + dtPlot
        % Save data for plotting
        tPlot = t + dtPlot - eps;
        tArray = [tArray t];
        xArray = [xArray x];
        xhatukfsArray = [xhatukfsArray xhatukfs];
    end
    
    Range_Accel_ = RangeAccel(iter);
    Pitch_ = Pitch(iter);
    
    u = (-Range_Accel_/tan(Pitch_)) - g;
    x(1) = x(1) + x(2)*dt + (u*dt)/2;
    x(2) = x(2) + u*dt;
    x = [x(1);x(2)] + diag(sqrt(Q)*randn);
    
    y_comp = x(1) + sqrt(R)*randn; %Altitude (ft)
    
    %% Test Cases for Variable Measurement Noise
    
    switch (method)
      case 'Test 1'
      dtMeas = 0.8; %%%%UPDATE RATE%%%%
        if t >= tMeas + dtMeas 
            tMeas = t + dtMeas - eps;
                Noise = 1000; %%%%MAGNITUDE%%%%
                deg = 0; %%%%FREQUENCY%%%%
                y_meas = Altitude(iter)+ (Noise*cos(deg*t)*randn);
            y_meaArr = [y_meaArr y_meas];
            tMArray = [tMArray t];
        end
      case 'Test 2'
      dtMeas = 1.5; %%%%UPDATE RATE%%%%
        if t >= tMeas + dtMeas 
            tMeas = t + dtMeas - eps;
                Noise = 1000; %%%%MAGNITUDE%%%%
                deg = 0; %%%%FREQUENCY%%%%
                y_meas = Altitude(iter)+ (Noise*cos(deg*t)*randn);
            y_meaArr = [y_meaArr y_meas];
            tMArray = [tMArray t];
        end
      case 'Test 3'
      dtMeas = 0.8; %%%%UPDATE RATE%%%%
        if t >= tMeas + dtMeas 
            tMeas = t + dtMeas - eps;
                Noise = 1000; %%%%MAGNITUDE%%%%
                deg = 50; %%%%FREQUENCY%%%%
                y_meas = Altitude(iter)+ (Noise*cos(deg*t)*randn);
            y_meaArr = [y_meaArr y_meas];
            tMArray = [tMArray t];
        end
      case 'Test 4'
      dtMeas = 1.5; %%%%UPDATE RATE%%%%
        if t >= tMeas + dtMeas 
            tMeas = t + dtMeas - eps;
                Noise = 1000; %%%%MAGNITUDE%%%%
                deg = 50; %%%%FREQUENCY%%%%
                y_meas = Altitude(iter)+ (Noise*cos(deg*t)*randn);
            y_meaArr = [y_meaArr y_meas];
            tMArray = [tMArray t];
        end
      case 'Test 5'
      dtMeas = 0.8; %%%%UPDATE RATE%%%%
        if t >= tMeas + dtMeas 
            tMeas = t + dtMeas - eps;
                Noise = 1000; %%%%MAGNITUDE%%%%
                deg = 5; %%%%FREQUENCY%%%%
                y_meas = Altitude(iter)+ (Noise*cos(deg*t)*randn);
            y_meaArr = [y_meaArr y_meas];
            tMArray = [tMArray t];
        end
      case 'Test 6'
      dtMeas = 0.8; %%%%UPDATE RATE%%%%
        if t >= tMeas + dtMeas 
            tMeas = t + dtMeas - eps;
                Noise = 1000; %%%%MAGNITUDE%%%%
                deg = 25; %%%%FREQUENCY%%%%
                y_meas = Altitude(iter)+ (Noise*cos(deg*t)*randn);
            y_meaArr = [y_meaArr y_meas];
            tMArray = [tMArray t];
        end
      case 'Test 7'
      dtMeas = 0.8; %%%%UPDATE RATE%%%%
        if t >= tMeas + dtMeas 
            tMeas = t + dtMeas - eps;
                Noise = 1000; %%%%MAGNITUDE%%%%
                deg = 500; %%%%FREQUENCY%%%%
                y_meas = Altitude(iter)+ (Noise*cos(deg*t)*randn);
            y_meaArr = [y_meaArr y_meas];
            tMArray = [tMArray t];
        end
      case 'Test 8'
      dtMeas = 0.8; %%%%UPDATE RATE%%%%
        if t >= tMeas + dtMeas 
            tMeas = t + dtMeas - eps;
                Noise = 100; %%%%MAGNITUDE%%%%
                deg = 5; %%%%FREQUENCY%%%%
                y_meas = Altitude(iter)+ (Noise*cos(deg*t)*randn);
            y_meaArr = [y_meaArr y_meas];
            tMArray = [tMArray t];
        end
      case 'Test 9'
      dtMeas = 0.8; %%%%UPDATE RATE%%%%
        if t >= tMeas + dtMeas 
            tMeas = t + dtMeas - eps;
                Noise = 5000; %%%%MAGNITUDE%%%%
                deg = 5; %%%%FREQUENCY%%%%
                y_meas = Altitude(iter)+ (Noise*cos(deg*t)*randn);
            y_meaArr = [y_meaArr y_meas];
            tMArray = [tMArray t];
        end
      case 'Test 10'
      dtMeas = 0.8; %%%%UPDATE RATE%%%%
        if t >= tMeas + dtMeas 
            tMeas = t + dtMeas - eps;
                Noise = 7000; %%%%MAGNITUDE%%%%
                deg = 5; %%%%FREQUENCY%%%%
                y_meas = Altitude(iter)+ (Noise*cos(deg*t)*randn);
            y_meaArr = [y_meaArr y_meas];
            tMArray = [tMArray t];
        end
      case 'WORST CASE'
      dtMeas = 1.5; %%%%UPDATE RATE%%%%
        if t >= tMeas + dtMeas 
            tMeas = t + dtMeas - eps;
            Noise = 7000; %%%%MAGNITUDE%%%%
            deg = 500; %%%%FREQUENCY%%%%
            if (100 < t) && (t < 140)
                y_meas = y_meas; %%%%STALE DATA%%%%
            elseif (75 < t) && (t < 100)
                y_meas = Altitude(iter)+ (Noise*cos(deg*t)*randn);
            elseif (50 < t) && (t < 75)
                y_meas = y_meas; %%%%STALE DATA%%%%
            elseif (0 < t) && (t < 50)
                y_meas = Altitude(iter)+ (Noise*cos(deg*t)*randn);
            else
                y_meas = y_meas; %%%%STALE DATA%%%%
            end
            y_meaArr = [y_meaArr y_meas];
            tMArray = [tMArray t];
        end
        
      otherwise
        %do nothing
    end
     
   %% Start of standard UKF - generate the UKF sigma points.
    [root,p] = chol(N*Pukfs); %sqrt(n*P) 
    for i = 1 : N
        sigmas(:,i) = xhatukfs + root(i,:)'; %xhatsigma = xhat + xtilde
        sigmas(:,i+N) = xhatukfs - root(i,:)'; %xhatsigma = xhat - xtilde
    end
    for i = 1 : 2*N
        xbreve(:,i) = sigmas(:,i); %collect the sigmas into xbreve's
    end
    % Standard UKF time update (Eqn 14.59)
    for i = 1 : 2*N
        for tau = dtt : dtt : T
            u = (-Range_Accel_/tan(Pitch_)) - g;
            xbreve(1,i) = xbreve(1,i) + xbreve(2,i)*dtt + (u*dtt)/2;
            xbreve(2,i) = xbreve(2,i) + u*dtt;
            xbreve(:,i) = [xbreve(1,i);xbreve(2,i)];
        end
    end
    %Eqn 14.60
    xhatukfs = zeros(N,1);
    for i = 1 : 2*N
        xhatukfs = xhatukfs + Ws(i) * xbreve(:,i); %(Eqn. 14.60)
    end
    %Eqn 14.61
    Pukfs = Q;
    for i = 1 : 2*N
        Pukfs = Pukfs + Ws(i) * (xbreve(:,i) - xhatukfs) * (xbreve(:,i) - xhatukfs)';
    end
        
    % generate the UKF sigma points.
    [root,p] = chol(N*Pukfs); %sqrt(n*P) (Eqn 14.62)
    for i = 1 : N
        sigmas(:,i) = xhatukfs + root(i,:)'; %xhatsigma = xhat + xtilde
        sigmas(:,i+N) = xhatukfs - root(i,:)'; %xhatsigma = xhat - xtilde
    end
    for i = 1 : 2*N
        xbreve(:,i) = sigmas(:,i); %collect the sigmas into breve's
    end
    
     %Eqn 14.63
    for i = 1 : 2*N
        yukf = xbreve(1,i); %Range (ft)
        zukf(:,i) = ((y_meas + y_comp)/2 - yukf);
    end
    %Eqn 14.64
    zhat = 0;
    for i = 1 : 2*N
        zhat = zhat + Ws(i) * zukf(:,i);
    end
    %Eqn 14.65 and 66
    Py = R;
    Pxy = zeros(N,1);
    for i = 1 : 2*N
        Py = Py + Ws(i) * (zukf(:,i) - [zhat]) * (zukf(:,i) - [zhat])';
        Pxy = Pxy + Ws(i) * (xbreve(:,i) - xhatukfs) * (zukf(:,i) - [zhat])';
    end
    %Eqn 14.67
    Kukf = Pxy * inv(Py); 
%     xhatukfs = xhatukfs + Kukf * ((y_meas - y_comp)/2 - [zhat]); %follows sim prop more closely
    xhatukfs = xhatukfs + Kukf * ((y_meas - y_comp) - [zhat]); %follows measurements more closely
    Pukfs = Pukfs - Kukf * Py * Kukf';
end

%% Plots

p1 = figure(1);
subplot(2,1,1);
plot(tArray, xhatukfsArray(1,:),'k-o');
hold on
plot(tMArray, y_meaArr,'mo');
plot(tArray, xArray(1,:),'c-');
title(['Unscented Kalman Filter: ',method], 'FontSize', 12);
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('Est Altitude (ft)');
legend('xhat','meas','x','Location','SouthEast');
xlim([0,140]);

subplot(2,1,2);
plot(tArray, xhatukfsArray(2,:),'k-');
hold on
plot(tArray, xArray(2,:),'c-');
title(['Unscented Kalman Filter: ',method], 'FontSize', 12);
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('Est Altitude Rate (ft/sec)');
legend('xhat','x');
xlim([0,140]);

AltErr = std(xArray(1,:) - xhatukfsArray(1,:));
disp(['UKF - RMS Altitude estimation error = ', num2str(AltErr)]);

AltdotErr = std(xArray(2,:) - xhatukfsArray(2,:));
disp(['UKF - RMS Altitude Rate estimation error = ', num2str(AltdotErr)]);

% state estimate errors
p2 = figure(2);
subplot(2,1,1);
plot(tArray, xArray(1,:) - xhatukfsArray(1,:),'k-');
title(['Estimation Errors-- Unscented Kalman Filter: ',method]);
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('Altitude(ft)');
xlim([0,140]);
text(min(xlim),max(ylim)-300,['UKF - RMS Altitude estimation error = ', num2str(AltErr)]);

subplot(2,1,2);
plot(tArray, xArray(2,:) - xhatukfsArray(2,:),'k-');
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('Altitude Rate (ft/sec)');
text(min(xlim),max(ylim)-3,['UKF - RMS Altitude Rate estimation error = ', num2str(AltdotErr)]);


saveas(p1,['UKF_state_estimates_',method,'.png']);
saveas(p2,['UKF_state_estimate_errors_',method,'.png']);