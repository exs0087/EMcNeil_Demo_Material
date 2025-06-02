clear all; close all;clc;
%% AE6505 Final Project -- Apollo 11 Final Descent: Particle Filter
% Use of an Particle Filter to improve estimates of the true altitude from unreliable
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
Q = eye(2)*1; %process noise

%% initial state  
x = [400; ... %altitude (ft)
    0.5;... %altitude rate (ft/sec)
    ]; % Initial state x0
xhatukfs = x;

y_meas = 400;

%% initialize particles
N = 1000; % number of particles
RoughParam = 0.05;
xhatPart = x;
% Initialize the particle filter
for j = 1 : N
    xpart(:,j) = x + sqrt(P) * [randn; randn];
end
randn('state',sum(100*clock)); % random number generator seed (investigate different randn settings)

%% Iterations
dt = 0.2;
tf = 150;

dtPlot = dt; % How often to plot results
tPlot = -inf;
tMeas = -inf;

% Initialize arrays for plotting at the end of the program
tArray = [];
xhatPartArr = [];
xArray = [];
y_meaArr = [];
tMArray = [];

stp = 0;
iter = 0;
for t = 0 : dt : tf
    iter = iter + 1;
    
    if t >= tPlot + dtPlot
        % Save data for plotting
        tPlot = t + dtPlot - eps;
        tArray = [tArray t];
        xArray = [xArray x];
        xhatPartArr = [xhatPartArr xhatPart]; 
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
     
    %% Particle filter
    xpartminus = xpart;
    for i = 1 : N
        xtemp = xpartminus(:,i);
        
        u = (-Range_Accel_/tan(Pitch_)) - g;
        xtemp(1) = xtemp(1) + xtemp(2)*dt + (u*dt)/2;
        xtemp(2) = xtemp(2) + u*dt;
        xpartminus(:,i) = [xtemp(1);xtemp(2)] + diag(sqrt(Q)*randn);
        
        y = xpartminus(1);
        ypart = y + sqrt(R)*randn;
        
        vhat(:,i) = (y_meas - y_comp) - ypart;
    end
    % Note that we need to scale all of the q(i) probabilities in a way
    % that does not change their relative magnitudes.
    % Otherwise all of the q(i) elements will be zero because of the
    % large value of the exponential.

    vhatscale = max(abs(vhat(1,i))) / 4;
    qsum = [0];
    for i = 1 : N
        q(1,i) = (1 / sqrt(abs(R(1,1))) / sqrt((2*pi)*N)) * exp(-(vhat(1,i)/vhatscale)^2 / 2 / R(1,1));
        qsum = qsum + q(:,i);
    end
    % Normalize the likelihood of each a priori estimate.
    for i = 1 : N
        q(1,i) = q(1,i) / qsum(1);
    end
    % Resample.
    for i = 1 : N
        u = rand; % uniform random number between 0 and 1
        qtempsum = [0];
        for j = 1 : N
            qtempsum(1) = qtempsum(1) + q(1,j);
            if qtempsum(1) >= u
                xpart(:,i) = xpartminus(:,j);
                % Use roughening to prevent sample impoverishment.
                E = max(xpartminus')' - min(xpartminus')';
                sigma = RoughParam * E * N^(-1/length(x));
                xpart(:,i) = xpart(:,i) + sigma .* [randn; randn];
                break;
            end
        end
    end
    % The particle filter estimate is the mean of the particles.
    xhatPart = 0;
    for i = 1 : N
        xhatPart = xhatPart + xpart(:,i);
    end
    xhatPart = xhatPart / N;
end

%% Plots

p1 = figure(1);
subplot(2,1,1);
plot(tArray, xhatPartArr(1,:),'k-o');
hold on
plot(tMArray, y_meaArr,'mo');
plot(tArray, xArray(1,:),'c-');
title(['Particle  Filter: ',method], 'FontSize', 12);
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('Est Altitude (ft)');
legend('xhat','meas','x','Location','SouthEast');
xlim([0,140]);

subplot(2,1,2);
plot(tArray, xhatPartArr(2,:),'k-');
hold on
plot(tArray, xArray(2,:),'c-');
title(['Particle  Filter: ',method], 'FontSize', 12);
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('Est Altitude Rate (ft/sec)');
legend('xhat','x');
xlim([0,140]);

AltErr = std(xArray(1,:) - xhatPartArr(1,:));
disp(['Particle - RMS Altitude estimation error = ', num2str(AltErr)]);

AltdotErr = std(xArray(2,:) - xhatPartArr(2,:));
disp(['Particle - RMS Altitude Rate estimation error = ', num2str(AltdotErr)]);

% state estimate errors
p2 = figure(2);
subplot(2,1,1);
plot(tArray, xArray(1,:) - xhatPartArr(1,:),'k-');
title(['Estimation Errors-- Particle  Filter: ',method]);
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('Altitude(ft)');
xlim([0,140]);
text(min(xlim),max(ylim)-300,['Particle - RMS Altitude estimation error = ', num2str(AltErr)]);

subplot(2,1,2);
plot(tArray, xArray(2,:) - xhatPartArr(2,:),'k-');
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('Altitude Rate (ft/sec)');
text(min(xlim),max(ylim)-3,['Particle - RMS Altitude Rate estimation error = ', num2str(AltdotErr)]);


saveas(p1,['Particle_state_estimates_',method,'.png']);
saveas(p2,['Particle_state_estimate_errors_',method,'.png']);