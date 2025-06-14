clear all; close all;clc;
%% AE6505 Final Project -- Apollo 11 Final Descent: Extended Kalman Filter
% Use of an EKF to improve estimates of the true altitude from unreliable
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
xhat = x;
y_meas = 400;

%% Iterations
dt = 0.2;
tf = 150;

dtPlot = dt; % How often to plot results
tPlot = -inf;
tMeas = -inf;


% Initialize arrays for plotting at the end of the program
tArray = [];
xhatArr = [];
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
        xhatArr = [xhatArr xhat]; 
    end
    
    Range_Accel_ = RangeAccel(iter);
    Pitch_ = Pitch(iter);
    
    u = (-Range_Accel_/tan(Pitch_)) - g;
    x(1) = x(1) + x(2)*dt + (u*dt)/2;
    x(2) = x(2) + u*dt;
    x = [x(1);x(2)] + diag(sqrt(Q)*randn);
    
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
        
%% Extended Kalman Filter
   
    A = [0 1; 0 0];
    H = [1 0];
    
    
    u = (-Range_Accel_/tan(Pitch_)) - g;
    xhat(1) = xhat(1) + xhat(2)*dt + (u*dt)/2;
    xhat(2) = xhat(2) + u*dt;
    xhat = [xhat(1);xhat(2)];
    
    y_comp = xhat(1); %Altitude (ft)
    
    Pdot = A*P + P*A' + Q;
    P = P + P*dt;

    K = P * H' * inv(H * P * H' + R);  
    xhat = xhat + K * (y_meas - y_comp); 
    P = (eye(2) - K*H)*P*(eye(2) - K*H)' + K*R*K'; 

end

%% Plots

p1 = figure(1);
subplot(2,1,1);
plot(tArray, xhatArr(1,:),'k-o');
hold on
plot(tMArray, y_meaArr,'mo');
plot(tArray, xArray(1,:),'c-');
title(['Extended Kalman Filter: ',method], 'FontSize', 12);
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('Est Altitude (ft)');
legend('xhat','meas','x','Location','SouthEast');
xlim([0,140]);

subplot(2,1,2);
plot(tArray, xhatArr(2,:),'k-');
hold on
plot(tArray, xArray(2,:),'c-');
title(['Extended Kalman Filter: ',method], 'FontSize', 12);
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('Est Altitude Rate (ft/sec)');
legend('xhat','x');
xlim([0,140]);

AltErr = std(xArray(1,:) - xhatArr(1,:));
disp(['EKF - RMS Altitude estimation error = ', num2str(AltErr)]);

AltdotErr = std(xArray(2,:) - xhatArr(2,:));
disp(['EKF - RMS Altitude Rate estimation error = ', num2str(AltdotErr)]);

% state estimate errors
p2 = figure(2);
subplot(2,1,1);
plot(tArray, xArray(1,:) - xhatArr(1,:),'k-');
title(['Estimation Errors-- Extended Kalman Filter: ',method]);
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('Altitude(ft)');
xlim([0,140]);
text(min(xlim),max(ylim)-300,['EKF - RMS Altitude estimation error = ', num2str(AltErr)]);

subplot(2,1,2);
plot(tArray, xArray(2,:) - xhatArr(2,:),'k-');
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('Altitude Rate (ft/sec)');
text(min(xlim),max(ylim)-3,['EKF - RMS Altitude Rate estimation error = ', num2str(AltdotErr)]);


saveas(p1,['EKF_state_estimates_',method,'.png']);
saveas(p2,['EKF_state_estimate_errors_',method,'.png']);