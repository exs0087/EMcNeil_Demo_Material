clear;close all;clc;
%% MEKF filter -- Spacecraft Attitude Estimation with Gyro Bias Compensation
% Author: Erin McNeil
% Date: 10/26/24

format long
%%
load('A5P7.mat') % load simulated gyro measurement data

%%  
sigma_theta = 5; %deg
sigma_Beta = 0.1/3600; %deg/sec
P_1 = (sigma_theta^2)*eye(3);
P_2 = (sigma_Beta^2)*eye(3);
P = blkdiag(P_1, P_2);

dt = 1;
sigma_v = 0.1/3600; %deg/sec
sigma_u = 0.1*(1/3600); %deg/sec standard deviation in all three axes
Q_1 = ((sigma_v^2)*dt + (1/3)*(sigma_u^2)*dt^3)*eye(3);
Q_2 = -(0.5*(sigma_u^2)*dt^2)*eye(3);
Q_3 = -(0.5*(sigma_u^2)*dt^2)*eye(3);
Q_4 = (sigma_v^2)*dt*eye(3);
Q = [Q_1 Q_2; 
     Q_3 Q_4];
  
sigma_M = 10; %deg
sigma_S = 1; %deg
R_1 = (sigma_M^2)*eye(3);
R_2 = (sigma_S^2)*eye(3);
R = blkdiag(R_1, R_2);

%% initial state  
q = [0.010293885233473;...
      -0.706932213120670; ...
      -0.017713746763756; ...
      -0.706984515498537]; % Initial state attitude
 
delta_x = [0;0;0;0;0;0];
dtheta = [0;0;0];
Beta_hat = [0;0;0];

sigx = 0;
sigy = 0;
sigz = 0;
sig_err_3a = 0;

%% MEKF Iterations
dt = 1;
tf = 3600-1;

dtPlot = dt; % How often to plot results
tPlot = -inf;

% Initialize arrays for plotting at the end of the program
delta_xArray = [];
tArray = [];
qArray = [];
KArray = [];
PArray = [];
BetaArray = [];
sigxArray = [];
sigyArray = [];
sigzArray = [];
sig_err_3aArray = [];

iter = 0;
for t = 0 : dt : tf
    iter = iter + 1;
    
    if t >= tPlot + dtPlot
        % Save data for plotting
        tPlot = t + dtPlot - eps;
        delta_xArray = [delta_xArray delta_x];
        qArray = [qArray q];
        PArray = [PArray P];
        BetaArray = [BetaArray Beta_hat];
        sigxArray = [sigxArray sigx];
        sigyArray = [sigyArray sigy];
        sigzArray = [sigzArray sigz];
        sig_err_3aArray = [sig_err_3aArray sig_err_3a];
        tArray = [tArray t];
    end
    
    % reset delta_x
    delta_x = [0;0;0;0;0;0]; %deg
    
    % gyro meas update
    w_meas = wm(:,iter)*180/pi; %deg/sec
    w_meas_hat = w_meas - Beta_hat;
    
    Phi_meas = (sin(0.5*norm(w_meas_hat)*dt)*w_meas_hat)/norm(w_meas_hat);

    Phi_meas_cross = [0 -Phi_meas(3) Phi_meas(2);
                      Phi_meas(3) 0 -Phi_meas(1);
                      -Phi_meas(2) Phi_meas(1) 0];

    theta_of_wm = [(cos(0.5*norm(w_meas_hat)*dt)*eye(3) - Phi_meas_cross) (Phi_meas); ...
                           (-Phi_meas') (cos(0.5*norm(w_meas_hat)*dt))];
    q = theta_of_wm*q;
    
    A_q = [(q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2) (2*(q(1)*q(2) + q(3)*q(4))) (2*(q(1)*q(3) - q(2)*q(4))); ...
         (2*(q(2)*q(1) - q(3)*q(4))) (-q(1)^2 + q(2)^2 - q(3)^2 + q(4)^2) (2*(q(2)*q(3) + q(1)*q(4))); ...
         (2*(q(3)*q(1) + q(2)*q(4))) (2*(q(3)*q(2) - q(1)*q(4))) (-q(1)^2 - q(2)^2 + q(3)^2 + q(4)^2)];


    F = -[0 -w_meas_hat(3) w_meas_hat(2) -1     0     0;
          w_meas_hat(3) 0 -w_meas_hat(1)  0    -1     0;
         -w_meas_hat(2) w_meas_hat(1) 0   0     0    -1;
          0 0 0                           0     0     0;
          0 0 0                           0     0     0;
          0 0 0                           0     0     0];

    wx = -[0 -w_meas(3) w_meas(2);
          w_meas(3) 0 -w_meas(1);
          -w_meas(2) w_meas(1) 0];
      
    phi_11 = eye(3) - wx*(sin(norm(wx)*dt)/norm(wx)) + (wx^2)*((1-cos(norm(wx)*dt))/norm(wx)^2);
    phi_12 = wx*((1-cos(norm(wx)*dt))/norm(wx)^2) - eye(3)*dt - (wx^2)*((norm(wx)*dt - sin(norm(wx)*dt))/norm(wx)^3);
    phi_21 = eye(3)*0;
    phi_22 = eye(3);
    phi = [phi_11 phi_12;
           phi_21 phi_22];%state transition matrix
       
    Y = eye(6); 
    
    P = phi*P'*phi' + Y*Q*Y'; % a priori state estimate covariance (propogate)
    delta_x = phi*delta_x; %deg => zero always
    
%   gather measurements
    b1 = b1m(:,iter);
    b2 = b2m(:,iter);
        
    h_q = [A_q*r1; ...
           A_q*r2];

    b1_cross = [0 -b1(3) b1(2);
                b1(3) 0 -b1(1);
               -b1(2) b1(1) 0];

    b2_cross = [0 -b2(3) b2(2);
                b2(3) 0 -b2(1);
               -b2(2) b2(1) 0];
   zero_matrix = eye(3)*0;

    H_q = [b1_cross zero_matrix; ...
            b2_cross zero_matrix];


    K = P * H_q' * inv(H_q * P * H_q' + R);  % calculate kalman gain
    
    y = [b1; b2];
    delta_x = delta_x + K*(y - h_q); %calculate posteriori state estimate (update)
    dtheta = delta_x(1:3); %deg
    dBeta = delta_x(4:6); %deg

    P = (eye(6) - K * H_q) * P; %calculate posteriori state est covariance (update)

    q1v = dtheta/2;
    q1s = 1;  
    q1 = [q1v;
          q1s];

    q2v = q(1:3);
    q2s = q(4);  
    q2 = [q2v;
          q2s];  

     qstar = [(q1v*q2s + q2v*q1s - cross(q1v,q2v)); ...
               (q1s*q2s - dot(q1v,q2v))];

    q = qstar/norm(qstar);
    
    Beta_hat = Beta_hat + dBeta;
    
% calcs for plots
    sigx = sqrt(P(1,1));
    sigy = sqrt(P(2,2));
    sigz = sqrt(P(3,3));
    
    sig_err_3a = sqrt(sigx^2 + sigy^2 + sigz^2);
    
end

%% part b

figure(1);
subplot(4,1,1);
plot(tArray, qArray(1,:),'c-');
sgtitle('Erin McNeil:', 'FontSize', 12);
title('Spacecraft Quaternion Estimation', 'FontSize', 12);
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('Est q1');xlabel('Secs');  

subplot(4,1,2);
plot(tArray, qArray(2,:),'c-');
title('Spacecraft Quaternion Estimation', 'FontSize', 12);
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('Est q2');xlabel('Secs'); 

subplot(4,1,3);
plot(tArray, qArray(3,:),'c-');
title('Spacecraft Quaternion Estimation', 'FontSize', 12);
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('Est q3');xlabel('Secs'); 

subplot(4,1,4);
plot(tArray, qArray(4,:),'c-');
title('Spacecraft Quaternion Estimation', 'FontSize', 12);
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('Est q4');xlabel('Secs'); 
set(gcf, 'Position', get(0, 'Screensize'));

saveas(gcf,'SC_Quaternion_Estimations.png');

%% part c

figure(2);
subplot(3,1,1);
plot(tArray, BetaArray(1,:),'c-');
sgtitle('Erin McNeil:', 'FontSize', 12);
title('Spacecraft Bias Estimation', 'FontSize', 12);
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('Est Beta x (deg)');xlabel('Secs');  

subplot(3,1,2);
plot(tArray, BetaArray(2,:),'c-');
title('Spacecraft Bias Estimation', 'FontSize', 12);
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('Est Beta y (deg)');xlabel('Secs'); 

subplot(3,1,3);
plot(tArray, BetaArray(3,:),'c-');
title('Spacecraft Bias Estimation', 'FontSize', 12);
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('Est Beta z (deg)');xlabel('Secs'); 
set(gcf, 'Position', get(0, 'Screensize'));

saveas(gcf,'SC_Gyro_Bias_Estimations.png');

