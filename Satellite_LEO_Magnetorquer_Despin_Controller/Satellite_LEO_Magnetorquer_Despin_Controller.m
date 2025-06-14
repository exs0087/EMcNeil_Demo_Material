clear;close all;clc;
%% AE 6356 - Satellite in LEO Magnetorquer Despin Controller

% Author: Erin McNeil
% Date: 09/23/24

format longG
%% Bdot control
cmode = 1;
tspan = [0:1:16000];
w0 = [0; 0.2; 0.2]; % rad/s
q0 = [0; 0; 0; 1];
y0 = [w0; q0];

[ysim] = ode4m(@eulerseqns2, tspan, y0,cmode);
tsim = tspan';


figure; 
plot(tsim,ysim(:,4),'r-');
hold on;
plot(tsim,ysim(:,5),'g-');
plot(tsim,ysim(:,6),'b-');
plot(tsim,ysim(:,7),'m-');
grid on
xlim([0 16000])
ylim([-1 1])
legend('q1','q2','q3','q4')
title({'Erin McNeil:','Quaternion Plot'});
xlabel('time (s)')
ylabel('quaternion values')
saveas(gcf,'quaternion_plot_cmode_1','png')

figure; 
plot(tsim,ysim(:,1),'r-*');
hold on;
plot(tsim,ysim(:,2),'b-');
plot(tsim,ysim(:,3),'m-');
grid on
xlim([0 16000])
ylim([-0.004 0.004])
legend('w1','w2','w3')
title({'Erin McNeil:','Angular Velocity Plot'});
xlabel('time (s)')
ylabel('angular velocity (rad/sec)')
saveas(gcf,'angular_velocity_plot_cmode_1','png')

mag_w = zeros(length(tsim),1);
for ii = 1:length(tsim)
    w_norm = sqrt(ysim(ii,1)^2 + ysim(ii,2)^2 + ysim(ii,3)^2);
    mag_w(ii,1) = w_norm; 
end

figure; 
plot(tsim,mag_w,'m-');
grid on
xlim([0 16000])
legend('magnitude of angular velocity')
title({'Erin McNeil:', 'Magnitude of Angular Velocity Plot'});
xlabel('time (s)')
ylabel('magnitude of angular velocity (rad/sec)')
saveas(gcf,'magnitude_angular_velocity_cmode_1','png')

mu = 398600; %km^3/s^2 = G*M
Re = 6371; % km, Earth's radius
altitude = 400; %% km, hard coded altitude
r = Re + altitude; % satellite position

T_orbit = (2*pi/sqrt(mu))*r^(3/2); %time per 1 orbit (sec)

BangBang_final_value = mag_w(end)*(T_orbit/2*pi); %rad/sec * (1)rev/(2*pi)rad * (#)sec/(1)orbit

global MHIST
figure; 
plot(tsim',MHIST(1,:),'r-*');
hold on;
plot(tsim',MHIST(2,:),'b-');
plot(tsim',MHIST(3,:),'m-');
grid on
xlim([0 16000])
% ylim([-0.004 0.004])
legend('mb1','mb2','mb3')
title({'Erin McNeil:','Magnetic Dipole Moment Plot'});
xlabel('time (s)')
ylabel('Magnetic Dipole Moment (A*m^2)')
saveas(gcf,'magnetic_dipole_moment_plot_cmode_1','png')



%% wxb control
cmode = 2;
tspan = [0:1:16000];
w0 = [0; 0.2; 0.2]; % rad/s
q0 = [0; 0; 0; 1];
y0 = [w0; q0];

[ysim] = ode4m(@eulerseqns2, tspan, y0,cmode);
tsim = tspan';


figure; 
plot(tsim,ysim(:,4),'r-');
hold on;
plot(tsim,ysim(:,5),'g-');
plot(tsim,ysim(:,6),'b-');
plot(tsim,ysim(:,7),'m-');
grid on
xlim([0 16000])
ylim([-1 1])
legend('q1','q2','q3','q4')
title({'Erin McNeil:','Quaternion Plot'});
xlabel('time (s)')
ylabel('quaternion values')
saveas(gcf,'quaternion_plot_cmode_2','png')

figure; 
plot(tsim,ysim(:,1),'r-*');
hold on;
plot(tsim,ysim(:,2),'b-');
plot(tsim,ysim(:,3),'m-');
grid on
xlim([0 16000])
ylim([-0.004 0.004])
legend('w1','w2','w3')
title({'Erin McNeil:','Angular Velocity Plot'});
xlabel('time (s)')
ylabel('angular velocity (rad/sec)')
saveas(gcf,'angular_velocity_plot_cmode_2','png')

mag_w = zeros(length(tsim),1);
for ii = 1:length(tsim)
    w_norm = sqrt(ysim(ii,1)^2 + ysim(ii,2)^2 + ysim(ii,3)^2);
    mag_w(ii,1) = w_norm; 
end

figure; 
plot(tsim,mag_w,'m-');
grid on
xlim([0 16000])
legend('magnitude of angular velocity')
title({'Erin McNeil:', 'Magnitude of Angular Velocity Plot'});
xlabel('time (s)')
ylabel('magnitude of angular velocity (rad/sec)')
saveas(gcf,'magnitude_angular_velocity_cmode_2','png')

mu = 398600; %km^3/s^2 = G*M
Re = 6371; % km, Earth's radius
altitude = 400; %% km, hard coded altitude
r = Re + altitude; % satellite position

T_orbit = (2*pi/sqrt(mu))*r^(3/2); %time per 1 orbit (sec)

wcrossb_final_value = mag_w(end)*(T_orbit/2*pi); %rad/sec * (1)rev/(2*pi)rad * (#)sec/(1)orbit

global MHIST
figure; 
plot(tsim',MHIST(1,:),'r-*');
hold on;
plot(tsim',MHIST(2,:),'b-');
plot(tsim',MHIST(3,:),'m-');
grid on
xlim([0 16000])
% ylim([-0.004 0.004])
legend('mb1','mb2','mb3')
title({'Erin McNeil:','Magnetic Dipole Moment Plot'});
xlabel('time (s)')
ylabel('Magnetic Dipole Moment (A*m^2)')
saveas(gcf,'magnetic_dipole_moment_plot_cmode_2','png')

