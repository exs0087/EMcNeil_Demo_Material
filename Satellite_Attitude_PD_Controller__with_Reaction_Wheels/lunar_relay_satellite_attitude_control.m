clear;close all;clc;
%% Spacecraft Control -- lunar relay satellite
% Author: Erin McNeil
% Date: 11/23/24

format long

%% determine desired quaternion for s/c to point to Lunar Base while orbiting

altitude = 200; %% km, altitude
Rm = 1737; % km, Moon's radius
angle_w_LOS = acos(Rm/(Rm + altitude));  
% when latitude is between -63 and -90 we have LOS to base

mu = 4904; %km^3/s^2 = G*M
r = Rm + altitude; % satellite position
T_orbit = (2*pi/sqrt(mu))*r^(3/2); %time per 1 orbit (sec)

LOS = 1;
init_lat = -90;

%orbital elements
a = r; %km,
e = 0;
i = 90*(pi/180); %rad

M0 = 0*(pi/180);%rad, 
omega = 90*(pi/180); %rad
w = 90*(pi/180); 
t0 = 0; %sec
alt = 200; %% km, satellite altitude
Rm = 1737; % km, Moon's radius
n = sqrt(mu/a^3);

qc = zeros(4,1);
angle_p_sec = (126*(pi/180))/(1146); %rad/sec
alpha = asin(Rm/r); %rad 
stp = 0;
for t = -573:1:573 
    stp = stp + 1;
        
    M = M0 + n*(t - t0);
    E = M; 

    r = a; %radial distance in km
    alt(stp) = r - Rm; %altitude in km

    x = a*(cos(E) - e);
    y = a*sqrt(1 - e^2)*sin(E);
    
    x_plot(stp) = x;
    y_plot(stp) = y;
    z_plot(stp) = 0;

    xdot = -((n*a^2)/r)*sin(E);
    ydot = ((n*a^2)/r)*sqrt(1-e^2)*cos(E);

    A11 = cos(omega)*cos(w) - sin(omega)*sin(w)*cos(i);
    A12 = sin(omega)*cos(w) + cos(omega)*sin(w)*cos(i);
    A13 = sin(w)*sin(i);
    A21 = -cos(omega)*sin(w)-sin(omega)*cos(w)*cos(i);
    A22 = -sin(omega)*sin(w)+cos(omega)*cos(w)*cos(i);
    A23 = cos(w)*sin(i);
    A31 = sin(omega)*sin(i);
    A32 = -cos(omega)*sin(i);
    A33 = cos(i);

    r_comp = [A11 A21; A12 A22; A13 A23]*[x;y];%radial distance x,y,z components in km

    rdot_comp = [A11 A21; A12 A22; A13 A23]*[xdot;ydot];%radial velocity components in km/sec
    rdot = norm(rdot_comp);%radial velocity in km/sec
    rdotdot = [0;0;0]; 

    o3I = -cross(r_comp,rdot_comp)/norm(cross(r_comp,rdot_comp));
    o1I = -r_comp/norm(r_comp);
    o2I = cross(o1I,o3I);
    
    % z-axis rotation DCM
    alpha = alpha - angle_p_sec; %rad
    alpha_plot(stp) = alpha*(180/pi); %deg
    
    A3 = [cos(alpha) sin(alpha) 0;
          -sin(alpha) cos(alpha) 0;
          0 0 1]; 

    A_IO = [o1I o2I o3I]*A3; %desired attitude earth pointing LVLH aligned with RTN at t=0 with correction to point at base

    theta = acos((trace(A_IO) - 1)/2);

    e_ = (1/(2*sin(theta)))*[(A_IO(2,3) - A_IO(3,2)); ...
                            (A_IO(3,1) - A_IO(1,3)); ...
                            (A_IO(1,2) - A_IO(2,1))];

    qc(:,stp) = [(e_*sin(theta/2)); cos(theta/2)]; %desired attitude earth pointing LVLH
    
    wc_cur = [0; ...
                  -norm(cross(r_comp,rdot_comp))/norm(r_comp)^2; ...
                  (norm(r_comp)*dot(o2I,rdotdot))/norm(cross(r_comp,rdot_comp))];

    wc(:,stp) = wc_cur;
end

%% satellite PD controller 

        kp = 0.01;
        kd = 0.025;
        
        % PD controller
        cmode = 1;

        tspan = [-573:1:573]; % seconds
        w0 = [0; 0; 0]; % rad/s
        q0 = [0; 0; 0; 1];
        y0 = [w0; q0];
        [ysim] = ode4fp(@eulerseqns, tspan, y0, cmode,kp,kd);

        stp = 0;
        for tt = 1:1:1147
            stp = stp + 1;

            q1 = ysim(tt,4);
            q2 = ysim(tt,5);
            q3 = ysim(tt,6);
            q4 = ysim(tt,7);

            w = ysim(tt,1:3);

            Xi=[qc(4,tt),-qc(3,tt),qc(2,tt);qc(3,tt),qc(4,tt),-qc(1,tt);-qc(2,tt),qc(1,tt),qc(4,tt);-qc(1,tt),-qc(2,tt),-qc(3,tt)];

            dq3=Xi'*[q1;q2;q3;q4];
            dq4=[q1;q2;q3;q4]'*qc(:,tt);

        %     theta_err < 1 deg
            theta_err(stp) = 2*asin(sqrt(dq3(1)^2 + dq3(2)^2 + dq3(3)^2))*180/pi;

            w_err(stp) = sqrt((w(1) - wc(1,tt)').^2 + (w(2) - wc(2,tt)').^2 + (w(3) - wc(3,tt)').^2 )*180/pi;
        end
        theta_err_tot = sum(theta_err);
        w_err_tot = sum(w_err);

%% 
        q1 = ysim(1,4);
        q2 = ysim(1,5);
        q3 = ysim(1,6);
        q4 = ysim(1,7);

        Xi=[qc(4,1),-qc(3,1),qc(2,1);qc(3,1),qc(4,1),-qc(1,1);-qc(2,1),qc(1,1),qc(4,1);-qc(1,1),-qc(2,1),-qc(3,1)];

        dq3=Xi'*[q1;q2;q3;q4];
        dq4=[q1;q2;q3;q4]'*qc(:,1);
        
        dq = [dq3;dq4];
        
%% 
        T = -573:1:573;
        figure;
        plot(T,qc(1,:),'r'); hold on;
        plot(T,qc(2,:),'c')
        plot(T,qc(3,:),'m')
        plot(T,qc(4,:),'g')
        title('Erin McNeil: quaternion attitude vs time', 'FontSize', 12);
        ylabel('quaternions');
        xlabel('time (sec)');
        legend('qBI_1','qBI_2','qBI_3','qBI_4')
        ylim([-1 1]);
        set(gcf, 'Position', get(0, 'Screensize'));
        saveas(gcf,'quaternions.png')



%%         
        tspan = [-573:1:573]'; % seconds
        figure
        plot(tspan,theta_err)
        ylabel('theta error (deg)');
        xlabel('time (sec)');
        title('Erin McNeil: theta error vs time', 'FontSize', 12);
        ylim([0 2]);
        grid on
        set(gcf, 'Position', get(0, 'Screensize'));
        saveas(gcf,'thetaerr.png')

        global LHIST

        figure;
        plot(tspan,LHIST(1,:),'r-*');
        hold on;
        plot(tspan,LHIST(2,:),'b-');
        plot(tspan,LHIST(3,:),'m-');
        grid on
        % xlim([-573:1:573])
        grid on
        title('Erin McNeil: Torque vs time', 'FontSize', 12);
        legend('L1','L2','L3')
        xlabel('time (s)')
        ylabel('Torque')
        set(gcf, 'Position', get(0, 'Screensize'));
        saveas(gcf,'Torque.png')
