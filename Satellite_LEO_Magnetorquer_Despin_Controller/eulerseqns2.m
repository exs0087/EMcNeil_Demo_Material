function dydt = eulerseqns2(t,y)
% AE 6356 Spacecraft Attitude
% Assignment 3 Problem 3
% Glenn Lightsey Fall 2024
%
global L % global variable 3x1 control torque
%
% define mass properties
mass=10;     % mass in kg
height=0.34;   % height in m  (x)
width=0.2;   % width in m (y)
depth=0.1;  % depth in m (z)
%
% rectangular prism spacecraft model
Jx=mass/12*(width^2+depth^2);
Jy=mass/12*(height^2+depth^2);
Jz=mass/12*(height^2+width^2);
%
% principal axis inertia tensor
J=[Jx,0,0;0,Jy,0;0,0,Jz];
%
% calculate inertia coefficients
c1=(Jy-Jz)/Jx; 
c2=(Jz-Jx)/Jy; 
c3=(Jx-Jy)/Jz; 
%
w1=y(1); w2=y(2); w3=y(3);
w=[w1;w2;w3];
q1=y(4); q2= y(5); q3=y(6); q4=y(7);
xi=[q4,-q3,q2;q3,q4,-q1;-q2,q1,q4;-q1,-q2,-q3];
qd=0.5*xi*w;
%
dydt = [    c1*w2*w3 + L(1)
            c2*w3*w1 + L(2)
            c3*w1*w2 + L(3)
            qd(1)
            qd(2)
            qd(3)
            qd(4)];

