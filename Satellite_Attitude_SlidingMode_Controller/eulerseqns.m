function dydt = eulerseqns(t,y)
% AE 6356 Assignment 6
% Solution by Glenn Lightsey
% November 12, 2024
%
global L % global variable 3x1 control torque
%
% define mass properties
%
Jx=100; Jy=90; Jz=120;
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
dydt = [    c1*w2*w3 + L(1)/Jx
            c2*w3*w1 + L(2)/Jy
            c3*w1*w2 + L(3)/Jz
            qd(1)
            qd(2)
            qd(3)
            qd(4)];

