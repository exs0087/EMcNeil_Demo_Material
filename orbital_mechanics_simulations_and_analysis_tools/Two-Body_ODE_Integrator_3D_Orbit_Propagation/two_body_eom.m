function dydt = two_body_eom(t,y)
% AE 6356 Spacecraft Attitude
% Assignment 1 Problem 4
% by Erin McNeil, 1/2025
%
% define params
mu = 3.986e5; %km^3/s^2 (Earth)

% define EOM
r1=y(1); r2=y(2); r3=y(3);
r=[r1;r2;r3];
a=(-mu/norm(r)^3)*r;
v1=y(4); v2= y(5); v3=y(6);
v=[v1;v2;v3];
dydt = [    v(1)
            v(2)
            v(3)
            a(1) 
            a(2)
            a(3)];