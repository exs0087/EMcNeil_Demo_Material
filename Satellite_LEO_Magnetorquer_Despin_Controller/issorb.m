function [latitude,longitude,altitude]=issorb(t)

% ISSORB -- find the latitude and longitude of a satellite in a circular
% orbit given time
%
% author: Erin McNeil

alpha = 52*pi/180;  %% rad,  inclination
theta = -pi/2;  %% rad,  initial latitude;
altitude = 400; %% km, altitude

mu = 398600; %km^3/s^2 = G*M
Re = 6371; % km, Earth's radius
r = Re + altitude; % satellite position

T_orbit = (2*pi/sqrt(mu))*r^(3/2); %time per 1 orbit (sec)

theta_dot = 2*pi/T_orbit; %rad/sec

theta = theta + theta_dot*t;
A3_theta = [cos(theta) sin(theta) 0; ...
        -sin(theta) cos(theta) 0; ...
        0 0 1];
x = r*cos(theta)*cos(alpha); 
y = r*sin(theta)*cos(alpha); 
z = r*sin(theta)*sin(alpha);

% ellipsoid constants:
a = r; 
e=0;

% calculations:
b   = sqrt(a^2*(1-e^2));
ep  = sqrt((a^2-b^2)/b^2);
p   = sqrt(x.^2+y.^2);
th  = atan2(a*z,b*p);
lon = atan2(y,x);
lat = atan2((z+ep^2.*b.*sin(th).^3),(p-e^2.*a.*cos(th).^3));

% return lon in range [0,2*pi)
lon = mod(lon,2*pi);

N   = a./sqrt(1-e^2.*sin(lat).^2);
alt = p./cos(lat)-N;

% correct for numerical instability in altitude near exact poles:
% (after this correction, error is about 2 millimeters, which is about
% the same as the numerical precision of the overall function)

k=abs(x)<1 & abs(y)<1;
alt(k) = abs(z(k))-b;

latitude = lat;
longitude = lon;
altitude = altitude;

return