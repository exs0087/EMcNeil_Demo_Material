function [latitude,longitude,altitude]=issorb_finalproj(t,theta_init)

% ISSORB -- find the latitude and longitude of a satellite in a circular
% orbit given time
%
% author: Erin McNeil

alpha = 90*pi/180;  %% rad, hard coded inclination
theta = theta_init*pi/180;  %% rad, hard coded initial latitude;
altitude = 200; %% km, hard coded altitude

mu = 4904; %km^3/s^2 = G*M
Rm = 1737; % km, Moon's radius
r = Rm + altitude; % satellite position

T_orbit = (2*pi/sqrt(mu))*r^(3/2); %time per 1 orbit (sec) , 7649 sec

theta_dot = 2*pi/T_orbit; %rad/sec

theta = theta + theta_dot*t;

x = r*cos(theta); ...
%     x = r*cos(theta)*cos(alpha); 
y = r*sin(theta); ...
%     y = r*sin(theta)*cos(alpha); 
z = 0;
% z = r*sin(theta)*sin(alpha);

% ellipsoid constants:
a = r; 
e=0;

% calculations:
b   = sqrt(a^2*(1-e^2));
ep  = sqrt((a^2-b^2)/b^2);
p   = sqrt(x.^2+y.^2);
th  = atan2(a*z,b*p);
lon = pi/2; ...
%     lon = atan2(y,x);
lat = atan2(y,x);
% lat = atan2((z+ep^2.*b.*sin(th).^3),(p-e^2.*a.*cos(th).^3));

% return lon in range [0,2*pi)
% lat = mod(lat,2*pi);

N   = a./sqrt(1-e^2.*sin(lat).^2);
alt = 200; ...alt = p./cos(lat)-N;

% correct for numerical instability in altitude near exact poles:
% (after this correction, error is about 2 millimeters, which is about
% the same as the numerical precision of the overall function)

k=abs(x)<1 & abs(y)<1;
alt(k) = abs(z(k))-b;

latitude = lat;
longitude = lon;
altitude = altitude;

return