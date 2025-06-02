function ecivec = ned2eci(nedvec,latitude,longitude,time)
%% NED to ECEF
BxNorth = nedvec(1); %nT
ByEast = nedvec(2); %nT
BzDown = nedvec(3); %nT

lat0 = latitude; %deg
lon0 = longitude; %deg

%longitudinal rotation about z-axis, latitudinal rotation about y-axis
Bx_ECEF = cosd(lon0) .* (cosd(lat0) .* -BzDown) - (sind(lat0) .* BxNorth) - sind(lon0) .* ByEast;
By_ECEF = sind(lon0) .* (cosd(lat0) .* -BzDown) - (sind(lat0) .* BxNorth) + cosd(lon0) .* ByEast;
Bz_ECEF = (sind(lat0) .* -BzDown) + (cosd(lat0) .* BxNorth);

%% ECEF to ECI 
omega_dot  = 360/86164; % degrees/second (sidereal rotation rate of the Earth)
omega = omega_dot*time; %deg

%rotation omega about z-axis 
A3_omega = [cosd(omega) sind(omega) 0; -sind(omega) cosd(omega) 0; 0 0 1];

ecivec = A3_omega*[Bx_ECEF; By_ECEF; Bz_ECEF];

end