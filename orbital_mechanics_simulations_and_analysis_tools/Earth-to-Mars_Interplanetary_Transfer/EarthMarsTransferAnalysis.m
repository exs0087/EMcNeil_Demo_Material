% =========================================================================
% Earth–Mars Transfer Analysis (Orbiter & Lander)
%
% Description:
%   (a) Solves Lambert's problem via universal variables for an Earth–Mars transfer.
%   (b) Computes departure injection ΔV from 200 km LEO and hyperbolic excess speeds.
%   (c) Determines lander entry speed and miss-distance precision for a desired flight-path angle.
%   (d) Computes orbiter Mars Orbit Insertion burn ΔV and targeting precision for circular orbit.
%
% Author: Erin McNeil
% Last updated: 2025-06-10
% =========================================================================

clear; clc; close all;
format longG;

%% --- Constants & Unit Conversions ---
AU      = 1.496e8;            % km
muSun   = 1.32712440018e11;   % km^3/s^2
muE     = 3.986e5;            % km^3/s^2
muM     = 4.282837e4;         % km^3/s^2
Re      = 6378;               % km
RM      = 3396;               % km

%% --- Mission Parameters ---
LEO_alt = 200;                       
r_LEO   = Re + LEO_alt;       % km
r_entry = 3522.2;             % km (Mars entry interface)
orb_alt = 400;                
r_orb   = RM + orb_alt;       % km

%% --- Heliocentric Vectors & Transfer Time ---
r1 = [0.707, -0.707, 0] * AU;
r2 = [0, 1.6655, 0]      * AU;
T  = 210 * 24 * 3600;         % s

%% Part (a): Lambert Solution
[v_dep, v_arr] = lambertSolver(r1, r2, T, muSun);
fprintf('Part (a): v_dep = [%.3f %.3f %.3f], v_arr = [%.3f %.3f %.3f] km/s\n', ...
    v_dep, v_arr);

%% Part (b): Departure ΔV & v∞
v_inf_dep = abs(norm(v_dep) - 29.78);                    
v_LEO     = sqrt(muE / r_LEO);                          
dV_dep    = sqrt(v_inf_dep^2 + v_LEO^2) - v_LEO;        
v_inf_arr = abs(norm(v_arr) - 24.077);                  
fprintf('Part (b): v∞_dep=%.3f, ΔV_dep=%.3f, v∞_arr=%.3f km/s\n', ...
    v_inf_dep, dV_dep, v_inf_arr);

%% Part (c): Lander Entry & Precision
v_entry    = sqrt(v_inf_arr^2 + 2*muM/r_entry);
gamma_tgt  = -11.5*pi/180;        
computeG    = @(rp) flightPathAngle(rp, v_inf_arr, muM, r_entry);
rp_tgt     = fzero(@(rp) computeG(rp)-gamma_tgt, [r_entry-200, r_entry-1]);
dgamma_dr  = derivative(computeG, rp_tgt);
tol_gamma  = 0.5*pi/180;         
miss_prec_L= tol_gamma/abs(dgamma_dr);
fprintf('Part (c): v_entry=%.3f km/s, rp=%.3f km, miss_prec=±%.3f km\n', ...
    v_entry, rp_tgt, miss_prec_L);

%% Part (d): Orbiter MOI & Precision
v_circ    = sqrt(muM / r_orb);     
v_peri    = @(rp) sqrt(v_inf_arr^2 + 2*muM/rp);
dV_MOI    = @(rp) abs(v_peri(rp)-v_circ);
d_dVdr    = derivative(dV_MOI, r_orb);
dvcirc_dr = 0.5*sqrt(muM)/r_orb^(3/2);
allowedDV = dvcirc_dr*4;      
miss_prec_O = allowedDV/abs(d_dVdr);
fprintf('Part (d): ΔV_MOI=%.3f km/s, miss_prec=±%.3f km\n', ...
    dV_MOI(r_orb), miss_prec_O);

%% --- Local Functions ---

function gamma = flightPathAngle(rp, v_inf, mu, r_e)
    e   = 1 + (rp*v_inf^2)/mu;
    h   = rp*sqrt(v_inf^2+2*mu/rp);
    cosn= (h^2/(mu*r_e)-1)/e;
    nu  = acos(max(min(cosn,1),-1));
    gamma = -atan((e*sin(nu))/(1+e*cos(nu)));
end

function val = derivative(f, x)
    dx = 1e-3; val = (f(x+dx)-f(x-dx))/(2*dx);
end

function [v1, v2] = lambertSolver(r1, r2, T, mu)
    tol=1e-6; maxit=200;
    r1n=norm(r1); r2n=norm(r2);
    th = acos(dot(r1,r2)/(r1n*r2n));
    A  = sin(th)*sqrt(r1n*r2n/(1-cos(th)));
    z=0; ratio=1; it=0;
    while abs(ratio)>tol && it<maxit
        [C,S]=stumpff(z);
        y=r1n+r2n+A*(z*S-1)/sqrt(C);
        x=sqrt(y/C);
        Tz=(x^3*S+A*sqrt(y))/sqrt(mu);
        dTdz=(Tz-T)/1; % approximate
        ratio=(Tz-T)/dTdz; z=z-ratio; it=it+1;
    end
    [C,S]=stumpff(z); y=r1n+r2n+A*(z*S-1)/sqrt(C);
    f=1-y/r1n; g=A*sqrt(y/mu); gdot=1-y/r2n;
    v1=(r2-f*r1)/g; v2=(gdot*r2-r1)/g;
end

function [C,S] = stumpff(z)
    if z>1e-6
        C=(1-cos(sqrt(z)))/z;
        S=(sqrt(z)-sin(sqrt(z)))/(z*sqrt(z));
    elseif z<-1e-6
        C=(1-cosh(sqrt(-z)))/z;
        S=(sinh(sqrt(-z))-sqrt(-z))/(-z*sqrt(-z));
    else
        C=1/2; S=1/6;
    end
end
