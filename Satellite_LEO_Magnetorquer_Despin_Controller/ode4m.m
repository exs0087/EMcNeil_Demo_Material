function Y = ode4m(odefun,tspan,y0,cmode,varargin)
% ODE4M modified ODE4 fixed step size integrator
% AE 6356 Spacecraft Attitude
%
% Solve differential equations with a non-adaptive method of order 4.
%   Y = ODE4(ODEFUN,TSPAN,Y0) with TSPAN = [T1, T2, T3, ... TN] integrates 
%   the system of differential equations y' = f(t,y) by stepping from T0 to 
%   T1 to TN. Function ODEFUN(T,Y) must return f(t,y) in a column vector.
%   The vector Y0 is the initial conditions at T0. Each row in the solution 
%   array Y corresponds to a time specified in TSPAN.
%
%   Y = ODE4(ODEFUN,TSPAN,Y0,P1,P2...) passes the additional parameters 
%   P1,P2... to the derivative function as ODEFUN(T,Y,P1,P2...). 
%
%   This is a non-adaptive solver. The step sequence is determined by TSPAN
%   but the derivative function ODEFUN is evaluated multiple times per step.
%   The solver implements the classical Runge-Kutta method of order 4.   
%
%   Example 
%         tspan = 0:0.1:20;
%         y = ode4(@vdp1,tspan,[2 0]);  
%         plot(tspan,y(:,1));
%     solves the system y' = vdp1(t,y) with a constant step size of 0.1, 
%     and plots the first component of the solution.   
%
%% set global variables
global L LHIST MHIST % global variables control torque vector and torque history

L=zeros(3,1); % L=[L1;L2;L3] total control torque per axis
%% check inputs
if ~isnumeric(tspan)
  error('TSPAN should be a vector of integration steps.');
end

if ~isnumeric(y0)
  error('Y0 should be a vector of initial conditions.');
end

h = diff(tspan);
if any(sign(h(1))*h <= 0)
  error('Entries of TSPAN are not in order.') 
end  

try
  f0 = feval(odefun,tspan(1),y0,varargin{:});
catch
  msg = ['Unable to evaluate the ODEFUN at t0,y0. ',lasterr];
  error(msg);  
end  

y0 = y0(:);   % Make a column vector.
if ~isequal(size(y0),size(f0))
  error('Inconsistent sizes of Y0 and f(t0,y0).');
end  

%% initialize params
neq = length(y0);
N = length(tspan);
Y = zeros(neq,N);
F = zeros(neq,4);
Y(:,1) = y0;

%% specify control mode
% cmode=2; % 1 = Bdot control, 2 = wxb control

%% initialize global variables
LHIST=zeros(3,N);   % history of control torque vectors
MHIST=zeros(3,N);   % history of magnetic dipole vectors

%% simulate Earth's magnetic field
Bx_ECI = zeros(16001,1);
By_ECI = zeros(16001,1);
Bz_ECI = zeros(16001,1);
B_ECI_mag = zeros(16001,1);
timey_wimey = zeros(16001,1);

stp = 0;
for dt = [0:1:16000]
    stp = stp + 1;
    
    [latitude,longitude,altitude]=issorb(dt); %rad, rad, km
    
    time = 739526;
    [Bx_North, By_East, Bz_Down] = igrf(time, latitude*180/pi, longitude*180/pi, altitude, 'geodetic');
    
    nedvec = [Bx_North,By_East,Bz_Down]; %nT, nT, nT
    ecivec = ned2eci(nedvec,latitude*180/pi,longitude*180/pi,dt);

    Bxi(stp,1) = ecivec(1); %nT
    Byi(stp,1) = ecivec(2); %nT
    Bzi(stp,1) = ecivec(3); %nT
    B_ECI_mag(stp,1) = sqrt(ecivec(1)^2 + ecivec(2)^2 + ecivec(3)^2);
    timey_wimey(stp,1) = dt;
end
%% specify static parameters
moment = 0.2;  % A m^2
sample_rate = 1;   % sec (1 Hz)

%% define mass properties
mass=10;     % mass in kg
height=0.34;   % height in m  (x)
width=0.2;   % width in m (y)
depth=0.1; % depth in m (z)
%
% rectangular prism spacecraft model
Jx=mass/12*(width^2+depth^2);
Jy=mass/12*(height^2+depth^2);
Jz=mass/12*(height^2+width^2);
%
% principal axis inertia tensor
J=[Jx,0,0;0,Jy,0;0,0,Jz];
%
% calculate coefficients
c1=(Jy-Jz)/Jx; d1=1/Jx;
c2=(Jz-Jx)/Jy; d2=1/Jy;
c3=(Jx-Jy)/Jz; d3=1/Jz;

%% controller
if cmode==1   % B-dot bang-bang control law
  % initialize parameters for calculation of B-dot
  % this section coded by Glenn Lightsey on February 11, 2022
  tlast=0; Bxlast=0; Bylast=0; Bzlast=0; mxlast=0; mylast=0; mzlast=0;
  %
  for i = 2:N
  %
    ti = tspan(i-1);
    hi = h(i-1);
    yi = Y(:,i-1);
  %
  % control algorithm goes here.
  % only calculate one control
  % distinguish between control sample steps and simulation steps
  %
    w1=yi(1); w2=yi(2); w3=yi(3);
    w=[w1;w2;w3];
    q1=yi(4); q2=yi(5); q3=yi(6); q4=yi(7);
    xi=[q4,-q3,q2;q3,q4,-q1;-q2,q1,q4;-q1,-q2,-q3];
    qd=0.5*xi*w;
  %
    C11=1-2*(q2^2+q3^2);
    C21=2*(q2*q1-q3*q4);
    C31=2*(q3*q1+q2*q4);
    C12=2*(q1*q2+q3*q4);
    C22=1-2*(q3^2+q1^2);
    C32=2*(q2*q3-q1*q4);
    C13=2*(q3*q1-q2*q4);
    C23=2*(q3*q2+q1*q4);
    C33=1-2*(q1^2+q2^2);
  %
    Bxeci1=Bxi(i-1);
    Byeci1=Byi(i-1);
    Bzeci1=Bzi(i-1);
    Bx_body=C11*Bxeci1+C12.*Byeci1+C13*Bzeci1;
    By_body=C21*Bxeci1+C22.*Byeci1+C23*Bzeci1;
    Bz_body=C31*Bxeci1+C32.*Byeci1+C33*Bzeci1;
  %
    Bvec_body=[Bx_body(1);By_body(1);Bz_body(1)]*1e-9; % units are in T
  %
    if (ti/sample_rate-floor(ti/sample_rate))<0.001, % test to see if this is a control update (1 msec test threshold)
      if (ti~=tlast),    % execute if this is not the first sample
        % update control
        % calculate new control
        Bdotx=(Bx_body(1)-Bxlast)/(ti-tlast); % units are in nanoT/s
        Bdoty=(By_body(1)-Bylast)/(ti-tlast);
        Bdotz=(Bz_body(1)-Bzlast)/(ti-tlast);
        mvec_body=-moment*[sign(Bdotx);sign(Bdoty);sign(Bdotz)]; % bang-bang control
        torque_body=cross(mvec_body,Bvec_body);
        mxlast=mvec_body(1);
        mylast=mvec_body(2);
        mzlast=mvec_body(3);
        % update BXxlast, Bylast, Bzlast, tlast for future reference
	    Bxlast=Bx_body(1);
	    Bylast=By_body(1);
	    Bzlast=Bz_body(1);
	    tlast=ti;
      else   % this is the first sample
        mvec_body=[0;0;0];     % zero control
	    torque_body=cross(mvec_body,Bvec_body);
        Bxlast=Bx_body(1);
        Bylast=By_body(1);
	    Bzlast=Bz_body(1);
	    tlast=ti;
      end
    else  % this is not a control sample
    % simulate with previous control
    % but torque is different, so need to calculate
      mvec_body=[mxlast;mylast;mzlast];
	  torque_body=cross(mvec_body,Bvec_body);
    end
    % pass torques to odefun
    % using global variables now, pass as function parameters possibly eventually
    L(1)=torque_body(1)*d1;
    L(2)=torque_body(2)*d2;
    L(3)=torque_body(3)*d3;
    LHIST(:,i)=L(1:3);   % store control torque in global history variable for plotting
    MHIST(:,i)=mvec_body(1:3); % store magnetic dipole in global variable for plotting
  %
  % integrator function calls happen below here
  % assumes control torque is the same, even though bfield changes at future times
  % if sample time is fast compared to change in bfield, should be okay.
  % to improve performance, calculate torque at each intermediate step
  % this is where you need interpolation step.
  %
    F(:,1) = feval(odefun,ti,yi,varargin{:});
    F(:,2) = feval(odefun,ti+0.5*hi,yi+0.5*hi*F(:,1),varargin{:});
    F(:,3) = feval(odefun,ti+0.5*hi,yi+0.5*hi*F(:,2),varargin{:});  
    F(:,4) = feval(odefun,tspan(i),yi+hi*F(:,3),varargin{:});
    Y(:,i) = yi + (hi/6)*(F(:,1) + 2*F(:,2) + 2*F(:,3) + F(:,4));
  end
end
%% this section coded by Erin McNeil on September 28, 2024
if cmode==2 % w cross b control
  % initialize parameters for calculation of B-dot
  % this section coded by Erin McNeil on September 28, 2024
  tlast=0; Bxlast=0; Bylast=0; Bzlast=0; mxlast=0; mylast=0; mzlast=0;
  %
  for i = 2:N
  %
    ti = tspan(i-1);
    hi = h(i-1);
    yi = Y(:,i-1);
  %
  % control algorithm goes here.
  % only calculate one control
  % distinguish between control sample steps and simulation steps
  %
    w1=yi(1); w2=yi(2); w3=yi(3);
    w=[w1;w2;w3];
    q1=yi(4); q2=yi(5); q3=yi(6); q4=yi(7);
    xi=[q4,-q3,q2;q3,q4,-q1;-q2,q1,q4;-q1,-q2,-q3];
    qd=0.5*xi*w;
  %
    C11=1-2*(q2^2+q3^2);
    C21=2*(q2*q1-q3*q4);
    C31=2*(q3*q1+q2*q4);
    C12=2*(q1*q2+q3*q4);
    C22=1-2*(q3^2+q1^2);
    C32=2*(q2*q3-q1*q4);
    C13=2*(q3*q1-q2*q4);
    C23=2*(q3*q2+q1*q4);
    C33=1-2*(q1^2+q2^2);
  %
    Bxeci1=Bxi(i-1);
    Byeci1=Byi(i-1);
    Bzeci1=Bzi(i-1);
    Bx_body=C11*Bxeci1+C12.*Byeci1+C13*Bzeci1;
    By_body=C21*Bxeci1+C22.*Byeci1+C23*Bzeci1;
    Bz_body=C31*Bxeci1+C32.*Byeci1+C33*Bzeci1;
  %
    Bvec_body=[Bx_body(1);By_body(1);Bz_body(1)]*1e-9; % units are in T
    
  %
    if (ti/sample_rate-floor(ti/sample_rate))<0.001 % test to see if this is a control update (1 msec test threshold)
      if (ti~=tlast)    % execute if this is not the first sample
        % update control
        % calculate new control
        Bdotx=(Bx_body(1)-Bxlast)/(ti-tlast); % units are in nanoT/s
        Bdoty=(By_body(1)-Bylast)/(ti-tlast);
        Bdotz=(Bz_body(1)-Bzlast)/(ti-tlast);

        B_body_norm = norm(Bvec_body);
        bx = Bx_body(1)/B_body_norm;
        by = By_body(1)/B_body_norm;
        bz = Bz_body(1)/B_body_norm;
        b = [bx;by;bz];
        k = 1;
        mvec_body = (k/B_body_norm)*cross(w,b);
        
        if mvec_body(1) > 0.2
            mvec_body(1) = 0.2;
        elseif mvec_body(1) < -0.2
            mvec_body(1) = -0.2;
        end
        
        if mvec_body(2) > 0.2
            mvec_body(2) = 0.2;
        elseif mvec_body(2) < -0.2
            mvec_body(2) = -0.2;
        end
        
        if mvec_body(3) > 0.2
            mvec_body(3) = 0.2;
        elseif mvec_body(3) < -0.2
            mvec_body(3) = -0.2;
        end
        
        torque_body=cross(mvec_body,Bvec_body);
        mxlast=mvec_body(1);
        mylast=mvec_body(2);
        mzlast=mvec_body(3);
        % update BXxlast, Bylast, Bzlast, tlast for future reference
	    Bxlast=Bx_body(1);
	    Bylast=By_body(1);
	    Bzlast=Bz_body(1);
	    tlast=ti;
      else   % this is the first sample
        mvec_body=[0;0;0];     % zero control
	    torque_body=cross(mvec_body,Bvec_body);
        Bxlast=Bx_body(1);
        Bylast=By_body(1);
	    Bzlast=Bz_body(1);
	    tlast=ti;
      end
    else  % this is not a control sample
    % simulate with previous control
    % but torque is different, so need to calculate
      mvec_body=[mxlast;mylast;mzlast];
	  torque_body=cross(mvec_body,Bvec_body);
    end
    % pass torques to odefun
    % using global variables now, pass as function parameters possibly eventually
    L(1)=torque_body(1)*d1;
    L(2)=torque_body(2)*d2;
    L(3)=torque_body(3)*d3;
    LHIST(:,i)=L(1:3);   % store control torque in global history variable for plotting
    MHIST(:,i)=mvec_body(1:3); % store magnetic dipole in global variable for plotting
  %
  % integrator function calls happen below here
  % assumes control torque is the same, even though bfield changes at future times
  % if sample time is fast compared to change in bfield, should be okay.
  % to improve performance, calculate torque at each intermediate step
  % this is where you need interpolation step.
  %
    F(:,1) = feval(odefun,ti,yi,varargin{:});
    F(:,2) = feval(odefun,ti+0.5*hi,yi+0.5*hi*F(:,1),varargin{:});
    F(:,3) = feval(odefun,ti+0.5*hi,yi+0.5*hi*F(:,2),varargin{:});  
    F(:,4) = feval(odefun,tspan(i),yi+hi*F(:,3),varargin{:});
    Y(:,i) = yi + (hi/6)*(F(:,1) + 2*F(:,2) + 2*F(:,3) + F(:,4));
  end
end
Y = Y.';
