function V_dot = Nonlinear_dynamics(airship, MAF)
 u = airship(1);
 v = airship(2);
 w = airship(3);
 p = airship(4);
 q = airship(5);
 r = airship(6);
 phi  = airship(7);
 theta = airship(8);
 psi = airship(9);
 
 XE = airship(10);
 YE = airship(11);
altitude = -airship(12); %linspace(0,12000, 500); %Altitude from 0 to 12000 m
%time = = airship(17);

%% Constants
g = 9.81;
m   = 5942; % mass
VB = 5153; %volume (m^3)
S = VB^(2/3) ; %surface area (m^2)

XB = 0; %x-origin point, m (boyancy)
ZB = 0; % z-origin point, m
XG = 7; %0.0 as given in [2], m
ZG = 5.1816; %m
rGxm = [XG-XB 0 ZG-ZB]';

L = 52; %Hull length
D = 14; %Hull max diameter
R = 5.6522; % Hull cross section radius at fins (m)
b = 14.13; %Hull centre line to fin tip length
delta = 0; %Elevator deflection 5 degrees
AR = 0.54;% Aspect ratio

keta = ketaparameters(L, D, R, b, AR, delta);

 % k1 = 0.08;
 % k2 = 0.89;
 % k_dash = 0.5;
 k1 = keta(1);
 k2 = keta(2);
 k_dash = keta(3);

etaf = 0.8;     % Example: Flap effectiveness
k44 = 0.1;      % Example: Roll damping coefficient

%eta = keta(1,4);
%etad = keta(1,5);
%etaf = keta(1,6);
%k44 = keta(1,7);
%k3D = keta(1,8);

Ixx = 200658; %moment of inerta Nm2
Iyy = 850900; %moment of inerta Nm2
Izz = 649699; %moment of inerta Nm2
Ixy =0;
Ixz =0;
Iyz =0;

%% Inertial forces (N) and moments (N-m)
I33 = eye(3);
rGx = [0 -rGxm(3) rGxm(2) ; rGxm(3) 0 -rGxm(1) ; -rGxm(2) rGxm(1) 0];
J = [ Ixx, Ixy, Ixz ; Ixy, Iyy, Iyz ; Ixz, Iyz, Izz];
Mrigid = [m*I33, -m*rGx; m*rGx,  J];
g_cap = [-sin(theta), cos(theta)*sin(phi), cos(theta)*cos(phi)]';
omg = [p q r]';
omgx = [0 -omg(3) omg(2) ; omg(3) 0 -omg(1) ; -omg(2) omg(1) 0];
V = [u v w]';

tau_I = [-m*omgx*V + m*omgx*rGx*omg ; -m*rGx*omgx*V - omgx*J*omg];


%% Gravitational forces and moments

%Atmosphere model
rho_0 = 1.225; %Density at sea level in kg/m3
H = 8400; %Scale height in meters;

%alitude range
%altitude = -airship(12); %linspace(0,12000, 500); %Altitude from 0 to 12000 m

density = rho_0 *exp(-altitude / H );
rho = density;
w = m*g;
B = VB*rho*g;
FG = w*g_cap; % Weight component in body frame (N)
tau_G = [FG; rGx*FG];

%% Aerostatic forces and moments

% first three are forces along x y z
% next three are moments along x y z 
rv = [0 0 0]'; %Since body frame is established at CV, the position ...
FAS = -rho*g*VB*g_cap; %Aerostatic forces N
tau_AS = [FAS; cross(rv,FAS)]; %Aerostatic forces and moments

%% Control forces and moments
% Constants

 V_bar  = sqrt(v^2 + u^2 + w^2) ; %true airspeed in m/s
 Sf = 65; % Surface area of the fin m^2
 c= 7 ; % Chord length of the fin
 cf = 3; %Flap chord
 
 S_FE = 65; %reference surface area m^2
 %
 delta = deg2rad(10);
 %lift curve slope
 CL_alpha = 5.7; %Change as needed
 
 %theoretical effectiveness factor tau
 theta_f = acos(2*cf/c -1);
 tau = 1- (theta_f - sin(theta_f)) /pi;
 
 %correction factor eta_d
 eta_d = 0.8; %simplified fucntion
 %3D effectiveness value
 k3D = 1.2; %example value
 %calc DeltaCL
 
 DeltaCL = CL_alpha *tau * eta_d * k3D * delta;
 
 DeltaCD = 1.7* (cf/c)^1.38 * (Sf / S_FE) * sin(delta)^2;
 
 %Moment calcs for each fin
 numerator = -2*sin(theta_f) - sin(2*theta_f);
 denominator = 8 * (pi - theta_f + sin(theta_f));
 
 ratio = numerator/denominator;
 
 DeltaCM1_4 = ratio*DeltaCL;
 
 %Deflection anges degrees
 theta = [45, 135, 225, 315];
 
 %initialize force and moment arrays
 
 Fx = zeros(1,4);
 Fy = zeros(1,4);
 Fz = zeros(1,4);
 
 L = zeros(1,4);
 M = zeros(1,4);
 N = zeros(1,4);
 
 %calc forces and moements for each fin
 
 for i = 1:4
     %conv theta to radians
     theta_rad = deg2rad(theta(i));
 
     %calc forces
     Fx(i) = -0.5 * rho * V_bar^2 *Sf*DeltaCD;
     Fy(i) = 0.5*rho*V_bar^2*Sf*(DeltaCL *sin(theta_rad)+ DeltaCD * cos(theta_rad));
     Fz(i) = -0.5 * rho * V_bar^2 *Sf*(DeltaCL *cos(theta_rad)- DeltaCD * sin(theta_rad));
 
     L(i) = Fz(i); %Assume y_i = 1 for simplicity
     M(i) = 0.5*rho*V_bar^2*Sf*c*DeltaCM1_4;
     N(i) = Fy(i); %Assume x_i = 1 for simplicity 
 end
 % Totals
 TotalFx = sum(Fx(:));
 TotalFy = sum(Fy(:));
 TotalFz = sum(Fz(:));
 
 TotalL = sum(L(:));
 TotalM = sum(M(:));
 TotalN = sum(N(:));

 tau_C = [TotalFx TotalFy TotalFz TotalL TotalM TotalN]';

%For natural response 
tau_C = [0 0 0 0 0 0]';

%% Airship thrust force and moments
%Duct fan
rTDF1 = [3.5 7.87 21]'; % Position of motor wrt CV gomes pg 179
rTDF2 = [3.5 -7.87 21]';
fTDF = 4000; %Thrust at sea level; N
etaTDF = 0.6;% 
FTDF1 = etaTDF * [fTDF 0 0 ]'; %Thrust force N
FTDF2 = etaTDF * [fTDF 0 0 ]';
FTDF = FTDF1+FTDF2;
MTDF1 = cross(rTDF1, FTDF1); %Moments N-m
MTDF2 = cross(rTDF2, FTDF2);
MTDF = MTDF1 + MTDF2;

%Turbo prop
rTTP = [0 0 20.1]'; % Position of motor wrt CV
fTTP = 2000;%600;
etaTP = 0.6;
FTTP = etaTP* [fTTP 0 0 ]';
MTTP = cross(rTTP, FTTP);

FT = FTTP + FTDF;
MT = MTTP + MTDF;
tau_Th = [FT; MT];


%% Estimation of added mass matrix
% MAH - added mass of airship hull
B = VB*rho*g; %Bouyancy force
m_dash = B/g;
%I_dash = sum(m_dash * (L.^2 + D^2) / 20); 
I_dash = (m_dash * (L.^2 + D^2) / 20); 

mH11 = k1*m_dash;
mH22 = k2*m_dash;
mH33 = mH22;
mH44 = 0;
mH55 = k_dash*I_dash;
mH66 = mH55;
mH = [mH11, mH22, mH33, mH44, mH55, mH66];

MAH = diag(mH); %6x6 diagonal matrix

%% MAF calcs Do once!
%{
syms x
format short
rho_dt = 1.158;
mF22 = etaf*int(rho_dt*pi*(x - (R.^2/x)).^2 ,[R b]);
mF33 = mF22;
mF35 = -etaf*int(rho_dt*pi*x*(x - (R.^2/x)).^2 ,[R b]);
mF26 = -mF35;
mF44 = etaf*int((2/pi)*k44*rho_dt*x^4 ,[R b]);
mF55 = etaf*int(x^2*rho_dt*pi*(x - (R.^2/x)).^2, [R b]);
mF66 = mF55;

MAF = zeros(6,6);
MAF(2,2) = mF22;
MAF(2,6) = mF26;
MAF(3,3) = mF33;
MAF(3,5) = mF35;
MAF(4,4) = mF44;
MAF(5,5) = mF55;
MAF(6,6) = mF66;
%}



%% Aerodynamic forces and moments
% symmetric added mass axis
% Aerodynamic forces N and moments N-m

%extract components from the syymetric added mass axis
M11 = MA_rigid(1:3, 1:3);
M12 = MA_rigid(1:3, 4:6);
M21 = MA_rigid(4:6, 1:3);
M22 = MA_rigid(4:6, 4:6);

% Calc aero forces and moments

forces_due_to_rotation = cross(omg, (M11 * V + M12 *omg));
moments_due_to_translation = cross(V, (M11*V+M12*omg));
moments_due_to_rotation = cross (omg, (M21*V + M22*omg));

%combine into aero forces and moments (sign change - to +)
tau_A_dash = [forces_due_to_rotation; moments_due_to_translation+moments_due_to_rotation];

%% Viscous forces
FVN = 0; %Viscous force
MVN = 0; %Viscous Moment
FV = FVN*[0 (-v/sqrt(v^2 + w^2)) (-w/sqrt(v^2 +w^2)) ]';
MV = MVN*[0 (w/sqrt(v^2 + w^2)) (-v/sqrt(v^2 +w^2)) ]';
tau_V = [FV; MV];

%tau_V = [0 0 0 0 0 0]'; %if slow speed

%% Return
V_dot = zeros(6,1);
MA_rigid = (MAH+MAF);


%                                       +tau_Imp for gusts
V_dot = (MA_rigid + Mrigid) \ (tau_I + tau_G + tau_C + tau_AS + tau_A_dash + tau_Th + tau_V);


end
