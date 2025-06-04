%%
t_s = 0.005;
%%
h=3.7495; %altitude in m
u0 = 4.0; %m/s
u0 = 0.0; %m/s
v0 = 0.0; %m/s
w0 = 0.0; %m/s

p0 = 0; %rad/s
q0 = 0; %rad/s
r0 = 0; %rad/s

phi0=0; %rad roll
theta0=1; %rad pitch
psi0=0; %rad yaw

xE0=0; %meters
yE0=0; %meters
zE0=-h;%meters
%% Atmosphere model
rho_0 = 1.225; %Density at sea level in kg/m3
H = 8400; %Scale height in meters;
rho_he = 0.175 ;%gas density (Helium)kg/m3
%% Viscous forces
%FVN = 5000; %Viscous force
FVN = [5, 5, 5];
MVN = [0.2, 10, 5]; %Viscous Moment
%%

g = 9.81;
%m   = 7.3147; %Cloudship %5942; % mass
m   = 7.31;
VB = 5.97; %Cloudship %5153; %volume (m^3)
S = VB^(2/3) ; %surface area (m^2)

XB = 0;%cloudship %0; %x-origin point, m (boyancy)
ZB = 0; %Cloudship % z-origin point, m
%XG = 0.0239; %Cloudship %0; %0.0 as given in [2], m
XG=0;
ZG = 0.4954; %Cloudship %5.1816; %m
rGxm = [XG-XB 0 ZG-ZB]';

L = 4.7; %Cloudship %52; %Hull length
D = 1.6; %Cloudship %14; %Hull max diameter
R = 0.8; %Cloudship %5.6522; % Hull cross section radius at fins (m)
b = 1.19;%Cloudship %14.13; %Hull centre line to fin tip length
delta_e = 0; %Elevator deflection
delta_r = 0; % rudder deflection
AR = L/(2*D);% Aspect ratio
I33 = eye(3);

%%
%keta = ketaparameters(L, D, R, b, AR, delta);
keta = [0.04662, 0.9147, 0.75567]; % k1, k2, k_dash Cloudship

 % k1 = 0.08;
 % k2 = 0.89;
 % k_dash = 0.5;
 k1 = keta(1);
 k2 = keta(2);
 k_dash = keta(3);

%etaf = 0.8;     % Example: Flap effectiveness
%k44 = 0.1;      % Example: Roll damping coefficient
etaf = 0.191395; %cloudship
k44 = 1; %%cloudship
%eta = keta(1,4);
%etad = keta(1,5);
%etaf = keta(1,6);
%k44 = keta(1,7);
%k3D = keta(1,8);

Ixx = 1.7972; %moment of inerta Nm2
Iyy = 4.9584; %moment of inerta Nm2
Izz = 7.9854; %moment of inerta Nm2

Ixy =0;
Ixz =0;
Iyz =0;
Inertia = [Ixx, Ixy, Ixz; Ixy, Iyy, Iyz; Ixz, Iyz, Izz];
%% control surfaces
 Sf = 0.18135; % Surface area of the fin m^2
 %Sf = 1;
 c= 0.465; %cloudship %7 ; % Chord length of the fin
 cf = 0.1395; %Flap chord

 %lift curve slope
 CL_alpha = 5.7; %Change as needed

%%

%Duct fan
rTDF1 = [0 0.5 0]'; % Position of left motor wrt CV gomes pg 179
rTDF2 = [0 -0.5 0]';% Position of right motor wrt CV gomes pg 179

fProp_max = 20;
fProp1 = 10; %4000; %Thrust at sea level; N
fProp2 = 10;
fProp1 = 0; %4000; %Thrust at sea level; N
fProp2 = 0;

Prop1_alpha = 0; %prop1 swivel angle, 0degrees is aftward, 90degrees is downward
Prop2_alpha = 0;
etaTDF = 1;% 

%Turbo prop
rTTP = [0 0 0]'; % Position of motor wrt CV
fTTP =0 ;%2000;%600;
etaTP = 0.6;
%%

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

