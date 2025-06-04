h=100; %altitude in m

u0 = 0.0; %m/s
v0 = 0.0; %m/s
w0 = 0.0; %m/s

p0 = 0; %rad/s
q0 = 0; %rad/s
r0 = 0; %rad/s

phi0=0; %rad
theta0=0; %rad
psi0=0; %rad

xE0=0; %meters
yE0=0; %meters
zE0=-h; %meters

%% Atmosphere model
rho_0 = 1.225; %Density at sea level in kg/m3
H = 8400; %Scale height in meters;

%% Viscous forces
FVN = 10; %Viscous force
MVN = 0.0010; %Viscous Moment
%%

g = 9.81;
m = 5942; % mass
VB = 4911; %5153; %volume (m^3)
S = VB^(2/3) ; %surface area (m^2)

XB = 0; %x-origin point, m (boyancy)
ZB = 0; % z-origin point, m
XG = 0; %0.0 as given in [2], m
ZG = 5.1816; %m
rGxm = [XG-XB 0 ZG-ZB]';

L = 52; %Hull length
D = 14; %Hull max diameter
R = 5.6522; % Hull cross section radius at fins (m)
b = 14.13; %Hull centre line to fin tip length
delta = -15; %Elevator deflection 5 degrees
%AR = 0.54;% Aspect ratio
AR = L/(2*D);% Aspect ratio
I33 = eye(3);

%%
keta = ketaparameters(L, D, R, b, AR, delta);

 k1 = 0.08;
 k2 = 0.89;
 k_dash = 0.5;
%k1 = keta(1);
%k2 = keta(2);
%k_dash = keta(3);

etaf = 0.8; % Example: Flap effectiveness
k44 = 0.1; % Example: Roll damping coefficient

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

%% control surfaces
 Sf = 65; % Surface area of the fin m^2
 c= 7; % Chord length of the fin
 cf = 3; %Flap chord

 %lift curve slope
 CL_alpha = 5.7; %Change as needed
%%

%Duct fan
rTDF1 = [3.5 7.87 21]'; % Position of motor wrt CV gomes pg 179
rTDF2 = [3.5 -7.87 21]';

fTDF = 0; %4000; %Thrust at sea level; N
etaTDF = 0.6;%

%Turbo prop
rTTP = [0 0 20.1]'; % Position of motor wrt CV
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