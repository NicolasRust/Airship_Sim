%% Trim and vehicle parameters
m   = 7.31;   % total mass (kg)
rho     = 1.225;    % air density (kg/m^3)
r = 0.8;
S       = pi*r^2;    % reference area for drag (m^2)
u0      = 4;    % trim speed (m/s)
CD_eff  = 1.6/(0.8 * S * u0^2 /2);    % effective drag coefficient at trim Cd = D/(r A V^2 /2)

%% Actuator (thrust) parameters
K_T     = 10;    % steady-state thrust gain (N per unit command)
tau_a   = 0.01;    % actuator time constant (s)

%% Linearize drag around trim
% Steady drag at trim
D0      = 0.5 * rho * S * CD_eff * u0^2;
% Aerodynamic damping derivative in surge
X_u     = rho * S * CD_eff * u0;

%% 4) Compute small-signal thrust trim
T0      = D0;      % at trim, thrust equals drag

%% Build individual transfer functions
% Surge dynamics (thrust to velocity)
num1    = K_T;
den1    = [m, X_u];
G_thrust2vel = tf(num1, den1);

% Actuator dynamics (command to thrust)
num2    = 1;
den2    = [tau_a, 1];
G_cmd2thrust = tf(num2, den2);

%% Combined transfer function (command to velocity)
G_cmd2vel = series(G_cmd2thrust, G_thrust2vel);

%% Display results
disp('Thrust-to-velocity TF (G_thrust2vel):');
G_thrust2vel
%disp('Command-to-thrust TF (G_cmd2thrust):');
%G_cmd2thrust
%disp('Command-to-velocity TF (G_cmd2vel):');
%G_cmd2vel

%% Bode plo
 %bode(G_thrust2vel);
 %grid on;
