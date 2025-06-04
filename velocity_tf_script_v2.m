%% Trim and vehicle parameters
m   = 7.31;   % total mass (kg)
rho = 1.225;  % air density (kg/m^3)
r   = 0.8;
S   = pi * r^2;  % frontal area (m^2)
u0  = 4;         % trim speed (m/s)

% More accurate drag coefficient
D_trim = 1.6;    % drag force at trim (N)
CD_eff = D_trim / (0.5 * rho * S * u0^2);  % Cd = D / (0.5*rho*S*u^2)

%% Actuator (thrust) parameters
K_T     = 10;      % thrust gain (N per command unit)
tau_a   = 0.01;    % actuator time constant (s)

%% Linearized model
D0  = 0.5 * rho * S * CD_eff * u0^2;  % should equal D_trim
X_u = rho * S * CD_eff * u0;

% Steady-state thrust
T0 = D0;

%% Transfer Functions
G_thrust2vel = tf(K_T, [m, X_u]);
G_cmd2thrust = tf(1, [tau_a, 1]);

% Optional: add aerodynamic lag
tau_aero = 0.2;  % guess; tune this value
G_aero_lag = tf(1, [tau_aero, 1]);

% Chain together
G_total = series(G_cmd2thrust, G_aero_lag);
G_cmd2vel = series(G_total, G_thrust2vel);

%% Plot
figure;
bode(G_cmd2vel);
grid on;
title('Command to Velocity Transfer Function with Aero Lag');

%% Show transfer functions
disp('Drag coefficient:');
disp(CD_eff);
disp('Linearized surge damping X_u:');
disp(X_u);
disp('Thrust-to-velocity TF:');
G_thrust2vel;
% disp('Command-to-velocity TF with lag:');
% G_cmd2vel;
