% calculate_pd_gains.m
% General-purpose PD controller design via pole placement

clc; clear;

%% User-defined plant (enter numerator and denominator coefficients)
% Example: G(s) = K / (s + a)
% numG = [0.13345];
% denG = [1, 0.667];
% numG = [-0.5335, 2.206];
% denG = [1, 0.6654, 0.004132];
numG =  [0.5955 1.938 0.5328];
denG = [1 0.9006 0.1659 0.0003742];
numG =[-0.986 1.461];
denG=[1 1.237 0.3881];
numG =[-0.3327 2.514];
denG =[1 0.7054 0.00177];
numG=[-2.053e-05 0.01595];
denG=[1 0.356 0.0007681];
ts = 20; %settling time

G = tf(numG, denG);
disp('Open-loop transfer function G(s):');
G

%% Desired closed-loop specifications
zeta = 0.707;        % Damping ratio (1 = critical damping)
omega_n = 4/ts;     % Natural frequency [rad/s] — affects response speed

% Desired second-order characteristic equation
desired_char_poly = [1, 2*zeta*omega_n, omega_n^2];
disp('Desired characteristic polynomial:');
disp(desired_char_poly);

%% PD controller form: C(s) = Kd*s + Kp

% Let C(s)*G(s) = (Kd*s + Kp)*numG / (s*denG)
% Characteristic equation: denG(s)*s + numG*(Kd*s + Kp) = 0

% Multiply out the terms:
% (s*denG) = [1 0.667] → s*(s + 0.667) = s^2 + 0.667s
% (numG * controller) = 0.13345*(Kd*s + Kp) = [0.13345*Kd, 0.13345*Kp]

% So the characteristic equation becomes:
% s^2 + (0.667 + 0.13345*Kd)s + 0.13345*Kp = desired_char_poly

a1 = denG(2); % 0.667
K = numG(1);  % 0.13345

% Solve equations:
% 0.667 + K*Kd = 2*zeta*omega_n
% K*Kp = omega_n^2

Kd = (2*zeta*omega_n - a1) / K;
Kp = omega_n^2 / K;

%% Display results
fprintf('\nPD Controller Gains:\n');
fprintf('Kp = %.4f\n', Kp);
fprintf('Kd = %.4f\n', Kd);

%% Create PD controller and plot step response
s = tf('s');
C = Kd*s + Kp;

T = feedback(C*G, 1);

figure;
step(T);
title('Closed-Loop Step Response with PD Controller');
grid on;
