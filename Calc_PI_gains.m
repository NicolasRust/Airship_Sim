% calculate_pi_gains.m
% This script computes PI gains for a first-order system to achieve fast, no-overshoot response

clc; clear;

%% Plant definition
K = 0.13345;       % Plant gain
a = 0.667;         % Denominator coefficient (s + a)
ts = 1; %settling time
% Transfer function: G(s) = K / (s + a)

%% Desired closed-loop specifications
zeta = 0.707;          % Critical damping (no overshoot)
omega_n = 4/ts;       % Natural frequency for fast response

% Desired characteristic polynomial: s^2 + 2*zeta*omega_n*s + omega_n^2
desired_char_poly = [1, 2*zeta*omega_n, omega_n^2];

% Coefficients from PI + plant (denominator of closed-loop system):
% Denominator: s^2 + (a + K*Kp)*s + K*Ki
% Match this with desired_char_poly: [1, 2*zeta*omega_n, omega_n^2]

%% Solve for Kp and Ki
% Let the PI controller be: C(s) = Kp + Ki/s

% Matching:
% Coeff of s: a + K*Kp = 2*zeta*omega_n
% Constant:   K*Ki = omega_n^2

Kp = (2*zeta*omega_n - a)/K;
Ki = omega_n^2 / K;

%% Display results
fprintf('PI Controller Gains: (might need to swap Kp and Ki)\n');
fprintf('Kp = %.4f\n', Kp);
fprintf('Ki = %.4f\n', Ki);

%% Optional: Verify step response
s = tf('s');
G = K / (s + a);
C = Kp + Ki/s;

T = feedback(C*G, 1);
figure;
step(T);
title('Closed-Loop Step Response');
grid on;