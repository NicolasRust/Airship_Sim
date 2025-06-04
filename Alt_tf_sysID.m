% identify_altitude_tf.m
%
% This script identifies a continuous‐time transfer function relating
% elevator deflection (input) to altitude (output) for your airship,
% using data logged in Simulink as “out.delta_e” and “out.Z”.  It assumes
% the simulation was linearized around a forward speed of 4 m/s and a
% neutral‐buoyant altitude of h = 3.7495 m.
%
% Requirements:
%  • Your Simulink run must have produced a timeseries struct “out” with
%    fields:
%       out.delta_e_ref.Time  – time vector (s)
%       out.delta_e_ref.Data  – elevator angle (rad or degrees)
%       out.Z.Time        – time vector (s), same as delta_e_ref.Time
%       out.Z.Data        – altitude (meters)
%  • The System Identification Toolbox (tfest, iddata, compare, etc.)
%    must be installed and on your MATLAB path.
%
% Usage:
%   1. Run your Simulink model (sim) so that “out” is in the workspace.
%   2. In the MATLAB Command Window, type:
%        >> identify_altitude_tf
%   3. The script will:
%       – Build an iddata object from (u, y, Ts)
%       – Detrend the data (remove any constant bias)
%       – Estimate a 2-pole, 1-zero transfer function using tfest
%       – Plot model‐versus‐data comparison
%
% Notes:
%  • You can change np (number of poles) and nz (number of zeros) below
%    if you think the dynamics require a different order.
%  • If your data still has a strong bias (e.g. small drifts in altitude),
%    consider detrending or high-pass filtering before identification.



%% 2. Extract input (u), output (y), and sampling time Ts
% Assuming the same Time vector for delta_e_ref and Z:

    t = out.delta_e_ref.Time;

% Input: elevator deflection (e.g., radians or degrees)
u = out.delta_e_ref.Data;
% Output: altitude in meters
y = out.Z.Data;

% Compute sampling time (assumes uniform sampling)
Ts = t(2) - t(1);
fprintf('Detected sampling time Ts = %.4f s\n', Ts);

%% 3. Form an iddata object
data = iddata(y, u, Ts);
data.Name = 'Airship Altitude Data';
data.InputName  = 'delta_e_ref';
data.InputUnit  = '(rad or deg)';
data.OutputName = 'Altitude';
data.OutputUnit = 'm';

%% 4. (Optional) Visualize raw data
figure;
subplot(2,1,1);
plot(t, u, 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('Elevator Angle');
title('Elevator Command vs Time');

subplot(2,1,2);
plot(t, y, 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('Altitude (m)');
title('Altitude Response vs Time');

%% 5. Detrend (remove constant bias) if needed
% If the altitude data has an offset or drift, detrending helps identification.
data_d = detrend(data);

%% 6. Choose model order
% Here we pick 2 poles, 1 zero (i.e. G(s) = (b1*s + b0)/(s^2 + a1*s + a0)).
np = 2;  % number of poles
nz = 1;  % number of zeros

fprintf('Estimating a transfer function with %d poles and %d zero...\n', np, nz);

%% 7. Estimate the continuous‐time transfer function
% You can add options to tfest if you want weighting, iterations, etc.
sys = tfest(data_d, np, nz);

%% 8. Display the identified model
disp('=== Identified Transfer Function G(s) ===');
disp(sys);

%% 9. Compare model output to measured data
figure;
compare(data, sys);
title('Measured Altitude vs. Model Prediction');

%% 10. Save results to .mat (optional)
save('identified_altitude_model.mat', 'sys', 'data', 'Ts');

% End of script
