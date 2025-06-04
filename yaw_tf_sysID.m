% Extract raw Simulink data (likely as timeseries or struct)

% If out.delta_e_ref and out.Z are timeseries:
if isa(out.delta_e_ref, 'timeseries')
    u = out.delta_e_ref.Data;
    t = out.delta_e_ref.Time;
else
    u = out.delta_e_ref;  % If already a numeric vector
    t = out.tout;         % Time vector
end

if isa(out.Z, 'timeseries')
    y = out.Z.Data;
else
    y = out.Z;
end

% Ensure column vectors
u = u(:);
y = y(:);
t = t(:);

% Confirm data lengths match
if length(u) ~= length(y)
    error('Input and output vectors must be the same length.');
end

% Sampling time (assuming uniform)
Ts = mean(diff(t));

% Create iddata object
data = iddata(y, u, Ts);

% Optional preprocessing
data = detrend(data);

% Estimate transfer function
np = 2;  % Number of poles
nz = 1;  % Number of zeros
sys_tf = tfest(data, np, nz);

% Display result
disp('Estimated Transfer Function:');
sys_tf

% Compare output
figure;
compare(data, sys_tf);
title('Model vs Measured Yaw Output');
