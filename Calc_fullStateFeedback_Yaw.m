% Define the transfer function
num = [-0.5335, 2.206];
den = [1, 0.6654, 0.004132];
G = tf(num, den);

% Convert to controllable canonical state-space
[A, B, C, D] = tf2ss(num, den);

% Display system matrices
disp('A matrix:'); disp(A);
disp('B matrix:'); disp(B);
disp('C matrix:'); disp(C);
disp('D matrix:'); disp(D);

% Check controllability
Co = ctrb(A, B);
if rank(Co) < size(A,1)
    error('System is not controllable. Pole placement not possible.');
end

% Choose desired closed-loop poles (example: faster, stable, well-damped)
p_desired = [-0.6 + 0.6j, -0.6 - 0.6j];

% Compute state feedback gain K using pole placement
K = place(A, B, p_desired);

disp('State feedback gain K:');
disp(K);

% Build closed-loop system (A - B*K)
A_cl = A - B*K;

% Simulate closed-loop response to initial condition
x0 = [0.5; 0];  % Example initial yaw angle and rate
t = 0:0.1:20;
[y, t, x] = initial(ss(A_cl, B, C, D), x0, t);

% Plot
figure;
plot(t, y, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Yaw Output');
title('Closed-Loop Response using Pole Placement Controller');
grid on;
