% Define plant transfer function
num =  [0.5955 1.938 0.5328];
den = [1 0.9006 0.1659 0.0003742];
G = tf(num, den);

% Design PID controller using pidtune
% You can specify the desired bandwidth, e.g., 1 rad/s
[C, info] = pidtune(G, 'PID');

% Display PID gains
Kp = C.Kp;
Ki = C.Ki;
Kd = C.Kd;

fprintf('PID Gains:\n');
fprintf('Kp = %.4f\n', Kp);
fprintf('Ki = %.4f\n', Ki);
fprintf('Kd = %.4f\n', Kd);

% Plot open-loop and closed-loop step response
T = feedback(C*G, 1);  % Closed-loop system

figure;
step(T);
title('Closed-loop Step Response with PID Controller');

% Bode plot of open-loop system
figure;
margin(C*G);
title('Open-loop Bode Plot with PID Controller');
