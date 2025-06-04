% Load the logged data
time = out.EA.Time; 
euler_angles = squeeze(out.EA.Data)';  % [phi, theta, psi]
positions = squeeze(out.XYZ.Data)';  % [X, Y, Z]
positions(:,3) = -positions(:,3);
pos = positions(1, :);

% Define airship shape (simple ellipsoid)
[x, y, z] = ellipsoid(0, 0, 0, 5, 2, 2, 20); % Scaled for visualization
% Create figure
figure;
hold on;
grid on;
xlabel('X Position');
ylabel('Y Position');
zlabel('Altitude');
title('Airship Animation');
view(3);
axis equal;

% Initialize airship patch
airship = surf(x, y, z, 'FaceColor', 'cyan', 'EdgeColor', 'none');
axis([-2 100 -30 30 -30 30]); % Adjust limits based on your airship's expected range

% Initialize body frame axes
x_axis = quiver3(pos(1), pos(2), pos(3), 1, 0, 0, 'r', 'LineWidth', 2);
y_axis = quiver3(pos(1), pos(2), pos(3), 0, 1, 0, 'g', 'LineWidth', 2);
z_axis = quiver3(pos(1), pos(2), pos(3), 0, 0, 1, 'b', 'LineWidth', 2);

% Animation loop
for i = 1:length(time)/10
    % Extract current position and orientation
    phi = euler_angles(i*10,1);
    theta = euler_angles(i*10,2);
    psi = euler_angles(i*10,3);
    pos = positions(i*10, :);

    
    % Create rotation matrix
    R = eul2rotm([psi, theta, phi]); % Convert Euler angles to rotation matrix

    % Transform airship vertices
    verts = [x(:), y(:), z(:)] * R'; % Rotate
    verts = verts + pos; % Translate

    % Update airship position
    set(airship, 'XData', reshape(verts(:,1), size(x)), ...
                 'YData', reshape(verts(:,2), size(y)), ...
                 'ZData', reshape(verts(:,3), size(z)));
% Update body frame axes
    transformed_axes = R * eye(3) * 10; % Scale for visibility
    set(x_axis, 'XData', pos(1), 'YData', pos(2), 'ZData', pos(3), ...
                'UData', transformed_axes(1,1), 'VData', transformed_axes(2,1), 'WData', transformed_axes(3,1));
    set(y_axis, 'XData', pos(1), 'YData', pos(2), 'ZData', pos(3), ...
                'UData', transformed_axes(1,2), 'VData', transformed_axes(2,2), 'WData', transformed_axes(3,2));
    set(z_axis, 'XData', pos(1), 'YData', pos(2), 'ZData', pos(3), ...
                'UData', transformed_axes(1,3), 'VData', transformed_axes(2,3), 'WData', transformed_axes(3,3));
    % Update title with orientation info
    title(sprintf('Airship Animation - Time: %.2f s | Phi: %.2f | Theta: %.2f | Psi: %.2f', ...
        time(i*10), phi, theta, psi));

    % Pause for animation effect
    pause(0.05);
end

hold off;
