% Define angles in radians
phi = rand();    % roll
theta = rand();  % pitch
psi = rand();    % yaw

% Compute trigonometric values
cphi = cos(phi); sphi = sin(phi);
ctheta = cos(theta); stheta = sin(theta);
cpsi = cos(psi); spsi = sin(psi);

% First representation
DCM1 = [
    cos(theta)*cos(psi), cos(theta)*sin(psi), -sin(theta);
    sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi), sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi), sin(phi)*cos(theta);
    cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi), cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi), cos(phi)*cos(theta)
]';
% R_phi = [1, 0, 0;
%          0, cos(phi), sin(phi);
%          0, -sin(phi), cos(phi)];
% 
% % Rotation matrix for pitch (theta)
% R_theta = [cos(theta), 0, -sin(theta);
%            0, 1, 0;
%           sin(theta), 0, cos(theta)];
% 
% % Rotation matrix for yaw (psi)
% R_psi = [cos(psi), sin(psi), 0;
%          -sin(psi), cos(psi), 0;
%          0, 0, 1];
% 
%  % Overall DCM (yaw -> pitch -> roll)
%  R_be = R_psi * R_theta * R_phi;

% R_be = [
%     cos(theta)*cos(psi), cos(theta)*sin(psi), -sin(theta);
%     sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi), sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi), sin(phi)*cos(theta);
%     cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi), cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi), cos(phi)*cos(theta)
% ];
% % Second representation
% R_phi = [1, 0, 0;
%          0, cos(phi), sin(phi);
%          0, -sin(phi), cos(phi)];
% 
% % Rotation matrix for pitch (theta)
% R_theta = [cos(theta), 0, -sin(theta);
%            0, 1, 0;
%           sin(theta), 0, cos(theta)];
% 
% % Rotation matrix for yaw (psi)
% R_psi = [cos(psi), sin(psi), 0;
%          -sin(psi), cos(psi), 0;
%          0, 0, 1];
% 
%  % Overall DCM (yaw -> pitch -> roll)
%  DCM2 = R_psi * R_theta * R_phi;

DCM2 =[cos(theta)*cos(psi) , sin(theta)*sin(phi)*cos(psi) - cos(theta)*sin(psi), cos(phi)*cos(psi)*sin(theta)+ sin(phi)*sin(psi); ...
    cos(theta)*sin(psi), sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi), cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi); ...
    -sin(theta) , sin(phi)*cos(theta) , cos(phi)*cos(theta)];
% Check for equivalence
difference = norm(DCM1 - DCM2);
disp(['Matrix difference: ', num2str(difference)]);