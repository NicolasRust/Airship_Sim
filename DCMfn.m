function DCM = DCMfn(dc)
phi = dc(1);
theta = dc(2);
psi = dc(3);

% Compute the DCM using the standard rotation matrices
% Rotation matrix for roll (phi)
R_phi = [1, 0, 0;
         0, cos(phi), sin(phi);
         0, -sin(phi), cos(phi)];

% Rotation matrix for pitch (theta)
R_theta = [cos(theta), 0, -sin(theta);
           0, 1, 0;
          sin(theta), 0, cos(theta)];

% Rotation matrix for yaw (psi)
R_psi = [cos(psi), sin(psi), 0;
         -sin(psi), cos(psi), 0;
         0, 0, 1];

% Overall DCM (yaw -> pitch -> roll)
DCM = R_psi * R_theta * R_phi;

end