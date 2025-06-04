function XYZ_dot = position_wrt_Earth(xyz)
u = xyz(1);
v = xyz(2);
w = xyz(3);

phi = xyz(7);
theta = xyz(8);
psi = xyz(9);

R_be = [cos(theta)*cos(psi) , sin(theta)*sin(phi)*cos(psi) - cos(theta)*sin(psi), cos(phi)*cos(psi)*sin(theta)+ sin(phi)*sin(psi); ...
    cos(theta)*sin(psi), sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi), cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi); ...
    -sin(theta) , sin(phi)*cos(theta) , cos(phi)*cos(theta)];

body_velocity = [u; v; w];
XYZ_dot = R_be*body_velocity;
end