function EA_dot = Euler_angle_estimation(Eu)
phi = Eu(1);
theta = Eu(2);
psi = Eu(3);

p = Eu(7);
q = Eu(8);
r = Eu(9);

E_dot = [1, sin(phi)*tan(theta) , cos(phi)*tan(theta);
         0, cos(phi) , -sin(phi);
         0, sin(phi)/cos(theta), cos(phi)/cos(theta)];

body_rates = [p; q; r];

EA_dot = E_dot*body_rates;

end