% ----------------------------------------------------------------
% General Geometry and Properties
% ----------------------------------------------------------------

% Geometry
% Hull Geometry
hull_length = 4.7;
hull_max_diameter = 1.6;

% Gondola Geometry
gondola_length = 0.4;
gondola_width = 0.17;
gondola_height = 0.1;

% Fin Geometry
fin_root_chord = 0.465;
fin_thickness = 0.01;
fin_span = 0.39;

% Center of volume measured from the nose
hull_cov = 2.22;

% Distances
% Gondola Distance from Nose and Height
d_gon_x = 2;
r_gon_y = 0;
r_gon_z = 0.8486;

% Payload Distance from Nose and Height
d_pl_x = 2.4;
r_pl_z = 0.784;

% Thrust Point Distance
d_p_x = 0.1;
r_p_z = 0.8486;

% Fin Position
d_fin_x_s = 3.635;
r_fin_rad = 0.5634;

% Actuator Properties
d_mot_y = 0.45;
mot_rod_r = 0.01;
mot_prop_r = 0.1651;
mot_prop_h = 0.02;
rud_prop_r = 0.0635;
rud_prop_h = 0.02;

% Masses
hull_mass = 1.7;
hull_volume = 5.97;
gondola_mass = 3.3;
pl_mass = 1.05;
fin_mass = 0.055;

% Material Properties
air_density = 1.225;
gas_density = 0.175;

% Motor Properties
mot_max_vel = 1500;
mot_max_thrust_kg = 1.62;
mot_tau = 0.05;
mot_max_angle = 135 * pi / 180;
rud_max_vel = 1500;
rud_max_thrust_kg = 0.24;
rud_tau = 0.025;

% Aerodynamic Properties
k_axial = 0.04662;
k_transverse = 0.9147;
k_rotation = 0.75567;
m_added_fins_k44 = 1;
m_added_fins_efficiency_factor = 0.191395;

% ----------------------------------------------------------------
% Positions
% ----------------------------------------------------------------

origin_height = 2;
origin_x = hull_cov - hull_length / 2;
origin_y = 0;
origin_z = origin_height + hull_max_diameter / 2 + gondola_height + 0.2;

% ----------------------------------------------------------------
% Calculations
% ----------------------------------------------------------------

% Gas Mass
gas_mass = gas_density * hull_volume;

% Total Mass
total_mass = hull_mass + gas_mass + gondola_mass + 4 * fin_mass + pl_mass;

% Distances (Using Hull Ellipsoid Geometry)
r_gon_x = hull_cov - d_gon_x;
r_pl_x = hull_cov - d_pl_x;
r_p_x = hull_cov - (d_p_x + d_gon_x - gondola_length / 2);
r_fin_x = hull_cov - (d_fin_x_s + fin_root_chord / 2);
r_fin_center_rad = r_fin_rad + fin_span / 2;

% Center of Gravity Coordinate
r_g_x = (gondola_mass * r_gon_x + 4 * fin_mass * r_fin_x + pl_mass * r_pl_x) / total_mass;
r_g_y = 0;
r_g_z = (gondola_mass * r_gon_z + pl_mass * r_pl_z) / total_mass;

% Moment of Inertia
hull_ixx = (1 / 3) * hull_mass * ((hull_max_diameter / 2) ^ 2 + (hull_max_diameter / 2) ^ 2);
hull_iyy = (1 / 3) * hull_mass * ((hull_max_diameter / 2) ^ 2 + (hull_length / 2) ^ 2);
hull_izz = (1 / 3) * hull_mass * ((hull_max_diameter / 2) ^ 2 + (hull_length / 2) ^ 2);

gas_ixx = (1 / 5) * gas_mass * ((hull_max_diameter / 2) ^ 2 + (hull_max_diameter / 2) ^ 2);
gas_iyy = (1 / 5) * gas_mass * ((hull_max_diameter / 2) ^ 2 + (hull_length / 2) ^ 2);
gas_izz = (1 / 5) * gas_mass * ((hull_max_diameter / 2) ^ 2 + (hull_length / 2) ^ 2);

gondola_ixx = (1 / 12) * gondola_mass * (gondola_height ^ 2 + gondola_width ^ 2) + gondola_mass * r_gon_x ^ 2;
gondola_iyy = (1 / 12) * gondola_mass * (gondola_height ^ 2 + gondola_length ^ 2) + gondola_mass * r_gon_y ^ 2;
gondola_izz = (1 / 12) * gondola_mass * (gondola_width ^ 2 + gondola_length ^ 2) + gondola_mass * r_gon_z ^ 2;

payload_ixx = pl_mass * r_pl_x ^ 2;
payload_iyy = 0;
payload_izz = pl_mass * r_pl_z ^ 2;

fins_ixx = 4 * ((1 / 12) * fin_mass * (fin_span ^ 2) + fin_mass * r_fin_x ^ 2);
fins_iyy = 2 * ((1 / 12) * fin_mass * fin_root_chord ^ 2 + fin_mass * r_fin_center_rad ^ 2);
fins_izz = fins_iyy;

ixx = hull_ixx + gas_ixx + gondola_ixx + fins_ixx + payload_ixx;
iyy = hull_iyy + gas_iyy + gondola_iyy + fins_iyy + payload_iyy;
izz = hull_izz + gas_izz + gondola_izz + fins_izz + payload_izz;

% Actuator Properties
sim_rotor_slow = 10;
mot_coeff = (mot_max_thrust_kg * 9.81) / (mot_max_vel ^ 2);
rud_coeff = (rud_max_thrust_kg * 9.81) / (rud_max_vel ^ 2);

% ----------------------------------------------------------------
% Printing Outputs
% ----------------------------------------------------------------

fprintf('Gas Mass: %.4f kg\n', gas_mass);
fprintf('Total Mass: %.4f kg\n', total_mass);
fprintf('Center of Gravity: (%.4f, %.4f, %.4f) m\n', r_g_x, r_g_y, r_g_z);
fprintf('Moment of Inertia (Ixx, Iyy, Izz): (%.4f, %.4f, %.4f) kg.m^2\n', ixx, iyy, izz);
fprintf('Motor Coefficient: %.4f\n', mot_coeff);
fprintf('Rudder Coefficient: %.4f\n', rud_coeff);
