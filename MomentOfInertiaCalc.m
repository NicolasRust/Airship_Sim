% Cloudship Moment of Inertia Calculation

% Geometry and Masses
hull_length = 4.7;
hull_max_diameter = 1.6;
hull_mass = 1.7;
hull_volume = 5.97;

gondola_length = 0.4;
gondola_width = 0.17;
gondola_height = 0.1;
gondola_mass = 3.3;

fin_root_chord = 0.465;
fin_span = 0.39;
fin_mass = 0.055;

pl_mass = 1.05;

gas_density = 0.175;
gas_mass = gas_density * hull_volume;

% Center locations (relative to hull center of volume)
hull_cov = 2.22; % center of volume from nose

% Gondola
d_gon_x = 2.0; % from nose
r_gon_x = hull_cov - d_gon_x;
r_gon_y = 0;
r_gon_z = 0.8486;

% Payload
d_pl_x = 2.4; % from nose
r_pl_x = hull_cov - d_pl_x;
r_pl_z = 0.784;

% Fins
d_fin_x_s = 3.635; % start of fin from nose
r_fin_x = hull_cov - (d_fin_x_s + fin_root_chord/2);
r_fin_rad = 0.5634;
r_fin_center_rad = r_fin_rad + fin_span/2;

% --- Hull Moments ---
hull_ixx = (1/3) * hull_mass * ((hull_max_diameter/2)^2 + (hull_max_diameter/2)^2);
hull_iyy = (1/3) * hull_mass * ((hull_max_diameter/2)^2 + (hull_length/2)^2);
hull_izz = hull_iyy;

% --- Gas Moments ---
gas_ixx = (1/5) * gas_mass * ((hull_max_diameter/2)^2 + (hull_max_diameter/2)^2);
gas_iyy = (1/5) * gas_mass * ((hull_max_diameter/2)^2 + (hull_length/2)^2);
gas_izz = gas_iyy;

% --- Gondola Moments ---
gondola_ixx = (1/12) * gondola_mass * (gondola_height^2 + gondola_width^2) + gondola_mass * r_gon_x^2;
gondola_iyy = (1/12) * gondola_mass * (gondola_height^2 + gondola_length^2) + gondola_mass * r_gon_y^2;
gondola_izz = (1/12) * gondola_mass * (gondola_width^2 + gondola_length^2) + gondola_mass * r_gon_z^2;

% --- Payload Moments ---
payload_ixx = pl_mass * r_pl_x^2;
payload_iyy = 0;
payload_izz = pl_mass * r_pl_z^2;

% --- Fins Moments ---
fins_ixx = 4 * ((1/12) * fin_mass * (fin_span^2) + fin_mass * r_fin_x^2);
fins_iyy = 2 * ((1/12) * fin_mass * fin_root_chord^2 + fin_mass * r_fin_center_rad^2 + ...
    (1/12) * fin_mass * (fin_root_chord^2 + fin_span^2) + fin_mass * (r_fin_center_rad^2 + 0 * r_fin_x^2));
fins_izz = fins_iyy;

% --- Total Moments of Inertia ---
ixx = hull_ixx + gas_ixx + gondola_ixx + fins_ixx + payload_ixx;
iyy = hull_iyy + gas_iyy + gondola_iyy + fins_iyy + payload_iyy;
izz = hull_izz + gas_izz + gondola_izz + fins_izz + payload_izz;

% Display results
fprintf('Cloudship Moments of Inertia:\n');
fprintf('ixx = %.4f kg·m²\n', ixx);
fprintf('iyy = %.4f kg·m²\n', iyy);
fprintf('izz = %.4f kg·m²\n', izz);
