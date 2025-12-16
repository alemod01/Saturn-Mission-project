function [t_vec_days, r_vec, v_vec] = planet_orbit_coplanar(planet_elements, jd_start, jd_end, jd_vec)
% This function computes the orbit of a planet or satellite in the
% J2000 frame from initial time jd_start to final time jd_end
% expressed as Julian Dates with the simplified hypotesis of zero
% inclination for every planet, so that they all live on the ecliptic
% plane.
%
% Parameters:
% - planet_elements: Structure containing the classical orbital elements 
%                     of the planet or satellite at the reference epoch.
%                     Inclination is assumed to be zero.
% - jd_start   : Initial Julian Date from which the orbital propagation begins.
% - jd_end     : Final Julian Date at which the propagation stops.
% - jd_vec     : Vector of Julian Dates at which the position and velocity 
%                 are to be evaluated.
% 
% Returns: 
% - t_vec_days : Vector of time values expressed in days from jd_start.
% - r_vec      : Matrix of position vectors (3×N) in the J2000 frame.
% - v_vec      : Matrix of velocity vectors (3×N) in the J2000 frame.

au = 149597870.7;               % au
mu_sun = 132712439935/au^3;     % au^3 s^-2

% Compute the mean motion of the planet
n = sqrt(mu_sun/(planet_elements.a^3));

if jd_vec == 0
    t_vec_days = (0 : 1 : jd_end - jd_start);   % days
else
    t_vec_days = jd_vec - jd_vec(1);
end

t_vec = t_vec_days*24*60*60;                % seconds

% Compute the mean anomaly M for each day from initial to final date

M0_rad = deg2rad(planet_elements.M);
M_vec_rad = M0_rad + n * t_vec;
M_vec_deg = rad2deg(M_vec_rad);

M_vec_deg_mod = mod(M_vec_deg, 360);

% Compute the eccentric anomaly solving the inverse kepler problem
E_vec_deg = kepler_inverse(planet_elements.e, M_vec_deg_mod, 0, 1e-8);
E_vec_rad = deg2rad(E_vec_deg);

% Compute true anomaly
nu = 2*atan(sqrt((1+planet_elements.e)/(1-planet_elements.e))*tan(E_vec_rad/2));

x_peri = planet_elements.a * (cos(E_vec_rad) - planet_elements.e);
y_peri = planet_elements.a * sqrt(1 - planet_elements.e^2)*sin(E_vec_rad);
z_peri = zeros(1, length(t_vec));

r_peri = [x_peri; y_peri; z_peri];

vx_peri = sqrt(mu_sun/(planet_elements.a*(1-planet_elements.e^2)))*(-sin(nu));
vy_peri = sqrt(mu_sun/(planet_elements.a*(1-planet_elements.e^2)))*(planet_elements.e*ones(1, length(t_vec)) + cos(nu));
vz_peri = zeros(1, length(t_vec));
v_peri = [vx_peri; vy_peri; vz_peri];

% Rotate along z-axis of raan + argp (Assuming coplanar orbits: inclination
% equal to zero)
alpha = deg2rad(planet_elements.raan + planet_elements.argp);

R = [cos(alpha), - sin(alpha), 0; sin(alpha), cos(alpha), 0; 0, 0, 1];
r_vec = R*r_peri;
v_vec = R*v_peri;
end