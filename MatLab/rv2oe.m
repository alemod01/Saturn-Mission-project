function elements = rv2oe(r, v, mu)
% Computes the orbital elements of a satellite given its position and velocity vectors
% in a central gravitational field with parameter mu.
%
% Input:
%   r    - Position vector [x, y, z] (either in km or au)
%   v    - Velocity vector [vx, vy, vz] (either in km/s or au/s)
%   mu   - Gravitational parameter of the central body (km^3/s^2)
%
% Output:
%   elements - Structure containing the orbital elements:
%              a    - Semi-major axis [km]
%              e    - Eccentricity
%              i    - Inclination [deg]
%              raan - Right Ascension of Ascending Node [deg]
%              argp - Argument of Periapsis [deg]
%              nu   - True Anomaly [deg]

% Magnitudes of position and velocity vectors
r_norm = norm(r);
v_norm = norm(v);

% Specific angular momentum vector (h) and its magnitude
h = cross(r, v);
h_norm = norm(h);

% Inclination (i) in degrees
i = acosd(h(3) / h_norm);

% Node line vector (n), perpendicular to h and z-axis
k = [0, 0, 1];
n = cross(k, h);
n_norm = norm(n);

% Eccentricity vector (e_vec) and eccentricity (e)
e_vec = (1/mu) * ((v_norm^2 - mu/r_norm) * r - dot(r, v) * v);
e = norm(e_vec);

% Semi-major axis (a)
a = 1 / ((2/r_norm) - (v_norm^2 / mu));

% Right Ascension of the Ascending Node (RAAN, Ω) in degrees
if n_norm ~= 0
    raan = acosd(n(1) / n_norm);
    if n(2) < 0
        raan = 360 - raan;  % Adjust quadrant
    end
else
    raan = 0;  % Undefined if the orbit is equatorial
end

% Argument of Periapsis (ω) in degrees
if e > 1e-10
    argp = acosd(dot(n, e_vec) / (n_norm * e));
    if e_vec(3) < 0
        argp = 360 - argp;  % Adjust quadrant
    end
else
    argp = 0;  % Undefined if the orbit is circular
end

% True Anomaly (ν) in degrees
if e > 1e-10
    nu = acosd(dot(e_vec, r) / (e * r_norm));
    if dot(r, v) < 0
        nu = 360 - nu;  % Adjust quadrant
    end
else
    % For circular orbits, true anomaly is angle between n and r
    nu = acosd(dot(n, r) / (n_norm * r_norm));
    if r(3) < 0
        nu = 360 - nu;
    end
end

% Output orbital elements as a structure
elements.a = a;
elements.e = e;
elements.i = i;
elements.raan = raan;
elements.argp = argp;
elements.nu = nu;
end
