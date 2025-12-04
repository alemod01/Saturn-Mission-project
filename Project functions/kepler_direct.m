function deltaT = kepler_direct(mu, a, e, E0, E)
% This function solves the direct kepler problem, given the gravitational
% parameter mu, the semi-major axis a, the eccentricity e, and the initial
% and final eccentric anomalies. It returns the smallest positive forward
% time from E0 to E by adding whole revolutions to E if necessary.
%
% Parameters:
% - mu: gravitational parameter [km^3/s^2]
% - a: semi-major axis [km]
% - e: orbital eccentricity (assumed e < 1 for elliptic)
% - E0: initial eccentric anomaly [deg]
% - E: final eccentric anomaly [deg]
% 
% Returns:
% - deltaT: time interval between the two observations [s]

    if e >= 1
        error('kepler_direct: only implemented for elliptical orbits (e < 1).');
    end

    % Convert degrees to radians and normalize to [0, 2*pi)
    E0_rad = mod(deg2rad(E0), 2*pi);
    E_rad  = mod(deg2rad(E),  2*pi);

    % Mean anomalies
    M0 = E0_rad - e .* sin(E0_rad);
    M  = E_rad  - e .* sin(E_rad);

    % Minimal positive difference in mean anomaly: add 2*pi*n to M until > M0
    deltaM = M - M0;
    if deltaM <= 0
        % smallest integer n so that M + 2*pi*n > M0
        n = ceil((M0 - M) / (2*pi));
        M = M + 2*pi * n;
        deltaM = M - M0;
    end

    % Time factor (sqrt(a^3/mu) has units of seconds)
    time_factor = sqrt(a.^3 ./ mu);

    % Final deltaT in seconds
    deltaT = time_factor .* deltaM;

    % numeric safety: enforce non-negative tiny tolerance
    tol = 1e-12;
    if deltaT < tol
        deltaT = 0;
    end
end