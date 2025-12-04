function [r_ecl, v_ecl] = eci_eq2ecl(r_eq, v_eq)
    % eci_eq2ecl Converts ECI Equatorial position and velocity to
    % ECI Ecliptic frame.
    %
    % Parameters:
    %   r_eq   - Satellite position in ECI [km]
    %   v_eq   - Satellite velocity in ECI [km/s]
    %
    % Returns:
    %   r_ecl - Satellite position in ECI Ecliptic [km]
    %   v_ecl - Satellite velocity in ECI Ecliptic [km/s]

    % Define obliquity of the ecliptic in radians
    epsilon = deg2rad(23.43928);

    % Rotation matrix from Equatorial to Ecliptic ECI 
    R_eq2ecl = [1, 0, 0;
                0, cos(epsilon), sin(epsilon);
                0, -sin(epsilon), cos(epsilon)];
            
    % Rotate Equatorial to ecliptic ECI
    r_ecl = R_eq2ecl * r_eq;
    v_ecl = R_eq2ecl * v_eq;
    
end