function planets_elements = elements_from_ephems(planets_list, jd_time)
% This function computes the actual orbital elements for the requested
% Julian Date for the planets specified in the input
jd_T = juliandate(datetime('2000-01-01 12:00:00', "TimeZone", 'UTC'));
T = (jd_time - jd_T)/36525;
planets_ephem0 = ephems(planets_list);

elements_names = ["a", "e", "i", "raan", "argp", "M"];
planets_elements = struct();

for planet = 1:length(planets_list)
 
    % Retrieve the base elements for this planet
    base_elements = planets_ephem0.(planets_list{planet});
    
    % Update semimajor axis a and other orbital elements
    a = base_elements('a') + base_elements('a_dot')*T;
    e = base_elements('e') + base_elements('e_dot')*T;
    i = base_elements('i') + base_elements('i_dot')*T;
    L = base_elements('L') + base_elements('L_dot')*T;
    long_peri = base_elements('long_peri') + base_elements('long_peri_dot')*T;
    raan = base_elements('long_node') + base_elements('long_node_dot')*T;
    
    % Compute the argument of perihelion and mean anomaly
    argp = long_peri - raan;
    M = L - long_peri;

    % Adjust mean anomaly to be between -180° and 180°
    M = mod(M, 360);
    if M > 180
        M = M - 360;
    end

    % Build the planet structure
    elements_values = [a, e, i, raan, argp, M];
    planets_elements.(planets_list{planet}) = struct(elements_names(1), elements_values(1), elements_names(2), elements_values(2), elements_names(3), elements_values(3), elements_names(4), elements_values(4), elements_names(5), elements_values(5), elements_names(6), elements_values(6));
end