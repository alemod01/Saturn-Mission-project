function planet_ephem0 = ephems(planet_list)
% Planets ephemirides elements and their rates taken from
% https://ssd.jpl.nasa.gov/planets/approx_pos.html
% with respect to J2000 and valid for time interval 1800 AD - 2050 AD

% Orbit elements
elements = ["a", "a_dot", "e", "e_dot", "i", "i_dot", "L", "L_dot", "long_peri", "long_peri_dot", "long_node", "long_node_dot"];

% Mercury
merc_values = [0.38709927, 0.00000037, 0.20563593, 0.00001906, 7.00497902, -0.00594749, 252.25032350, 149472.67411175, 77.45779628, 0.16047689, 48.33076593, -0.12534081];
mercury0 = dictionary(elements, merc_values);

% Venus
ven_values = [0.72333566, 0.00000390, 0.00677672, -0.00004107, 3.39467605, -0.00078890, 181.97909950, 58517.81538729, 131.60246718, 0.00268329, 76.67984255, -0.27769418];
venus0 = dictionary(elements, ven_values);

% Earth
earth_values = [1.00000261, 0.00000562, 0.01671123, -0.00004392, -0.00001531, -0.01294668, 100.46457166, 35999.37244981, 102.93768193, 0.32327364, 0.0, 0.0];
earth0 = dictionary(elements, earth_values);

% Mars
mar_values = [1.52371034, 0.00001847, 0.09339410, 0.00007882, 1.84969142, -0.00813131, -4.55343205, 19140.30268499, -23.94362959, 0.44441088, 49.55953891, -0.29257343];
mars0 = dictionary(elements, mar_values);

% Jupiter
jup_values = [5.20288700, -0.00011607, 0.04838624, -0.00013253, 1.30439695, -0.00183714, 34.39644051, 3034.74612775, 14.72847983, 0.21252668, 100.47390909, 0.20469106];
jupiter0 = dictionary(elements, jup_values);


% Create a dictionary with all these dictionaries
planets_ephem0_dict = dictionary('mercury', mercury0, 'venus', venus0, 'earth', earth0, 'mars', mars0, 'jupiter', jupiter0);

% Initialize output as an empty cell array
planet_ephem0 = {};
    
% Loop over input planet_list and retrieve corresponding ephemeris data
for planet = 1:length(planet_list)
    planet_name = lower(planet_list{planet});  % Convert to lowercase for matching
        
    % Check if the planet exists in the 'planets' struct
    if isKey(planets_ephem0_dict, planet_name)
        % Add the corresponding ephemeris data to the output
        planet_ephem0.(planet_name) = planets_ephem0_dict(planet_name);
    else
        error(['Planet ', planet_name, ' not recognized.']);
    end
end
end                              