close all
clc

% Astronomical Unit
au = 149597870.7;       % km

% Planet radius
R_mars = 3389.5;        % km  
R_jupiter = 71492;      % km
R_saturn= 58.232;       %km

% Gravitational parameters
mu_sun = 1.32712440018e11; % [km^3/s^2]
mu_earth = 3.986004418e5;  % [km^3/s^2]
mu_mars = 42828.37;        % [km^3/s^2]
mu_jupiter = 1.26686534e8; % [km^3/s^2]

mu_sun_au = mu_sun/au^3;         % [au^3/s^2]
mu_earth_au = mu_earth/au^3;     % [au^3/s^2]

% Ecliptic plane inclination with respect to earth equatorial plane
i_ecl = 23.43928;

%% Satellite initial orbit
% Define orbital parameters
sat.orbit0.a = 7500;           % km
sat.orbit0.e = 0;              
sat.orbit0.i = 0;              % deg
sat.orbit0.raan = 0;           % deg
sat.orbit0.argp = 0;           % deg
sat.orbit0.nu = -240.4281;     % deg

%% PHASE 0
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% IMPORT THE EPHEMERIDES AND BUILD THE SOLAR SYSTEM ENVIRONMENT
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%% Initialization
% Select starting date and convert it in Julian Date
timezone = 'UTC';
start_date = datetime('2005-03-04 12:00:00', "TimeZone", timezone);
mars_fb_date = datetime('2007-02-25 12:00:00', "TimeZone", timezone);
end_date = datetime('2014-10-19 12:00:00', "TimeZone", timezone);

jd_start = juliandate(start_date);
jd_mars_fb = juliandate(mars_fb_date);
jd_end = juliandate(end_date);

% Select planets to visualize
planets = {'venus', 'earth', 'mars', 'jupiter', 'saturn'};

% Compute planets orbital elements at starting epoch time, and consider
% them constant throughout the entire mission (except for the mean anomaly)
planets_elements = elements_from_ephems(planets, jd_start);

% Compute planets states from start to end date with a step of a day or
% some hours
step = 0.25; %6 ore (0.25giorno)
jd_vec = jd_start:step:jd_end;
[day_vec, r_venus, v_venus] = planet_orbit_coplanar(planets_elements.venus, jd_start, jd_end, jd_vec);
t_interp = jd_vec*24*60*60;

% Find closer idx to jd_mars_fb e jd_jupiter_arrival in jd_vec
[~, idx_mars_fb] = min(abs(jd_vec - jd_mars_fb));

% Find position and velocity vectors
[~, r_earth, v_earth] = planet_orbit_coplanar(planets_elements.earth, jd_start, jd_mars_fb, [jd_start, jd_mars_fb]);


[~, r_mars_fb, v_mars_fb] = planet_orbit_coplanar(planets_elements.mars, jd_start, jd_mars_fb, [jd_start, jd_mars_fb]);
r_mars_fb = r_mars_fb(:, end); 
v_mars_fb = v_mars_fb(:, end);

%% PHASE 1
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% ORBIT CHANGE OF PLANE, IN ORDER TO HAVE ORBIT PLANE CONCIDING WITH THE
% ECLIPTICAL PLANE. THEN DO THE ORBITAL COORDINATES CHANGE FROM EQUATORIAL
% TO ECLIPTICAL EARTH-CENTERED FRAME
% -----------------------------------------------------------------------------------------------------------------------------------------------------

% Compute initial position and velocity in ECI equatorial frame
[r0_eci, v0_eci] = oe2rv(sat.orbit0.a, sat.orbit0.e, sat.orbit0.i, sat.orbit0.raan, sat.orbit0.argp, sat.orbit0.nu);

% Compute delta V for plane change and new state after maneuver
% Velocity on the circular parking orbit
v_c = sqrt(mu_earth / sat.orbit0.a);

%change from sat.orbit0.i to i_ecl
delta_inc_rad = deg2rad(i_ecl - sat.orbit0.i);
delta_v_plane_change = 2 * v_c * sin(abs(delta_inc_rad) / 2); 
% Salviamo il dato nella struttura del satellite e aggiorniamo i 
sat.deltaV_plane_change = delta_v_plane_change;
sat.orbit0.i = i_ecl; 

[r_post_manoeuver, v_post_manoeuver] = oe2rv(sat.orbit0.a, sat.orbit0.e, sat.orbit0.i, sat.orbit0.raan, sat.orbit0.argp, sat.orbit0.nu);


% Change coordinates into ecliptic ECI frame
[r_ecl, v_ecl] = eci_eq2ecl(r_post_manoeuver, v_post_manoeuver);

%% PHASE 2
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% COMPUTATION OF THEORETICAL VALUES V_sp AND V_sa FOR EARTH AND MARS. 
% THEORETIC DESIGN OF THE TRAJECTORY WITH GRAVITY ASSIST.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%% Part 1: Escape trajectory from Earth to Mars for a fly-by
% Compute trajectory from Earth to Jupiter: orbital parameters of the
% transfer trajectory and initial ad final velocities
deltaT_earth_mars = (jd_mars_fb - jd_start)*24*60*60;

% Compute initial and final velocities in au/s
[v_earth_sp_appr, v_mars_sa_appr, ~, exitflag] = lambert(r_earth(:, 1)', r_mars_fb', deltaT_earth_mars, 0, mu_sun_au, 'au', 'sec');
v_earth_sp_appr = v_earth_sp_appr'; %velocità eliocentrica di partenza 
v_mars_sa_appr = v_mars_sa_appr'; %velocità eliocentrica di arrivo

% Compute transfer trajectory orbital parameters
sat.orbit_transfer_earth_mars = rv2oe(r_earth(:, 1), v_earth_sp_appr, mu_sun_au);

%% PHASE 3
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% COMPUTE DELTA-V AND TIME OF MANEUVER TO ESCAPE FROM EARTH SOI
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%% Part 1: Escape manoeuver from Earth - Calcolo Orbita di fuga 
% Calculate the required deltaV of the manoeuver so that we
% reach V_sp at Earth SoI, and the exact time of the maneuver
v_inf = v_earth_sp_appr - v_earth(:, 1);
v_c = sqrt(mu_earth/sat.orbit0.a);

% Compute hyperbola parameters
sat.orbit_escape_earth.a = - mu_earth/(norm(v_inf)*au)^2;
sat.orbit_escape_earth.e = 1 + (norm(v_inf)*au)^2/v_c^2;

% Compute phase angle theta_f and alpha, and maneuver anomaly on the
% initial circular orbit
theta_f = pi + acos(1/sat.orbit_escape_earth.e); %angolo asintoto 
alpha = atan2(v_earth(2, 1), v_earth(1, 1));
%alpha = atan2(v_inf(2), v_inf(1));
sat.orbit0.nu_manoeuver = rad2deg(alpha + theta_f - 2*pi); %punto esatto accensione motori 

% Compute hyperbola perigee position vector equal to position at
% nu_maneuver position along circular orbit. Compute velocity direction.
[r_manoeuver_eci, v_manoeuver_circ_eci] = oe2rv(sat.orbit0.a, sat.orbit0.e, sat.orbit0.i, sat.orbit0.raan, sat.orbit0.argp, sat.orbit0.nu_manoeuver);
[r_manoeuver_ecl, v_manoeuver_circ_ecl] = eci_eq2ecl(r_manoeuver_eci, v_manoeuver_circ_eci);
r_esc_ecl_eci = r_manoeuver_ecl;
v_direction = v_manoeuver_circ_ecl / norm(v_manoeuver_circ_ecl);

% Compute datetime of maneuver (kepler_direct from nu to nu_maneuver)
delta_t_wait = kepler_direct(mu_earth, sat.orbit0.a, sat.orbit0.e, sat.orbit0.nu, sat.orbit0.nu_manoeuver);
jd_esc_maneuver = jd_start + delta_t_wait / 60/60/24;

% Compute the velocity vector on the escape trajectory (after maneuver)
v_hyperbola_perigee_mag = sqrt(mu_earth * (2 / norm(r_esc_ecl_eci) - 1 / sat.orbit_escape_earth.a)); 
v_esc = v_hyperbola_perigee_mag * v_direction;

% Delta V richiesto
sat.deltaV_escape_earth = v_hyperbola_perigee_mag - v_c;

% Apply small corrections

% Compute Earth escape deltaV

%% PHASE 4
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% PROPAGATION OF COMPUTED INITIAL VALUES AFTER ESCAPE MANOEUVER FROM EARTH
% -----------------------------------------------------------------------------------------------------------------------------------------------------

%% Propagation from Earth orbit to Earth SoI limit
% SOI-Sphere of Influence radius
soi_earth = 924600;     % km
soi_mars = 577400;      % km
soi_jupiter = 48223000; % km
soi_saturn = 54432000;  % km

% Propagate from Earth since the time of escape manoeuver and find the time
% at which satellite is outside the Earth SoI (Stop Condition)
options = odeset('RelTol', 2.22045e-14, 'AbsTol', 1e-18, 'Events', @(t, y) stopCondition(t, y, soi_earth));
state0_sat_earth_escape = [r_esc_ecl_eci; v_esc]; %punto di partenza simulazione (calcolato prima) 
gg = 30;
t_vec_escape = linspace(jd_esc_maneuver*24*60*60, (jd_esc_maneuver+gg)*24*60*60, gg*24*60/30);
[t_vec_escape, state_sat_earth_escape] = ode45(@(t, y) satellite_ode(t, y, mu_earth), t_vec_escape, state0_sat_earth_escape, options);
r_sat_earth_escape = state_sat_earth_escape(:, 1:3)';
v_sat_earth_escape = state_sat_earth_escape(:, 4:6)';

% Compute jd time and Earth position when satellite is at Earth SoI
jd_earth_sp = t_vec_escape(end)/24/60/60; % momento esatto in cui si ha uscita Earth SoI
earth_soi_date = datetime(jd_earth_sp,'convertfrom','juliandate','Format','d-MMM-y HH:mm:ss', 'TimeZone', timezone);
[~, r_earth_sp, v_earth_sp] = planet_orbit_coplanar(planets_elements.earth, jd_start, jd_earth_sp, [jd_start, jd_earth_sp]);
r_earth_sp = r_earth_sp(:, end); % posizione e velocità della terra al momento in cui sat buca SoI 
v_earth_sp = v_earth_sp(:, end);

% Convert from ECI to J2000 absolute frame
r_sat_earth_sp = r_sat_earth_escape(:, end)/au + r_earth_sp;
v_sat_earth_sp = v_sat_earth_escape(:, end)/au + v_earth_sp;


%% Propagate from outside Earth SoI to Mars SoI
% Propagate from outside Earth SoI till the satellite enters the Mars SoI

% Check if Mars SoI has been reached

% Compute jd time and Mars position when satellite is entering Mars SoI limit


% Convert from J2000 absolute frame to Mars-Centered frame