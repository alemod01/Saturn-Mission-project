close all
clear all
clc

% Astronomical Unit
au = 149597870.7;       % km

% Planet radius
R_earth = 6371;         % km
R_mars = 3389.5;        % km  
R_jupiter = 71492;      % km
R_saturn= 58232;        % km

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
start_date = datetime('2003-01-02 12:00:00', "TimeZone", timezone);%rosetta parte il 2004-03-02 12:00:00
earth_fb_date = datetime('2005-03-04 12:00:00', "TimeZone", timezone);
mars_fb_date = datetime('2007-02-25 12:00:00', "TimeZone", timezone);
earth_fb1_date = datetime('2007-11-13 12:00:00', "TimeZone", timezone);
earth_fb2_date = datetime('2009-11-13 12:00:00', "TimeZone", timezone);
end_date = datetime('2014-10-19 12:00:00', "TimeZone", timezone);

jd_start = juliandate(start_date);
jd_earth_fb = juliandate(earth_fb_date);
jd_mars_fb = juliandate(mars_fb_date);
jd_earth_fb1 = juliandate(earth_fb1_date);
jd_earth_fb2 = juliandate(earth_fb2_date);
jd_end = juliandate(end_date);

% Select planets to visualize
planets = {'venus', 'earth', 'mars', 'jupiter', 'saturn'};

% Compute planets orbital elements at starting epoch time, and consider
% them constant throughout the entire mission (except for the mean anomaly)
planets_elements = elements_from_ephems(planets, jd_start);


% Compute planets states from start to end date with a step of a day or
% some hours
% step = 1; %step  giorno 
% jd_fine=jd_start+100;
% plot_trajectory(planets_elements.mars, jd_start, jd_fine , step)

% % Find closer idx to jd_mars_fb e jd_jupiter_arrival in jd_vec
% [~, idx_mars_fb] = min(abs(jd_vec - jd_mars_fb));

% Find position and velocity vectors

[~, r_earth, v_earth] = planet_orbit_coplanar(planets_elements.earth, jd_start, jd_earth_fb, [jd_start, jd_earth_fb]);
r_earth_fb_0 = r_earth(:, end); 
v_earth_fb_0 = v_earth(:, end);

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
% COMPUTATION OF THEORETICAL VALUES V_sp AND V_sa FOR EARTH AND EARTH. 
% THEORETIC DESIGN OF THE TRAJECTORY WITH GRAVITY ASSIST.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%% Part 1: Escape trajectory from Earth to Earth for a fly-by
% Compute trajectory from Earth to Jupiter: orbital parameters of the
% transfer trajectory and initial ad final velocities
deltaT_earth_earth = (jd_earth_fb - jd_start)*24*60*60;

% Compute initial and final velocities in au/s
[v_earth_sp_appr, v_earth_sa_appr, ~, exitflag] = lambert(r_earth(:, 1)', r_earth_fb_0', deltaT_earth_earth, 1, mu_sun_au, 'au', 'sec');

if dot(v_earth_sp_appr, v_earth(:,1)) < 0
    fprintf('Traiettoria retrograda rispetto alla Terra → non fisica');
end

v_earth_sp_appr = v_earth_sp_appr'; %velocità eliocentrica di partenza 
v_earth_sa_appr = v_earth_sa_appr'; %velocità eliocentrica di arrivo

% Compute transfer trajectory orbital parameters
sat.orbit_transfer_earth_earth_lambert = rv2oe(r_earth(:, 1), v_earth_sp_appr, mu_sun_au);

%% PHASE 3
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% COMPUTE DELTA-V AND TIME OF MANEUVER TO ESCAPE FROM EARTH SOI
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%% Part 1: Escape manoeuver from Earth - Calcolo Orbita di fuga 
% Calculate the required deltaV of the manoeuver so that we
% reach V_sp at Earth SoI, and the exact time of the maneuver
v_inf = v_earth_sp_appr - v_earth(:, 1) ;   % [AU/s]
v_inf_kms = v_inf * au;                   % [km/s]

v_c = sqrt(mu_earth / sat.orbit0.a);      % [km/s]       

% Compute hyperbola parameters
sat.orbit_escape_earth.a = - mu_earth / norm(v_inf_kms)^2;
sat.orbit_escape_earth.e = 1 + (sat.orbit0.a * norm(v_inf_kms)^2) / mu_earth;

% Compute phase angle theta_f and alpha
theta_f = pi + acos(1 / sat.orbit_escape_earth.e);   % [rad]
alpha = atan2(v_inf_kms(2), v_inf_kms(1));           % [rad]
%alpha = atan2(v_earth(2, 1), v_earth(1, 1));

sat.orbit0.nu_manoeuver = rad2deg(alpha + theta_f - 2*pi);
sat.orbit0.nu_manoeuver = mod(sat.orbit0.nu_manoeuver, 360);

% Position at maneuver point
[r_manoeuver_eci, v_manoeuver_circ_eci] = oe2rv( sat.orbit0.a, sat.orbit0.e, sat.orbit0.i, sat.orbit0.raan, sat.orbit0.argp, sat.orbit0.nu_manoeuver);

[r_manoeuver_ecl, ~] = eci_eq2ecl(r_manoeuver_eci, v_manoeuver_circ_eci);
r_esc_ecl_eci = r_manoeuver_ecl;

% Tangential velocity direction 
nu = deg2rad(sat.orbit0.nu_manoeuver);
v_direction = [-sin(nu); cos(nu); 0];
v_direction = v_direction / norm(v_direction);

% Time of maneuver
delta_t_wait = kepler_direct( mu_earth, sat.orbit0.a, sat.orbit0.e, sat.orbit0.nu, sat.orbit0.nu_manoeuver );
jd_esc_maneuver = jd_start + delta_t_wait / 86400;

% Velocity at hyperbolic perigee
v_hyperbola_perigee_mag = sqrt(mu_earth * (2 / norm(r_esc_ecl_eci) - 1 / sat.orbit_escape_earth.a));
v_esc = v_hyperbola_perigee_mag * v_direction;

% Circular orbit velocity vector
v_circ = v_c * v_direction;

% Delta V richiesto
sat.deltaV_escape_earth = norm(v_esc - v_circ);

%% DEBUG – Escape diagnostics

fprintf('\n================ EARTH ESCAPE DIAGNOSTICS ================\n');

% --- Velocity magnitudes ---
fprintf('||v_inf||              = %.3f km/s\n', norm(v_inf_kms));
fprintf('||v_c (parking)||      = %.3f km/s\n', v_c);
fprintf('||v_perigee_hyper||    = %.3f km/s\n', v_hyperbola_perigee_mag);
fprintf('DeltaV escape          = %.3f km/s\n', sat.deltaV_escape_earth);

% --- Direction checks ---
v_inf_dir   = v_inf_kms / norm(v_inf_kms);
v_esc_dir   = v_esc / norm(v_esc);
v_earth_dir = (v_earth(:,1)*au) / norm(v_earth(:,1)*au);

% Dot products (cosines of angles)
cos_esc_earth = dot(v_esc_dir, v_earth_dir);
cos_inf_earth = dot(v_inf_dir, v_earth_dir);
cos_esc_inf   = dot(v_esc_dir, v_inf_dir);

fprintf('\nDirection cosines:\n');
fprintf('cos(v_esc , v_earth)   = %.4f\n', cos_esc_earth);
fprintf('cos(v_inf , v_earth)   = %.4f\n', cos_inf_earth);
fprintf('cos(v_esc , v_inf)     = %.4f\n', cos_esc_inf);

% --- Angles in degrees (more intuitive)
ang_esc_earth = acosd(cos_esc_earth);
ang_inf_earth = acosd(cos_inf_earth);
ang_esc_inf   = acosd(cos_esc_inf);

fprintf('\nAngles [deg]:\n');
fprintf('angle(v_esc , v_earth) = %.2f deg\n', ang_esc_earth);
fprintf('angle(v_inf , v_earth) = %.2f deg\n', ang_inf_earth);
fprintf('angle(v_esc , v_inf)   = %.2f deg\n', ang_esc_inf);

% --- Interpretation ---
fprintf('\nInterpretation:\n');

if cos_inf_earth > 0
    fprintf('- v_inf is roughly ALIGNED with Earth velocity (leading escape)\n');
else
    fprintf('- v_inf is roughly OPPOSITE to Earth velocity (trailing escape)\n');
end

if abs(ang_esc_inf) < 5
    fprintf('- v_esc direction is consistent with asymptotic v_inf\n');
else
    fprintf('- WARNING: v_esc and v_inf are NOT well aligned\n');
end

fprintf('===========================================================\n');


% -------------------------------------------------------------------------
% CORREZIONE MANUALE DEL TIRO (TARGETING)
% -------------------------------------------------------------------------
% Il problema: Partiamo dalla SOI, non dal centro della Terra.
% La soluzione: Correggiamo leggermente la velocità e l'angolo di uscita.

% 1. Definiamo i fattori di correzione

% Correzioni fly by marte con partenza 2005 e lambert che completa un'orbita
% k_vel = 0.997929;        % Moltiplicatore di velocità (es. 0.999 o 1.001)
k_vel = 0.9979289;       % Moltiplicatore di velocità (es. 0.999 o 1.001)
delta_angle = 0;   % Correzione angolo in gradi (es. +0.5 o -0.5)


% 2. Applichiamo la correzione alla Magnitudine
v_esc_mag_corr = norm(v_esc) * k_vel;

% 3. Applichiamo la correzione alla Direzione
% Ruotiamo il vettore v_direction di 'delta_angle' gradi nel piano dell'Eclittica
theta_corr = deg2rad(delta_angle);
R_corr = [cos(theta_corr), -sin(theta_corr), 0;
          sin(theta_corr),  cos(theta_corr), 0;
          0,                0,               1];
      
v_dir_corr = (R_corr * v_direction)'; % Ruotiamo

% 4. Ricalcoliamo il vettore di fuga finale CORRETTO
v_esc = (v_esc_mag_corr * v_dir_corr)';

% (Stampa di debug per vedere cosa stiamo facendo)
fprintf('\n--- TARGETING CORRECTION APPLIED ---\n');
fprintf('Velocità scalata di: %.4f\n', k_vel);
fprintf('Angolo ruotato di:   %.2f deg\n', delta_angle);

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
options = odeset('RelTol', 2.22045e-14, 'AbsTol', 1e-18, 'Events', @(t, y) stopCondition(t, y, soi_earth, 'exit'));
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

% Compute transfer trajectory orbital parameters
sat.orbit_transfer_earth_earth = rv2oe(r_sat_earth_sp, v_sat_earth_sp, mu_sun_au);

%% Propagate from outside Earth SoI to Earth SoI - interplanetary 
% Propagate from outside Earth SoI till the satellite enters the Earth SoI

state0_interplanetary_earth_earth = [r_sat_earth_sp; v_sat_earth_sp];

% parametro correttivo gg interplanetary
kg = 0;
% calcolo il tempo per arrivare da fuori SoI Terra a dentro SoI terra
t_cruise_total_earth_earth = jd_earth_fb*24*60*60 - t_vec_escape(end); %durata in secondi del viaggio 

t_cruise_total_earth_earth=t_cruise_total_earth_earth+kg*24*60*60;
t_vec_cruise_earth_earth= linspace(t_vec_escape(end), jd_earth_fb*24*60*60  + (kg*24*60*60) ,t_cruise_total_earth_earth/3600);
% propagazione
options_cruise_earth_earth = odeset('RelTol', 2.22045e-14, 'AbsTol', 1e-18, 'Events', @(t, y) stopCondition_interplanetary(t, y, jd_start , planets_elements.earth , soi_earth));
[t_vec_cruise_earth_earth, state_cruise_earth_earth] = ode45(@(t, y) satellite_ode(t, y, mu_sun_au), t_vec_cruise_earth_earth, state0_interplanetary_earth_earth, options_cruise_earth_earth);

% Estrazione Risultati finali (Stato Eliocentrico all'arrivo)
r_sat_interplanetary_earth_earth = state_cruise_earth_earth(end, 1:3)'; % Posizione rispetto al Sole (AU)
v_sat_interplanetary_earth_earth = state_cruise_earth_earth(end, 4:6)'; % Velocità rispetto al Sole (AU/s)

% Calcolo Data Esatta di Arrivo a SoI terra (Julian Date)
jd_earth_arrival_actual = t_vec_cruise_earth_earth(end)/24/60/60;

% 7. Dove si trova la Terra in quel momento preciso?
% Calcoliamo posizione e velocità della terra usando le effemeridi
[~, r_earth_arr, v_earth_arr] = planet_orbit_coplanar(planets_elements.earth, jd_start, jd_earth_arrival_actual, [jd_start, jd_earth_arrival_actual]);
r_earth_arr = r_earth_arr(:, end); % Posizione Terra (AU)
v_earth_arr = v_earth_arr(:, end); % Velocità terra (AU/s)


% Convert from J2000 absolute frame to Earth-Centered frame
% Calcoliamo il vettore relativo (Sonda rispetto alla Terra)
r_rel_au = r_sat_interplanetary_earth_earth - r_earth_arr; % Vettore distanza in AU
v_rel_au = v_sat_interplanetary_earth_earth - v_earth_arr; % Vettore velocità relativa in AU/s

% Convertiamo in km e km/s per il controllo finale
r_sat_earth_km = r_rel_au * au; 
v_sat_earth_km = v_rel_au * au;

% check if Earth SoI has been reached
dist_from_earth = norm(r_sat_earth_km);

fprintf('\n--- INTERPLANETARY EARTH-EARTH CRUISE REPORT ---\n');
fprintf('Durata viaggio: %.2f giorni\n', (t_vec_cruise_earth_earth(end) - t_vec_cruise_earth_earth(1)) / 86400);
fprintf('Distanza finale dalla Terra: %.2f km\n', dist_from_earth);
fprintf('Raggio SOI Terra: %.2f km\n', soi_earth);

if dist_from_earth < soi_earth
    fprintf('SUCCESS: La sonda è entrata nella SOI della Terra!\n');
else
    fprintf('WARNING: La sonda è arrivata vicina, ma fuori dalla SOI.\n');
    fprintf('Errore di puntamento: %.2f km\n', dist_from_earth - soi_earth);
end

%% Propagate inside Earth SoI during fly-by
% Propagate inside the Earth SoI till the satellite exits the Earth SoI
options_earth_fb = odeset('RelTol', 2.22045e-14, 'AbsTol', 1e-18, 'Events', @(t, y) stopCondition(t, y, soi_earth, 'exit'));
state0_sat_earth_fb = [r_sat_earth_km; v_sat_earth_km]; % punto di partenza simulazione (calcolato prima) 
t_vec_earth_escape = linspace(t_vec_cruise_earth_earth(end), t_vec_cruise_earth_earth(end)+gg*24*60*60, gg*24*60);
[t_vec_earth_escape, state_sat_earth_fb] = ode45(@(t, y) satellite_ode(t, y, mu_earth), t_vec_earth_escape, state0_sat_earth_fb, options_earth_fb);
r_sat_earthfb_escape = state_sat_earth_fb(:, 1:3)';
v_sat_earthfb_escape = state_sat_earth_fb(:, 4:6)';

% Compute jd time and Earth position when satellite is exiting earth SoI limit
jd_earth_sp = t_vec_earth_escape(end)/24/60/60; % momento esatto in cui si ha uscita Earth SoI
earth_fb_soi_date = datetime(jd_earth_sp,'convertfrom','juliandate','Format','d-MMM-y HH:mm:ss', 'TimeZone', timezone);
[~, r_earthfb_sp, v_earthfb_sp] = planet_orbit_coplanar(planets_elements.earth, jd_start, jd_earth_sp, [jd_start, jd_earth_sp]);
r_earthfb_sp = r_earthfb_sp(:, end);  
v_earthfb_sp = v_earthfb_sp(:, end);

% Convert from Earth-Centered to J2000 absolute frame
r_sat_earthfb_sp = r_sat_earthfb_escape(:, end)/au + r_earthfb_sp;
v_sat_earthfb_sp = v_sat_earthfb_escape(:, end)/au + v_earthfb_sp;

% Compute fly-by deltaV [km/s]

deltaV_earthfb = ( norm(v_sat_earthfb_sp) - norm(v_sat_interplanetary_earth_earth) ) * au;

pericenter = norm(r_sat_earthfb_escape(:,end));
fprintf('Pericentro della traiettoria: %.2f km\n', pericenter);
fprintf('Delta V del Fly-by per Terra: %.2f km/s \n', deltaV_earthfb);

plot_mars_soi(t_vec_earth_escape, r_sat_earthfb_escape, soi_earth, v_earthfb_sp, R_earth);

sat.orbit_post_fb_earth = rv2oe(r_sat_earthfb_sp, v_sat_earthfb_sp, mu_sun_au);