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
mu_sun = 1.32712440018e11;  % [km^3/s^2]
mu_earth = 3.986004418e5;   % [km^3/s^2]
mu_mars = 42828.37;         % [km^3/s^2]
mu_jupiter = 1.26686534e8;  % [km^3/s^2]
mu_saturn = 3.7931187e7;    % [km^3/s^2]

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
start_date = datetime('2034-08-25 12:00:00', "TimeZone", timezone); 
earth_fb_date = datetime('2036-11-10 12:00:00', "TimeZone", timezone);
mars_fb_date = datetime('2040-04-14 12:00:00', "TimeZone", timezone);
saturn_arrival_date = datetime('2046-09-15 12:00:00', "TimeZone", timezone);
end_date = datetime('2046-11-01 12:00:00', "TimeZone", timezone);

jd_start = juliandate(start_date);
jd_earth_fb = juliandate(earth_fb_date);
jd_mars_fb = juliandate(mars_fb_date);
jd_saturn_arrival = juliandate(saturn_arrival_date);
jd_end = juliandate(end_date);

% Select planets to visualize
planets = {'venus', 'earth', 'mars', 'jupiter', 'saturn'};

% Compute planets orbital elements at starting epoch time, and consider
% them constant throughout the entire mission (except for the mean anomaly)
planets_elements = elements_from_ephems(planets, jd_start);

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
fprintf('================ EARTH ESCAPE DIAGNOSTICS =================\n');
fprintf('||v_inf||              = %.3f km/s\n', norm(v_inf_kms));
fprintf('||v_c (parking)||      = %.3f km/s\n', v_c);
fprintf('||v_perigee_hyper||    = %.3f km/s\n', v_hyperbola_perigee_mag);
fprintf('DeltaV escape          = %.3f km/s\n', sat.deltaV_escape_earth);
fprintf('===========================================================\n');


% -------------------------------------------------------------------------
% CORREZIONE MANUALE DEL TIRO (TARGETING)
% -------------------------------------------------------------------------
% Il problema: Partiamo dalla SOI, non dal centro della Terra.
% La soluzione: Correggiamo leggermente la velocità e l'angolo di uscita.

% 1. Definiamo i fattori di correzione

% Correzioni fly by Terra
k_vel = 0.997339;       % Moltiplicatore di velocità (es. 0.999 o 1.001)
delta_angle = -0.1334401;   % Correzione angolo in gradi (es. +0.5 o -0.5)


% 2. Applichiamo la correzione alla Magnitudine
v_old = v_esc; 
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

deltaV_cost_earth_esc_km= norm (v_esc - v_old); %[Km/s]
% (Stampa di debug per vedere cosa stiamo facendo)
fprintf('\n============== TARGETING CORRECTION APPLIED ===============\n');
fprintf('Earth Departure Corrections\n');
fprintf('Velocità scalata di: %.4f\n', k_vel);
fprintf('Angolo ruotato di:   %.2f deg\n', delta_angle);
fprintf('-----------------------------------------------------------\n');
fprintf('COSTO MANOVRA (Delta V): %.4f km/s\n', deltaV_cost_earth_esc_km);
fprintf('===========================================================\n');


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


fprintf('\n============= EARTH-EARTH CRUISE PHASE REPORT =============\n');
fprintf('Time of Flight: %.2f days\n', (t_vec_cruise_earth_earth(end) - t_vec_cruise_earth_earth(1)) / 86400);
if dist_from_earth < soi_earth
    fprintf('SUCCESS: Spacecraft successfully entered Earth''s SOI!\n');
else
    fprintf('WARNING: Spacecraft arrived outside the SOI.\n');
    fprintf('Miss distance from SOI boundary: %.2f km\n', dist_from_earth - soi_earth);
    return; 
end
fprintf('===========================================================\n');


%% PROPAGATE FLY-BY EARTH
% Propagate inside the Earth SoI till the satellite exits the Earth SoI
state0_sat_earth_fb = [r_sat_earth_km; v_sat_earth_km]; % punto di partenza simulazione (calcolato prima) 

t_vec_earth_escape = linspace(t_vec_cruise_earth_earth(end), t_vec_cruise_earth_earth(end)+gg*24*60*60, gg*24*60);
ode_stop_mode_fb_Earth='exit' ;  % Select stop mode: 'exit' or 'pericenter'
options_earth_fb = odeset('RelTol', 2.22045e-14, 'AbsTol', 1e-18, 'Events', @(t, y) stopCondition(t, y, soi_earth, ode_stop_mode_fb_Earth));
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


fprintf('\n=================== FLY BY EARTH REPORT ===================\n');
if strcmp(ode_stop_mode_fb_Earth, 'exit')
    dists_earthfb = vecnorm(state_sat_earth_fb(:, 1:3), 2, 2); % Distanza punto per punto
    [min_dist_earthfb, ~] = min(dists_earthfb);
    fprintf('Closest Approach: %.2f km from Earth surface\n', min_dist_earthfb - R_earth);
    fprintf('Earths Flyby Delta-v: %.2f km/s\n', deltaV_earthfb);
    v_exit_earthfb_helio_km = norm(v_sat_earthfb_sp) * au; %[km/s]
    fprintf('Spacecraft Heliocentric Exit Velocity: %.2f km/s\n', v_exit_earthfb_helio_km);
elseif strcmp(ode_stop_mode_fb_Earth, 'pericenter')
    fprintf('Simulation terminated at Pericenter (Closest Approach).\n');
    fprintf('Altitude: %.2f km from Earth surface\n', norm(r_sat_earthfb_escape(:, end)) - R_earth);
    return; 
end
fprintf('===========================================================\n');
plot_flyBy(t_vec_earth_escape, r_sat_earthfb_escape, soi_earth, v_earthfb_sp, R_earth);


% -------------------------------------------------------------------------
% CORREZIONE MANUALE DEL TIRO (TARGETING) - POST EARTH FLYBY
% -------------------------------------------------------------------------
% Correggiamo leggermente la velocità e l'angolo di uscita dalla SOI della
% Terra
% per compensare eventuali errori di puntamento verso la Terra
% 1. Definiamo i fattori di correzione
 
k_vel_earthfb = 0.972249;        % Moltiplicatore di velocità (es. 0.999 o 1.001)
delta_angle_earthfb = 0.01308;   % Correzione angolo in gradi (es. +0.5 o -0.5)

% 2. Applichiamo la correzione alla Magnitudine
v_old = v_sat_earthfb_sp; 
v_earthfb_sp_mag_corr = norm(v_sat_earthfb_sp) * k_vel_earthfb;

% 3. Applichiamo la correzione alla Direzione
% Direzione attuale della velocità
v_earthfb_sp_dir = v_sat_earthfb_sp / norm(v_sat_earthfb_sp);

% Ruotiamo il vettore direzione di 'delta_angle_mars' gradi nel piano dell'Eclittica
theta_corr_mars = deg2rad(delta_angle_earthfb);
R_corr_mars = [cos(theta_corr_mars), -sin(theta_corr_mars), 0;
               sin(theta_corr_mars),  cos(theta_corr_mars), 0;
               0,                     0,                    1];
      
v_earthfb_sp_dir_corr = (R_corr_mars * v_earthfb_sp_dir)'; % Ruotiamo

% 4.Ricalcoliamo il vettore di velocità finale CORRETTO
v_sat_earthfb_sp = (v_earthfb_sp_mag_corr * v_earthfb_sp_dir_corr)';

deltaV_cost_earth_fb_au= norm (v_sat_earthfb_sp - v_old); %[Au/s]
deltaV_cost_earth_fb_km= deltaV_cost_earth_fb_au*au; %[Km/s]

fprintf('\n============== TARGETING CORRECTION APPLIED ===============\n');
fprintf('Post fly by Earth Departure Corrections\n');
fprintf('Velocità scalata di: %.4f\n', k_vel_earthfb);
fprintf('Angolo ruotato di:   %.2f deg\n', delta_angle_earthfb);
fprintf('-----------------------------------------------------------\n');
fprintf('COSTO MANOVRA (Delta V): %.4f km/s\n', deltaV_cost_earth_fb_km);
fprintf('===========================================================\n');

sat.orbit_post_fb_earth = rv2oe(r_sat_earthfb_sp, v_sat_earthfb_sp, mu_sun_au);

%% Propagate from outside Earth SoI to Mars SoI - interplanetary 
% Propagate from outside Earth SoI till the satellite enters the Mars SoI

state0_interplanetary_earth_mars = [r_sat_earthfb_sp; v_sat_earthfb_sp];

% parametro correttivo gg interplanetary
kg = 0;
% calcolo il tempo per arrivare da fuori SoI Terra a dentro SoI marte 
t_cruise_total_earth_mars = jd_mars_fb*24*60*60 - t_vec_earth_escape(end); %durata in secondi del viaggio 

t_cruise_total_earth_mars=t_cruise_total_earth_mars+kg*24*60*60;
t_vec_cruise_earth_mars= linspace(t_vec_earth_escape(end), jd_mars_fb*24*60*60  + (kg*24*60*60) ,t_cruise_total_earth_mars/3600);
% propagazione
options_cruise_earth_mars = odeset('RelTol', 2.22045e-14, 'AbsTol', 1e-18, 'Events', @(t, y) stopCondition_interplanetary(t, y, jd_start , planets_elements.mars , soi_mars));
[t_vec_cruise_earth_mars, state_cruise_earth_mars] = ode45(@(t, y) satellite_ode(t, y, mu_sun_au), t_vec_cruise_earth_mars, state0_interplanetary_earth_mars, options_cruise_earth_mars);

% Estrazione Risultati finali (Stato Eliocentrico all'arrivo)
r_sat_interplanetary_earth_mars = state_cruise_earth_mars(end, 1:3)'; % Posizione rispetto al Sole (AU)
v_sat_interplanetary_earth_mars = state_cruise_earth_mars(end, 4:6)'; % Velocità rispetto al Sole (AU/s)

% Calcolo Data Esatta di Arrivo a SoI marte (Julian Date)
jd_mars_arrival_actual = t_vec_cruise_earth_mars(end)/24/60/60;

% Dove si trova Marte in quel momento preciso?
% Calcoliamo posizione e velocità di Marte usando le effemeridi
[~, r_mars_arr, v_mars_arr] = planet_orbit_coplanar(planets_elements.mars, jd_start, jd_mars_arrival_actual, [jd_start, jd_mars_arrival_actual]);
r_mars_arr = r_mars_arr(:, end); % Posizione Marte (AU)
v_mars_arr = v_mars_arr(:, end); % Velocità Marte (AU/s)


% Convert from J2000 absolute frame to Mars-Centered frame
% Calcoliamo il vettore relativo (Sonda rispetto a Marte)
r_rel_au = r_sat_interplanetary_earth_mars - r_mars_arr; % Vettore distanza in AU
v_rel_au = v_sat_interplanetary_earth_mars - v_mars_arr; % Vettore velocità relativa in AU/s

% Convertiamo in km e km/s per il controllo finale
r_sat_mars_km = r_rel_au * au; 
v_sat_mars_km = v_rel_au * au;

% heck if Mars SoI has been reached
dist_from_mars = norm(r_sat_mars_km);

fprintf('\n============= EARTH-MARS CRUISE PHASE REPORT ==============\n');
fprintf('Time of Flight: %.2f days\n', (t_vec_cruise_earth_mars(end) - t_vec_cruise_earth_mars(1)) / 86400);
fprintf('Distanza finale da Marte: %.2f km\n', dist_from_mars);
fprintf('Raggio SOI Marte: %.2f km\n', soi_mars);

if dist_from_mars < soi_mars
    fprintf('SUCCESS: Spacecraft successfully entered Earth''s SOI!\n');
else
    fprintf('WARNING: Spacecraft arrived outside the SOI.\n');
    fprintf('Miss distance from SOI boundary: %.2f km\n', dist_from_mars - soi_mars);
    return; 
end
fprintf('===========================================================\n');


%% PROPAGATE FLY-BY MARS
% Propagate inside the Mars SoI till the satellite exits the Mars SoI
state0_sat_mars_escape = [r_sat_mars_km; v_sat_mars_km]; % punto di partenza simulazione (calcolato prima)

t_vec_mars_escape = linspace(t_vec_cruise_earth_mars(end), t_vec_cruise_earth_mars(end)+gg*24*60*60, gg*24*60);
ode_stop_mode= 'exit'; %scegli tra 'exit' o 'pericenter'
options_mars_fb = odeset('RelTol', 2.22045e-14, 'AbsTol', 1e-18, 'Events', @(t, y) stopCondition(t, y, soi_mars, ode_stop_mode)); 
[t_vec_mars_escape, state_sat_mars_fb] = ode45(@(t, y) satellite_ode(t, y, mu_mars), t_vec_mars_escape, state0_sat_mars_escape, options_mars_fb);
r_sat_mars_escape = state_sat_mars_fb(:, 1:3)';
v_sat_mars_escape = state_sat_mars_fb(:, 4:6)';

% Compute jd time and Mars position when satellite is exiting Mars SoI limit
jd_mars_sp = t_vec_mars_escape(end)/24/60/60; % momento esatto in cui si ha uscita Earth SoI
mars_soi_date = datetime(jd_mars_sp,'convertfrom','juliandate','Format','d-MMM-y HH:mm:ss', 'TimeZone', timezone);
[~, r_mars_sp, v_mars_sp] = planet_orbit_coplanar(planets_elements.mars, jd_start, jd_mars_sp, [jd_start, jd_mars_sp]);
r_mars_sp = r_mars_sp(:, end);  
v_mars_sp = v_mars_sp(:, end);

% Convert from Mars-Centered to J2000 absolute frame
r_sat_marsfb_sp = r_sat_mars_escape(:, end)/au + r_mars_sp;
v_sat_marsfb_sp = v_sat_mars_escape(:, end)/au + v_mars_sp;

% Compute fly-by deltaV [km/s]
deltaV_marsfb = ( norm(v_sat_marsfb_sp) - norm(v_sat_interplanetary_earth_mars) ) * au;

fprintf('\n=================== FLY BY MARS REPORT ====================\n');
if strcmp(ode_stop_mode, 'exit')
    dists = vecnorm(state_sat_mars_fb(:, 1:3), 2, 2); % Distanza punto per punto
    [min_dist, idx_min] = min(dists);
    fprintf('Closest Approach: %.2f km from Mars surface\n', min_dist - R_mars);
    fprintf('Mars Flyby Delta-v: %.2f km/s\n', deltaV_marsfb);
    v_exit_marsfb_helio_km = norm(v_sat_marsfb_sp) * au; %[km/s]
    fprintf('Spacecraft Heliocentric Exit Velocity: %.2f km/s\n', v_exit_marsfb_helio_km);
elseif strcmp(ode_stop_mode, 'pericenter')
    fprintf('Simulation terminated at Pericenter (Closest Approach).\n');
    fprintf('Altitude: %.2f km from Mars surface\n', norm(r_sat_mars_escape(:, end)) - R_mars);
    return; 
end
fprintf('===========================================================\n');
plot_flyBy(t_vec_mars_escape, r_sat_mars_escape, soi_mars, v_mars_sp, R_mars);

% -------------------------------------------------------------------------
% CORREZIONE MANUALE DEL TIRO (TARGETING) - POST MARS FLYBY
% -------------------------------------------------------------------------
% Correggiamo leggermente la velocità e l'angolo di uscita dalla SOI di Marte
% per compensare eventuali errori di puntamento verso la Terra
% 1. Definiamo i fattori di correzione
 k_vel_marsfb = 1.11835;        % Moltiplicatore di velocità (es. 0.999 o 1.001)
 delta_angle_earthfb = -0.1463;   % Correzione angolo in gradi (es. +0.5 o -0.5)

% 2. Applichiamo la correzione alla Magnitudine
v_old = v_sat_marsfb_sp;
v_marsfb_sp_mag_corr = norm(v_sat_marsfb_sp ) * k_vel_marsfb;

% 3. Applichiamo la correzione alla Direzione
% Direzione attuale della velocità
v_marsfb_sp_dir = v_sat_marsfb_sp  / norm(v_sat_marsfb_sp);

% Ruotiamo il vettore direzione di 'delta_angle_mars' gradi nel piano dell'Eclittica
theta_corr_mars = deg2rad(delta_angle_earthfb);
R_corr_mars = [cos(theta_corr_mars), -sin(theta_corr_mars), 0;
               sin(theta_corr_mars),  cos(theta_corr_mars), 0;
               0,                     0,                    1];
      
v_marsfb_sp_dir_corr = (R_corr_mars * v_marsfb_sp_dir)'; % Ruotiamo
 
% 4.Ricalcoliamo il vettore di velocità finale CORRETTO
v_sat_marsfb_sp = (v_marsfb_sp_mag_corr * v_marsfb_sp_dir_corr)';

deltaV_cost_mars_fb_au= norm (v_sat_marsfb_sp - v_old); %[Au/s]
deltaV_cost_mars_fb_km= deltaV_cost_mars_fb_au*au; %[Km/s]

fprintf('\n============== TARGETING CORRECTION APPLIED ===============\n');
fprintf('Post fly by Earth Departure Corrections\n');
fprintf('Velocità scalata di: %.4f\n', k_vel_marsfb);
fprintf('Angolo ruotato di:   %.2f deg\n', delta_angle_earthfb);
fprintf('-----------------------------------------------------------\n');
fprintf('COSTO MANOVRA (Delta V): %.4f km/s\n', deltaV_cost_mars_fb_km);
fprintf('===========================================================\n');
sat.orbit_post_mars_fb = rv2oe(r_sat_marsfb_sp, v_sat_marsfb_sp, mu_sun_au);

%% Propagate from outside Mars SoI to Saturn SoI - interplanetary 
% Propagate from outside Mars SoI till the satellite enters the Saturn SoI

state0_interplanetary_mars_saturn = [r_sat_marsfb_sp; v_sat_marsfb_sp];

% parametro correttivo gg interplanetary
kg = 0;
% calcolo il tempo per arrivare da fuori SoI Marte a dentro SoI Saturno 
t_cruise_total_mars_saturn = jd_saturn_arrival*24*60*60 - t_vec_mars_escape(end); %durata in secondi del viaggio 

t_cruise_total_mars_saturn=t_cruise_total_mars_saturn+kg*24*60*60;
t_vec_cruise_mars_saturn= linspace(t_vec_mars_escape(end), jd_saturn_arrival*24*60*60  + (kg*24*60*60) ,t_cruise_total_mars_saturn/3600);
% propagazione
options_cruise_mars_saturn = odeset('RelTol', 2.22045e-14, 'AbsTol', 1e-18, 'Events', @(t, y) stopCondition_interplanetary(t, y, jd_start , planets_elements.saturn , soi_saturn));
[t_vec_cruise_mars_saturn, state_cruise_mars_saturn] = ode45(@(t, y) satellite_ode(t, y, mu_sun_au), t_vec_cruise_mars_saturn, state0_interplanetary_mars_saturn, options_cruise_mars_saturn);

% Estrazione Risultati finali (Stato Eliocentrico all'arrivo)
r_sat_interplanetary_mars_saturn = state_cruise_mars_saturn(end, 1:3)'; % Posizione rispetto al Sole (AU)
v_sat_interplanetary_mars_saturn = state_cruise_mars_saturn(end, 4:6)'; % Velocità rispetto al Sole (AU/s) V_sa

% Calcolo Data Esatta di Arrivo a SoI marte (Julian Date)
jd_saturn_arrival_actual = t_vec_cruise_mars_saturn(end)/24/60/60;

% Dove si trova Marte in quel momento preciso?
% Calcoliamo posizione e velocità di Marte usando le effemeridi
[~, r_saturn_arr, v_saturn_arr] = planet_orbit_coplanar(planets_elements.saturn, jd_start, jd_saturn_arrival_actual, [jd_start, jd_saturn_arrival_actual]);
r_saturn_arr = r_saturn_arr(:, end); % Posizione Marte (AU)
v_saturn_arr = v_saturn_arr(:, end); % Velocità Marte (AU/s)


% Convert from J2000 absolute frame to Saturn-Centered frame
% Calcoliamo il vettore relativo (Sonda rispetto a Saturno)
r_rel_au = r_sat_interplanetary_mars_saturn - r_saturn_arr; % Vettore distanza in AU
v_rel_au = v_sat_interplanetary_mars_saturn - v_saturn_arr; % Vettore velocità relativa in AU/s

% Convertiamo in km e km/s per il controllo finale
r_sat_saturn_km = r_rel_au * au; 
v_sat_saturn_km = v_rel_au * au;

% heck if Mars SoI has been reached
dist_from_saturn = norm(r_sat_saturn_km);

fprintf('\n============= MARS-SATURN CRUISE PHASE REPORT ==============\n');
fprintf('Time of Flight: %.2f days\n', (t_vec_cruise_mars_saturn(end) - t_vec_cruise_mars_saturn(1)) / 86400);

if dist_from_saturn < soi_saturn
    fprintf('SUCCESS: Spacecraft successfully entered Saturn''s SOI!\n');
else
    fprintf('WARNING: Spacecraft arrived outside the SOI.\n');
    fprintf('Miss distance from SOI boundary: %.2f km\n', dist_from_saturn - soi_saturn);
    return; 
end
fprintf('===========================================================\n');

%% PROPAGATE INSIDE SATURN SOI
% Propagate inside the Saturn SoI till the satellite reaches the desired
% orbit
state0_sat_saturn_soi = [r_sat_saturn_km; v_sat_saturn_km]; % punto di partenza simulazione (calcolato prima)

gg = 200;
t_vec_saturn_soi = linspace(t_vec_cruise_mars_saturn(end), t_vec_cruise_mars_saturn(end)+gg*24*60*60, gg*24*60);
ode_stop_mode = 'pericenter'; %scegli tra 'exit' o 'pericenter'
options_saturn_soi = odeset('RelTol', 2.22045e-14, 'AbsTol', 1e-18, 'Events', @(t, y) stopCondition(t, y, soi_saturn, ode_stop_mode)); 
[t_vec_saturn_soi, state_sat_saturn_soi] = ode45(@(t, y) satellite_ode(t, y, mu_saturn), t_vec_saturn_soi, state0_sat_saturn_soi, options_saturn_soi);
r_sat_saturn_soi = state_sat_saturn_soi(:, 1:3)';
v_sat_saturn_soi = state_sat_saturn_soi(:, 4:6)';

% Compute jd time and Mars position when satellite is exiting Mars SoI limit
jd_saturn_sp = t_vec_saturn_soi(end)/24/60/60; % momento esatto in cui si ha uscita Earth SoI
saturn_soi_date = datetime(jd_saturn_sp,'convertfrom','juliandate','Format','d-MMM-y HH:mm:ss', 'TimeZone', timezone);
[~, r_saturn_sp, v_saturn_sp] = planet_orbit_coplanar(planets_elements.saturn, jd_start, jd_saturn_sp, [jd_start, jd_saturn_sp]);
r_saturn_sp = r_saturn_sp(:, end);  
v_saturn_sp = v_saturn_sp(:, end);

% Convert from saturn-Centered to J2000 absolute frame
r_sat_saturn_sp = r_sat_saturn_soi(:, end)/au + r_saturn_sp;
v_sat_saturn_sp = v_sat_saturn_soi(:, end)/au + v_saturn_sp;

% Compute fly-by deltaV [km/s]
deltaV_saturn = ( norm(v_sat_saturn_sp) - norm(v_sat_interplanetary_mars_saturn) ) * au;
deltaV_inside_soi = norm(v_sat_saturn_soi(:, end)) - norm(v_sat_saturn_km);

fprintf('\n=================== SATURN SOI REPORT ====================\n');
if strcmp(ode_stop_mode, 'exit')
    dists = vecnorm(state_sat_saturn_soi(:, 1:3), 2, 2); % Distanza punto per punto
    [min_dist, idx_min] = min(dists);
    fprintf('Closest Approach: %.2f km from Saturn surface\n', min_dist - R_saturn);
    fprintf('Saturn Flyby Delta-v: %.2f km/s\n', deltaV_saturn);
    v_exit_saturn_helio_km = norm(v_sat_saturn_sp) * au; %[km/s]
    fprintf('Spacecraft Heliocentric Exit Velocity: %.2f km/s\n', v_exit_saturn_helio_km);
elseif strcmp(ode_stop_mode, 'pericenter')
    fprintf('Simulation terminated at Pericenter (Closest Approach).\n');
    fprintf('Altitude: %.2f km from Saturn surface\n', norm(r_sat_saturn_soi(:, end)) - R_saturn);
    fprintf('Speed at pericenter compared to Saturn: %.2f km/s \n', norm(v_sat_saturn_soi(:, end)));
    fprintf('DeltaV at pericenter compared to Saturn: %.2f km/s \n', deltaV_inside_soi);
    fprintf('Actual ending date of the mission: %s \n', saturn_soi_date);
end
fprintf('===========================================================\n');
plot_flyBy(t_vec_saturn_soi, r_sat_saturn_soi, soi_saturn, v_saturn_sp, R_saturn);

%% Capture Maneuver
% Manovra per rimanere in orbita circolare intorno a Saturno, applicata
% quando arriviamo nel punto più vicino a Saturno nella propagazione
% precedente
R_Titan = 1222000;          % km
R_Encelado = 238000;        % km
e_des = (R_Titan - R_Encelado)/(R_Titan + R_Encelado) + 0.1;
v_cattura = sqrt( (mu_saturn/norm(r_sat_saturn_soi(:,end))) * (1 + e_des));  % velocità necessaria per rimanere in orbita circolare
% Calcolo del Delta-V necessario
% Dobbiamo frenare: DeltaV = V_attuale - V_necessaria
deltaV_capture = v_cattura - norm(v_sat_saturn_soi(:,end));
% 4. Calcoliamo il nuovo vettore velocità post-manovra
% La direzione rimane la stessa (tangente all'orbita), cambia solo il modulo
v_p_dir = v_sat_saturn_soi(:,end) / norm(v_sat_saturn_soi(:,end)); % Versore direzione velocità
v_post_maneuver = v_cattura * v_p_dir;
sat.orbit_post_capture = rv2oe(r_sat_saturn_soi(:,end), v_post_maneuver, mu_saturn);
% Se vuoi propagare l'orbita circolare per vederla graficamente:
options_final_orbit = odeset('RelTol', 2.22045e-14, 'AbsTol', 1e-18);
period_saturn = 2 * pi * sqrt(norm(r_sat_saturn_soi(:,end))^3 / mu_saturn);
t_orbit = linspace(0, 10*period_saturn, 1000); % Propaghiamo per 2 periodi
    
state0_final = [r_sat_saturn_soi(:,end); v_post_maneuver];
[~, state_final] = ode45(@(t, y) satellite_ode(t, y, mu_saturn), t_orbit, state0_final, options_final_orbit);
    
% Plot dell'orbita di cattura (aggiungi al plot esistente o nuovo)
% Usiamo la tua funzione plot_flyBy o plot manuale
figure('Name', 'Saturn Capture Orbit');
plot3(state_final(:,1), state_final(:,2), state_final(:,3), 'g', 'LineWidth', 2);
hold on; grid on; axis equal;
[xs, ys, zs] = sphere(50);
surf(xs*R_saturn, ys*R_saturn, zs*R_saturn, 'FaceColor', [0.8 0.6 0.4]); % Saturno
% --- Orbite circolari di Encelado e Titano (assunte nel piano XY, z=0) ---
theta = linspace(0, 2*pi, 1000);
% Encelado
x_enc = R_Encelado * cos(theta);
y_enc = R_Encelado * sin(theta);
z_enc = zeros(size(theta));
plot3(x_enc, y_enc, z_enc, 'b--', 'LineWidth', 1.5);
% Titano
x_tit = R_Titan * cos(theta);
y_tit = R_Titan * sin(theta);
z_tit = zeros(size(theta));
plot3(x_tit, y_tit, z_tit, 'r--', 'LineWidth', 1.5);
% Aggiorna legenda (include anche le orbite)
title('Orbita Finale di Cattura attorno a Saturno');
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');
legend('Orbita Finale', 'Saturno', 'Orbita Encelado', 'Orbita Titano');