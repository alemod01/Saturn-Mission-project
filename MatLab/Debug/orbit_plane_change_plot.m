%% ========================================================================
%  PHASE 1 DEBUG VISUALIZATION: PLANE CHANGE MANEUVER
%  ========================================================================
fprintf('\n--- AVVIO VISUALIZZAZIONE GRAFICA FASE 1 ---\n');
fprintf('Generazione del grafico 3D in corso...\n');
% 1. SETUP DELLA FIGURA
figure('Name', 'Phase 1 Debug: Plane Change', 'Color', 'w', 'Position', [100, 100, 1000, 800]);
hold on; grid on; axis equal;
view(3); % Vista 3D standard
xlabel('X ECI [km]', 'FontWeight', 'bold');
ylabel('Y ECI [km]', 'FontWeight', 'bold');
zlabel('Z ECI [km]', 'FontWeight', 'bold');
title({'Debug Manovra Cambio Piano', 'Sistema di Riferimento Inerziale Geocentrico (ECI)'}, 'FontSize', 14);
% --- A. DISEGNA LA TERRA ---
R_earth_km = 6378.1;
[Xe, Ye, Ze] = sphere(20); % Genera una sfera
surf(Xe*R_earth_km, Ye*R_earth_km, Ze*R_earth_km, ...
    'FaceColor', [0 0.5 1], 'EdgeColor', 'none', 'FaceAlpha', 0.6, 'DisplayName', 'Terra');
% --- B. DISEGNA I PIANI DI RIFERIMENTO ---
% Raggio dei dischi per i piani (un po' più grande dell'orbita)
plane_radius = sat.orbit0.a * 1.3;
theta_plane = linspace(0, 2*pi, 100);
x_circ = plane_radius * cos(theta_plane);
y_circ = plane_radius * sin(theta_plane);
z_zero = zeros(size(x_circ));
% B1. Piano Equatoriale (ECI XY Plane - Z=0) - Colore GRIGIO CHIARO
fill3(x_circ, y_circ, z_zero, [0.8 0.8 0.8], 'FaceAlpha', 0.3, ...
    'EdgeColor', 'k', 'LineStyle', ':', 'DisplayName', 'Piano Equatoriale (i=0°)');
% B2. Piano dell'Eclittica (Tilted Plane) - Colore GIALLO/ARANCIO
% Per visualizzare il piano eclittico nel sistema ECI, dobbiamo ruotare
% il piano XY di +i_ecl gradi attorno all'asse X.
rot_matrix_ecl_vis = [1, 0, 0; 0, cosd(i_ecl), -sind(i_ecl); 0, sind(i_ecl), cosd(i_ecl)];
ecl_points = rot_matrix_ecl_vis * [x_circ; y_circ; z_zero];
fill3(ecl_points(1,:), ecl_points(2,:), ecl_points(3,:), [1 0.7 0.2], 'FaceAlpha', 0.3, ...
    'EdgeColor', 'r', 'LineStyle', '--', 'DisplayName', sprintf('Piano Eclittica (i=%.1f°)', i_ecl));
% --- C. PROPAGAZIONE E DISENGO DELLE ORBITE ---
% Calcoliamo il periodo orbitale per propagare per un giro completo
T_period = 2 * pi * sqrt(sat.orbit0.a^3 / mu_earth);
t_span_vis = linspace(0, T_period, 200); % Vettore tempo per un periodo
options_vis = odeset('RelTol', 1e-5, 'AbsTol', 1e-7); % Tolleranza per la grafica
% C1. Orbita PRE-Manovra (Equatoriale, i=0) - BLU
state0_pre = [r0_eci; v0_eci]; % Stato calcolato all'inizio della PHASE 1
[~, state_traj_pre] = ode45(@(t,y) satellite_ode(t,y,mu_earth), t_span_vis, state0_pre, options_vis);
% Disegna traiettoria BLU
plot3(state_traj_pre(:,1), state_traj_pre(:,2), state_traj_pre(:,3), ...
    'b-', 'LineWidth', 2, 'DisplayName', 'Orbita Pre-Manovra (i=0)');
% Disegna posizione del satellite (Pallino BLU)
plot3(r0_eci(1), r0_eci(2), r0_eci(3), ...
    'bo', 'MarkerSize', 12, 'MarkerFaceColor', 'b', 'DisplayName', 'Sat Pos Pre-Manovra');
% C2. Orbita POST-Manovra (Inclinata, i=23.4) - ROSSA TRATTEGGIATA
% Usiamo lo stato ricalcolato dopo aver aggiornato sat.orbit0.i
state0_post = [r_post_manoeuver; v_post_manoeuver];
[~, state_traj_post] = ode45(@(t,y) satellite_ode(t,y,mu_earth), t_span_vis, state0_post, options_vis);
% Disegna traiettoria ROSSA
plot3(state_traj_post(:,1), state_traj_post(:,2), state_traj_post(:,3), ...
    'r--', 'LineWidth', 2, 'DisplayName', sprintf('Orbita Post-Manovra (i=%.1f)', i_ecl));
% Disegna posizione del satellite (Stella ROSSA)
plot3(r_post_manoeuver(1), r_post_manoeuver(2), r_post_manoeuver(3), ...
    'r*', 'MarkerSize', 14, 'LineWidth', 2, 'DisplayName', 'Sat Pos Post-Manovra');
% --- D. DETTAGLI FINALI ---
% Linea dal centro al punto di manovra per evidenziare il raggio vettore
plot3([0, r0_eci(1)], [0, r0_eci(2)], [0, r0_eci(3)], 'k-', 'LineWidth', 1, 'DisplayName', 'Raggio Vettore al Nodo');
legend('show', 'Location', 'northeastoutside', 'FontSize', 10);
axis vis3d; % Blocca le proporzioni per la rotazione 3D
fprintf('Visualizzazione completata. Ruota il grafico per esplorare.\n');