%% PLOT 2D: Traiettoria Earth-CENTERED (Fly-by in km)
try
    % SETUP FIGURA
    figure('Name', 'Fly-by Earth (Earth-Centered)', 'Color', 'w'); clf; hold on
    axis equal; grid on;
    xlabel('X [km]'); ylabel('Y [km]');
    title('Fly-by Earth');

    % DISEGNA Terra (Corpo fisico) 
    theta_body = linspace(0, 2*pi, 100);
    
    x_body = R_earth * cos(theta_body);
    y_body = R_earth * sin(theta_body);
    
    % Usiamo fill: X, Y, Colore
    fill(x_body, y_body, [0 0 1], ...
         'EdgeColor', 'none', ...
         'DisplayName', 'Terra (Body)');
    
    % DISEGNA SOI TERRA (Confine)
    theta = linspace(0, 2*pi, 200);
    plot(soi_earth * cos(theta), soi_earth * sin(theta), 'b--', 'LineWidth', 1, 'DisplayName', 'Earth SOI');

    % TRAIETTORIA SATELLITE 
    sat_pos_km = state_sat_earth_fb(:, 1:3);
    
    plot(sat_pos_km(:,1), sat_pos_km(:,2), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Traiettoria Sonda');

    % PUNTI CHIAVE
    % Ingresso (Primo punto)
    plot(sat_pos_km(1,1), sat_pos_km(1,2), 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Ingresso SOI');
    % Fine Propagazione (Ultimo punto)
    plot(sat_pos_km(end,1), sat_pos_km(end,2), 'kx', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Uscita SOI');

    % VETTORI VELOCITÀ (Quiver)
    try
        % A. Direzione moto di Marte 
        % Conversione au/s -> km/s 
        v_km_s = v_earthfb_sp * au;
        
        % Scaliamo la freccia per renderla visibile nel grafico (es. lunga 1/3 della SOI)
        scale_arrow = soi_earth / 2;
        
        % Vettore unitario direzione Marte
        dir_v = v_km_s / norm(v_km_s);
        
        % Disegnamo la freccia partendo dal centro (0,0)
        quiver(0, 0, dir_v(1)*scale_arrow, dir_v(2)*scale_arrow, 0, ...
               'Color', 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5, 'DisplayName', 'Direzione Moto Terra (Elio)');
           
        % B. Velocità Sonda all'ingresso
        v_sat_entry = state_sat_earth_fb(1, 4:6); % Velocità relativa in km/s
        dir_sat = v_sat_entry / norm(v_sat_entry);
        
        quiver(sat_pos_km(1,1), sat_pos_km(1,2), dir_sat(1)*scale_arrow, dir_sat(2)*scale_arrow, 0, ...
               'Color', 'g', 'LineWidth', 1.5, 'MaxHeadSize', 0.5, 'DisplayName', 'V_{sat} Relativa');
           
    catch
        % Ignora se mancano dati vettoriali
    end
    legend('Location', 'bestoutside');
    
    % ANALISI AL VOLO (Stampa in console)
    dists = vecnorm(sat_pos_km, 2, 2); % Distanza punto per punto
    [min_dist, idx_min] = min(dists);
    fprintf('\n--- FLY-BY REPORT ---\n');
    fprintf('Altitudine al Periasse: %.2f km\n', min_dist - R_earth);
  
catch ME
    % Corretto l'errore di sintassi nel warning
    warning('MATLAB:PlotError', 'Impossibile creare il plot Earth-Centered: %s', ME.message);
end








