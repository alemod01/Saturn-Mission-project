%% PLOT 2D: Traiettoria MARS-CENTERED (Fly-by)
try
    % SETUP FIGURA
    figure('Name', 'Fly-by Marte (Mars-Centered)', 'Color', 'w'); clf; hold on
    axis equal; grid on;
    xlabel('X [km]'); ylabel('Y [km]');
    title('Fly-by Marte');

    % DISEGNA MARTE (Corpo fisico) 
    theta_body = linspace(0, 2*pi, 100);
    
    x_body = R_mars * cos(theta_body);
    y_body = R_mars * sin(theta_body);
    
    % Usiamo fill: X, Y, Colore
    fill(x_body, y_body, [0.8 0.4 0.4], ...
         'EdgeColor', 'none', ...
         'DisplayName', 'Marte (Body)');
    
    % SOI MARTE (Confine)
    theta = linspace(0, 2*pi, 200);
    plot(soi_mars * cos(theta), soi_mars * sin(theta), 'm--', 'LineWidth', 1, 'DisplayName', 'Mars SOI');

    % TRAIETTORIA SATELLITE 
    sat_pos_km = state_sat_mars_fb(:, 1:3);
    
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
        v_km_s = v_mars_sp * au;
        
        % Scaliamo la freccia per renderla visibile nel grafico (es. lunga 1/3 della SOI)
        scale_arrow = soi_mars / 2;
        
        % Vettore unitario direzione Marte
        dir_v = v_km_s / norm(v_km_s);
        
        % Disegnamo la freccia partendo dal centro (0,0)
        quiver(0, 0, dir_v(1)*scale_arrow, dir_v(2)*scale_arrow, 0, ...
               'Color', 'm', 'LineWidth', 2, 'MaxHeadSize', 0.5, 'DisplayName', 'Direzione Moto Marte (Elio)');
           
        % B. Velocità Sonda all'ingresso
        v_sat_entry = state_sat_mars_fb(1, 4:6); % Velocità relativa in km/s
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
    fprintf('Altitudine al Periasse: %.2f km\n', min_dist - R_mars);
  
catch ME
    % Corretto l'errore di sintassi nel warning
    warning('MATLAB:PlotError', 'Impossibile creare il plot Mars-Centered: %s', ME.message);
end








