%% PLOT 2D: Traiettoria SATURN-CENTERED (Approach + Capture + Moons)
try
    % ---  SETUP FIGURA ---
    figure('Name', 'Saturn Arrival & Capture', 'Color', 'w'); clf; hold on;
    axis equal; grid on;
    xlabel('X [km] (Saturn-Centered)'); ylabel('Y [km] (Saturn-Centered)');
    title({'KRONOS Mission: Saturn Arrival', 'Orbit Insertion & Moons Context'});
    
    % Colore "Giallo Oro" per Saturno e SOI
    saturn_color = [0.92 0.75 0.2]; 
    
    theta_circle = linspace(0, 2*pi, 300);
    
    % ---  DISEGNA SATURNO E ANELLI ---

    % Anelli (Grigio chiaro, per contesto)
    % Raggio interno ~1.2 R_sat, Esterno ~2.3 R_sat
    x_rings_out = (R_saturn * 2.3) * cos(theta_circle);
    y_rings_out = (R_saturn * 2.3) * sin(theta_circle);
    x_rings_in  = (R_saturn * 1.2) * cos(theta_circle);
    y_rings_in  = (R_saturn * 1.2) * sin(theta_circle);
    
    % Disegniamo l'anello esterno
    fill(x_rings_out, y_rings_out, [0.7 0.7 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.4, 'HandleVisibility', 'off');
    fill(x_rings_in, y_rings_in, 'w', 'EdgeColor', 'none', 'HandleVisibility', 'off');
    
    % Corpo di Saturno 
    x_sat = R_saturn * cos(theta_circle);
    y_sat = R_saturn * sin(theta_circle);
    fill(x_sat, y_sat, saturn_color, 'EdgeColor', 'none', 'DisplayName', 'Saturno');
    
    % ---  DISEGNA ORBITE LUNE (Titano & Encelado) ---
    % Orbita Encelado (Linea punteggiata sottile)
    plot(R_Encelado * cos(theta_circle), R_Encelado * sin(theta_circle), ':', ...
         'Color', [0.4 0.4 0.4], 'LineWidth', 0.8, 'DisplayName', 'Orbita Encelado');
     
    % Orbita Titano (Linea tratteggiata)
    plot(R_Titan * cos(theta_circle), R_Titan * sin(theta_circle), '--', ...
         'Color', [0.4 0.4 0.4], 'LineWidth', 0.8, 'DisplayName', 'Orbita Titano');
    
    % ---  DISEGNA SOI SATURNO (Giallo Oro) ---
    plot(soi_saturn * cos(theta_circle), soi_saturn * sin(theta_circle), '--', ...
         'Color', saturn_color, 'LineWidth', 1.5, 'DisplayName', 'Saturn SOI');
    

    % ---  TRAIETTORIE ---
    % A. Fase di Avvicinamento (Verde)
    r_approach = state_sat_saturn_soi(:, 1:3);
    plot(r_approach(:,1), r_approach(:,2), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Approach Trajectory');
    
    % B. Orbita di Cattura (Blu)
    r_capture = state_final(:, 1:3);
    plot(r_capture(:,1), r_capture(:,2), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Final Captured Orbit');
    
    % ---  PUNTI CHIAVE & MARKER ---
    % Ingresso SOI
    plot(r_approach(1,1), r_approach(1,2), 'mo', 'MarkerFaceColor', 'm', 'MarkerSize', 6, 'DisplayName', 'SOI Entry');
    
    % Punto di Manovra (Insertion Burn)
    r_maneuver_sat = r_approach(end, 1:3);
    plot(r_maneuver_sat(1), r_maneuver_sat(2), 'p', 'MarkerSize', 14, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'DisplayName', 'Insertion Burn (\DeltaV)');
     
    % ---  VETTORI VELOCITÀ ---
    try
        scale_arrow = soi_saturn / 3;
       
        % Velocità Saturno
        v_saturn_vec = v_saturn_sp; % [Au/s]
        dir_saturn = v_saturn_vec / norm(v_saturn_vec);
        quiver(0, 0, dir_saturn(1)*scale_arrow, dir_saturn(2)*scale_arrow, 0, 'Color', [0.8 0.6 0.2], 'LineWidth', 2, 'MaxHeadSize', 0.5, 'DisplayName', 'V_{Saturn} (Helio)');
        
        % Velocità Sonda all'ingresso
        v_entry_sat = state_sat_saturn_soi(1, 4:6); %[km/s]
        dir_entry = v_entry_sat / norm(v_entry_sat);
        quiver(r_approach(1,1), r_approach(1,2), dir_entry(1)*scale_arrow, dir_entry(2)*scale_arrow, 0, 'Color', 'g', 'LineWidth', 1.5, 'MaxHeadSize', 0.5, 'DisplayName', 'V_{entry} Relativa');
    catch
    end
    
    legend('Location', 'bestoutside');
    
    % ---  ZOOM (Opzionale) ---
    % Se vuoi vedere bene le lune e l'orbita, scommenta la riga sotto per limitare gli assi.
    %xlim([-R_Titan*1.5, R_Titan*1.5]); ylim([-R_Titan*1.5, R_Titan*1.5]); % Focus su Titano
    

catch ME
    warning('MATLAB:PlotError', 'Impossibile creare il plot Saturn-Centered: %s', ME.message);
end