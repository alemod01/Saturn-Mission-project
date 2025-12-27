%% ANIMAZIONE 2D/3D: Fly-by Earth 
try
    % --- SETUP PARAMETRI ---
    sat_pos_km = state_sat_earth_fb(:, 1:3); 
    sat_vel_km = state_sat_earth_fb(:, 4:6);
    
    dists = vecnorm(sat_pos_km, 2, 2);
    min_dist = min(dists); 
    
    % Zoom ampio per vedere la curvatura    
    zoom_mult = 15; 
    visual_radius = min_dist * zoom_mult; 
    
    % --- SCALING  ---
    % Il pianeta non può essere più grande della distanza minima della sonda,
    % altrimenti la sonda ci entra dentro!
    % Calcoliamo il fattore massimo consentito (lasciando un 10% di margine)
    max_safe_scale = (min_dist * 0.90) / R_earth;
    
    % ingrandimento per 3 
    planet_visual_scale = min(3, max_safe_scale); 
    
    % Se la sonda passa vicinissima, scale sarà ~1 (dimensione reale)
    % Se passa lontana, scale sarà 3 (ingrandito per vederlo meglio)
    if planet_visual_scale < 1, planet_visual_scale = 1; end 
    R_scaled = R_earth * planet_visual_scale;
    fprintf('Planet Scale Factor: %.2f (Limitato dalla periapside)\n', planet_visual_scale);
    
    % --- SELEZIONE DATI ---
    idx_inside = find(dists <= visual_radius);
    if isempty(idx_inside), idx_start = 1; idx_end = size(sat_pos_km, 1);
    else
        idx_start = max(1, idx_inside(1) - 10);
        idx_end   = min(size(sat_pos_km, 1), idx_inside(end) + 10);
    end
    
    pos_anim = sat_pos_km(idx_start:idx_end, :);
    vel_anim = sat_vel_km(idx_start:idx_end, :);
    N_anim = size(pos_anim, 1);

    % --- SETUP GRAFICO ---
    figure('Name', 'Earth Fly-by', 'Color', 'k', 'Position', [100, 100, 900, 900]); 
    clf; hold on; axis equal; 
    
    limit_ax = visual_radius * 1.2;
    xlim([-limit_ax, limit_ax]);
    ylim([-limit_ax, limit_ax]);
    zlim([-limit_ax, limit_ax]); 
    
    set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w', 'GridColor', 'w', 'GridAlpha', 0.15);
    xlabel('X [km]', 'Color', 'w'); ylabel('Y [km]', 'Color', 'w'); zlabel('Z [km]', 'Color', 'w');
    title({'KRONOS: Earth Gravity Assist'},'Color', 'w', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'tex');

    % --- EARTH- --
    [Xm, Ym, Zm] = sphere(50); 
    Xm = Xm * R_scaled; 
    Ym = Ym * R_scaled; 
    Zm = Zm * R_scaled;
    
    % Manteniamo 'interp' per avere la sfumatura
    surf(Xm, Ym, Zm, 'EdgeColor', 'none', 'FaceColor', 'interp', 'DisplayName', 'Terra');
   
    n_colors = 100;
    ocean_map = [zeros(n_colors,1), ...                   % Rosso
                 linspace(0, 0.5, n_colors)', ...         % Verde
                 linspace(0.3, 1, n_colors)'];            % Blu
                 
    colormap(gca, ocean_map); 
    caxis([-R_scaled, R_scaled]); 
    
    % Illuminazione
    light('Position', [1 -1 1], 'Style', 'infinite'); 
    lighting gouraud; 
    material shiny;
    
    % VISTA: 2 per 2D, 3 per 3D
    view(2); 

    % --- 5. ELEMENTI DI RIFERIMENTO ---
    theta = linspace(0, 2*pi, 200);
    plot3(visual_radius * cos(theta), visual_radius * sin(theta), zeros(size(theta)),  'b--',  'LineWidth', 1, 'HandleVisibility', 'off');

    % Punti chiave 
    p_in_raw = sat_pos_km(idx_inside(1), :); 
    p_marker_in = (p_in_raw(1:2) / norm(p_in_raw(1:2))) * visual_radius;
    p_out_raw = sat_pos_km(idx_inside(end), :);
    p_marker_out = (p_out_raw(1:2) / norm(p_out_raw(1:2))) * visual_radius;

    plot3(p_marker_in(1), p_marker_in(2), 0, 'co', 'MarkerSize', 7, 'MarkerEdgeColor', 'c', 'LineWidth', 1.5, 'DisplayName', 'Ingresso SOI');
    plot3(p_marker_out(1), p_marker_out(2), 0, 'x', 'MarkerSize', 7, 'MarkerEdgeColor', 'c', 'LineWidth', 1.5, 'DisplayName', 'Uscita SOI');

    % --- ANIMAZIONE ---
    h_probe = plot3(pos_anim(1,1), pos_anim(1,2), pos_anim(1,3), 'o', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'w', 'MarkerSize', 7, 'DisplayName', 'Traiettoria Sonda');
    
    h_trail = animatedline('Color', [0 1 0], 'LineWidth', 2.5, 'HandleVisibility', 'off');
    
    h_text_vel = text(0.5, 0.01, '', 'Units', 'normalized', 'Color', 'w', 'FontSize', 12, 'FontName', 'Courier', 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom'); 

    lgd = legend('Location', 'northeast');
    set(lgd, 'Color', 'k', 'TextColor', 'w', 'EdgeColor', 'w');

    % --- LOOP ---
    target_frames = 1000; 
    step_vis = max(1, floor(N_anim / target_frames));
    
    fprintf('Avvio animazione (%d frames)...\n', floor(N_anim/step_vis));
    
    for i = 1:step_vis:N_anim
        curr_x = pos_anim(i, 1);
        curr_y = pos_anim(i, 2);
        curr_z = pos_anim(i, 3);
        
        set(h_probe, 'XData', curr_x, 'YData', curr_y, 'ZData', curr_z);
        addpoints(h_trail, curr_x, curr_y, curr_z);
        
        set(h_text_vel, 'String', sprintf('Rel Vel: %.2f km/s', norm(vel_anim(i, :))));
        
        drawnow;
    end
    
    % Finale
    set(h_probe, 'XData', pos_anim(end,1), 'YData', pos_anim(end,2), 'ZData', pos_anim(end,3));
    drawnow;

catch ME
    warning('MATLAB:AnimationError','Errore animazione: %s', ME.message);
end