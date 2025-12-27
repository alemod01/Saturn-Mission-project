%% ANIMAZIONE 2D/3D: Saturn Arrival & Orbit Capture (SAMPLED SIMPLE)
try
    % ---  PREPARAZIONE DATI (Unione Approccio + Orbita) ---
    pos_app = state_sat_saturn_soi(:, 1:3);
    vel_app = state_sat_saturn_soi(:, 4:6);
    pos_orb = state_final(:, 1:3);
    vel_orb = state_final(:, 4:6);
    
    % Totale
    pos_total = [pos_app; pos_orb];
    vel_total = [vel_app; vel_orb];
    
    dists = vecnorm(pos_total, 2, 2);
    min_dist = min(dists); 
    
    % --- ZOOM & SCALING ---
    zoom_mult = 15; 
    visual_radius = max(min_dist * zoom_mult, R_Titan * 1.1); 
    
    max_safe_scale = (min_dist * 0.90) / R_saturn;
    planet_visual_scale = min(3, max_safe_scale);

    if planet_visual_scale < 1
        planet_visual_scale = 1; 
    end 
    R_scaled = R_saturn * planet_visual_scale;
    %fprintf('Planet Scale Factor: %.2f\n', planet_visual_scale);
    
    % ---  SELEZIONE DATI E CAMPIONAMENTO ---
    idx_inside = find(dists <= visual_radius);
    
    if isempty(idx_inside)
        idx_start = 1; 
        idx_end = size(pos_total, 1);
    else
        idx_start = max(1, idx_inside(1) - 20);
        idx_end   = min(size(pos_total, 1), idx_inside(end) + 20);
    end
    
    % Ritaglio sulla zona visibile
    pos_cropped = pos_total(idx_start:idx_end, :);
    vel_cropped = vel_total(idx_start:idx_end, :);
    
    k_step = 15; 
    pos_anim = pos_cropped(1:k_step:end, :);
    vel_anim = vel_cropped(1:k_step:end, :);
    N_anim = size(pos_anim, 1);
    

    % ---  SETUP GRAFICO ---
    figure('Name', 'Saturn Capture Animation', 'Color', 'k', 'Position', [100, 100, 900, 900]); 
    clf; hold on; axis equal; 
    
    limit_ax = visual_radius * 1.1;
    xlim([-limit_ax, limit_ax]);
    ylim([-limit_ax, limit_ax]);
    zlim([-limit_ax, limit_ax]); 
    
    set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w', 'GridColor', 'w', 'GridAlpha', 0.15);
    xlabel('X [km]', 'Color', 'w'); ylabel('Y [km]', 'Color', 'w'); zlabel('Z [km]', 'Color', 'w');
    title({'KRONOS: Saturn Orbit Insertion', '\fontsize{10}\it(Approach + Parking Orbit)'},'Color', 'w', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'tex');

    % ---  SATURNO ---
    [Xm, Ym, Zm] = sphere(50); 
    Xm = Xm * R_scaled; Ym = Ym * R_scaled; Zm = Zm * R_scaled;
    
    surf(Xm, Ym, Zm, 'EdgeColor', 'none', 'FaceColor', 'interp', 'DisplayName', 'Saturno');
   
    n_colors = 100;
    gold_map = [linspace(0.2, 0.92, n_colors)', linspace(0.1, 0.75, n_colors)', linspace(0.0, 0.2, n_colors)'];
    colormap(gca, gold_map); 
    caxis([-R_scaled, R_scaled]); 
    
    light('Position', [1 -1 0.5], 'Style', 'infinite'); 
    lighting gouraud; 
    material dull; 
    
    % ---  ANELLI E LUNE (Statici) ---
    theta_static = linspace(0, 2*pi, 200);
    z_static = zeros(size(theta_static));
     
    % Anelli
    r_in_vis  = (R_scaled * 1.2);
    r_out_vis = (R_scaled * 2.3);
    
    x_r_out = r_out_vis * cos(theta_static); y_r_out = r_out_vis * sin(theta_static);
    x_r_in  = r_in_vis * cos(theta_static);  y_r_in  = r_in_vis * sin(theta_static);
     
    fill3(x_r_out, y_r_out, z_static, [0.7 0.7 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.3, 'HandleVisibility', 'off');
    fill3(x_r_in, y_r_in, z_static, 'k', 'EdgeColor', 'none', 'HandleVisibility', 'off'); 
       
    % Lune
    if ~exist('R_Titan', 'var'), R_Titan = 1222000; end
    if ~exist('R_Encelado', 'var'), R_Encelado = 238000; end
    
    plot3(R_Titan * cos(theta_static), R_Titan * sin(theta_static), z_static, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'DisplayName', 'Orbita Titano');
    plot3(R_Encelado * cos(theta_static), R_Encelado * sin(theta_static), z_static, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.8, 'DisplayName', 'Orbita Encelado');

    % C. Marker Manovra
    plot3(pos_app(end,1), pos_app(end,2), pos_app(end,3), 'p', 'MarkerSize', 12, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'w', 'DisplayName', 'Insertion Burn');
    
    view(2); 
    
    % ---  ANIMAZIONE ---
    h_probe = plot3(pos_anim(1,1), pos_anim(1,2), pos_anim(1,3), 'o', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'w', 'MarkerSize', 6, 'DisplayName', 'Sonda KRONOS');
    h_trail = animatedline('Color', [0 1 0], 'LineWidth', 2, 'HandleVisibility', 'off');
    h_text_vel = text(0.5, 0.07, '', 'Units', 'normalized', 'Color', 'w', 'FontSize', 12, 'FontName', 'Courier', 'FontWeight', 'bold','HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom'); 
    lgd = legend('Location', 'northeast');
    set(lgd, 'Color', 'k', 'TextColor', 'w', 'EdgeColor', 'w');

    % --- LOOP ---
    target_frames = 1000; 
    step_vis = max(1, floor(N_anim / target_frames));
    
    fprintf('Avvio animazione...\n');
    
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