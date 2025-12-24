%% ANIMAZIONE 2D: Sonda + Pianeti Sincronizzati
try
    close all;

    % --- 0. IMPOSTAZIONI VIDEO ---
    save_video = false;  % <--- True per salvare video 
    video_name = 'Kronos_Mission.mp4';
    video_skip = 40;

    % --- 1. PREPARAZIONE DATI SINCRO ---
    k_step = 100; %per la parte denro la SoI 
    t_vec_saturn_soi_red = t_vec_saturn_soi(1:k_step:end);

    full_time_seconds = [t_vec_cruise_earth_earth; t_vec_cruise_earth_mars; t_vec_cruise_mars_saturn; t_vec_saturn_soi_red];
    N = size(t_vec_saturn_soi , 1);
    
    % Converto tutto in Date Giuliane (JD) assolute
    jd_start_sim=t_vec_earth_escape(end)/60/60/24; 
    
    full_jd_vector = (full_time_seconds / 86400)';
    
    % --- 2. CALCOLO POSIZIONE PIANETI (SUI TEMPI DELLA SONDA) ---
    % Nota: passiamo 'full_jd_vector' come ultimo argomento invece di un vettore creato a mano.
    % La funzione restituirà vettori della stessa lunghezza della sonda (es. 91421 punti).
    
    [~, r_earth_sync, ~]  = planet_orbit_coplanar(planets_elements.earth,  0, 0, full_jd_vector);
    [~, r_mars_sync, ~]   = planet_orbit_coplanar(planets_elements.mars,   0, 0, full_jd_vector);
    [~, r_saturn_sync, v_saturn_sync] = planet_orbit_coplanar(planets_elements.saturn, 0, 0, full_jd_vector);
    
    % --- 3. COSTRUZIONE TRAIETTORIA E VELOCITÀ SONDA ---

    %fase interplanetaria 
    % Posizione (AU)
    traj_part1 = state_cruise_earth_earth(:, 1:3);
    traj_part2 = state_cruise_earth_mars(:, 1:3);
    traj_part3 = state_cruise_mars_saturn(:, 1:3);

    % Velocità (AU/s) -> Cocertita in km/s
    vel_part1_km = state_cruise_earth_earth(:, 4:6) * au;
    vel_part2_km = state_cruise_earth_mars(:, 4:6) * au;
    vel_part3_km = state_cruise_mars_saturn(:, 4:6) * au;

    % Fase SOI Saturno 
    state_sat_saturn_soi_red = state_sat_saturn_soi(1:k_step:end, :);
    N_soi = size(t_vec_saturn_soi_red, 1);

    %estrazione posizone velocità Saturno 
    r_saturn_soi_portion = r_saturn_sync(:, end - N_soi + 1 : end)'; % (Nx3)
    v_saturn_soi_portion = v_saturn_sync(:, end - N_soi + 1 : end)'; % (Nx3) AU/s

    %posizone sat 
    state_sat_saturn_soi_au = state_sat_saturn_soi_red(:, 1:3)/au + r_saturn_soi_portion;
    %velocità eliocentrica 
    v_rel_km = state_sat_saturn_soi_red(:, 4:6);
    v_saturn_planet_km = v_saturn_soi_portion * au;
    vel_part4_km = v_rel_km + v_saturn_planet_km;

    % Concatenazione totale 
    full_trajectory_au = [traj_part1; traj_part2; traj_part3; state_sat_saturn_soi_au];
    full_velocity_km   = [vel_part1_km; vel_part2_km; vel_part3_km; vel_part4_km];
    
    % Modula velocità per ogni istante 
    full_velocity_mag = sqrt(sum(full_velocity_km.^2, 2));


    % (Opzionale) Calcolo orbite statiche per sfondo 
    jd_static = linspace(min(full_jd_vector), max(full_jd_vector), 500);
    [~, r_earth_bg, ~] = planet_orbit_coplanar(planets_elements.earth, 0, 0, jd_static);
    [~, r_mars_bg, ~]  = planet_orbit_coplanar(planets_elements.mars, 0, 0, jd_static);
    [~, r_saturn_bg, ~]= planet_orbit_coplanar(planets_elements.saturn, 0, 0, jd_static);

    % --- 4. SETUP GRAFICO ---
    figure('Name', 'Interplanetary Trajectory Analysis: Earth - Mars - Saturn' , 'Color', 'k', 'Position', [100, 100, 1200, 600]); 
    clf; hold on;  axis off;
    % % Impostiamo assi e testo bianchi su sfondo nero
    % set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'GridColor', 'w', 'GridAlpha', 0.3);
    % xlabel('X [AU]', 'Color', 'w'); 
    % ylabel('Y [AU]', 'Color', 'w');
    title('KRONOS Mission', 'Color', 'w');
    xlim([-8, 2]); ylim([-11, 2]);axis manual; 


    % --- ELEMENTI STATICI (Sfondo) ---
    plot(0,0,'y*','MarkerSize',14,'DisplayName','Sole');
    
    % Orbite pianeti (linee sottili tratteggiate)
    plot(r_earth_bg(1,:), r_earth_bg(2,:), '--', 'Color', [0.6 0.6 1], 'LineWidth',0.5, 'HandleVisibility','off');
    plot(r_mars_bg(1,:),  r_mars_bg(2,:),  '--', 'Color', [1 0.6 0.6], 'LineWidth',0.5, 'HandleVisibility','off');
    plot(r_saturn_bg(1,:),r_saturn_bg(2,:),'--', 'Color', [1 0.8 0.4], 'LineWidth',0.5, 'HandleVisibility','off');
    
    % Scia Sonda (tutto il percorso in grigio chiaro)
    plot(full_trajectory_au(:,1), full_trajectory_au(:,2), '-', 'Color', [0.8 0.9 0.8], 'LineWidth', 1, 'DisplayName', 'Percorso');

    % --- ELEMENTI DINAMICI ---
    
    % 1. Terra (Blu)
    h_earth = plot(r_earth_sync(1,1), r_earth_sync(2,1), 'o','MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none', 'MarkerSize', 6, 'DisplayName', 'Terra');
    % 2. Marte (Rosso)
    h_mars = plot(r_mars_sync(1,1), r_mars_sync(2,1), 'o','MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none', 'MarkerSize', 6, 'DisplayName', 'Marte');  
    % 3. Saturno (Arancione/Oro)
    h_saturn = plot(r_saturn_sync(1,1), r_saturn_sync(2,1), 'o','MarkerFaceColor', [1 0.6 0], 'MarkerEdgeColor', 'none', 'MarkerSize', 8, 'DisplayName', 'Saturno');
    % 4. Sonda (Rosso acceso con bordo nero)
    h_sat = plot(full_trajectory_au(1,1), full_trajectory_au(1,2), 'o','MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'MarkerSize', 5, 'DisplayName', 'Sonda');       
    % Scia sonda
    h_trail = animatedline('Color', [0 0.8 0], 'LineWidth', 1.5, 'HandleVisibility', 'off');
    
    % --- DATA ---
    % Posizioniamo a 0.05 (5% sinistra) e 0.95 (95% altezza).
    h_text_date = text(0.05, 0.95, 'Date: ...', ...
        'Units', 'normalized', ...
        'Color', 'w', ...
        'FontSize', 10, ...         
        'FontName', 'Courier', ...
        'FontWeight', 'bold', ...
        'VerticalAlignment', 'top'); 

    % --- VELOCITÀ  ---
    h_text_vel = text(0.5, 0, 'Vel: 0.00 km/s', ...
        'Units', 'normalized', ...
        'Color', 'w', ...
        'FontSize', 12, ...
        'FontName', 'Courier', ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center'); 

    % --- PRE-CALCOLO STRINGHE EVENTI FLY-BY ---
    str_earth_fb = string(datetime(jd_earth_arrival_actual, 'ConvertFrom', 'juliandate'), 'yyyy-MM-dd');
    str_mars_fb  = string(datetime(jd_mars_arrival_actual,  'ConvertFrom', 'juliandate'), 'yyyy-MM-dd');
    str_sat_arr  = string(datetime(jd_saturn_arrival_actual,'ConvertFrom', 'juliandate'), 'yyyy-MM-dd');

    % Legenda
    lgd = legend([h_sat, h_earth, h_mars, h_saturn], {'Sonda','Terra','Marte','Saturno'});
    set(lgd, 'Position', [0.87, 0.85, 0.12, 0.10]);
    set(lgd, 'Color', 'k', 'TextColor', 'w', 'EdgeColor', 'w');
    
    % --- APERTURA VIDEO ---
    if save_video
        v = VideoWriter(video_name, 'MPEG-4');
        v.FrameRate = 30; % 30 FPS standard
        v.Quality = 90;
        open(v);
        fprintf('Registrazione video avviata (Salvataggio 1 frame ogni %d cicli)...\n', video_skip);
    end

    
    % --- 5. CICLO DI ANIMAZIONE ---
    % Step animazione: aumentalo se va troppo lento (es. 20 o 50)
    step_anim = 5; 
    N = size(full_trajectory_au, 1);
    
    loop_counter = 0;
    fprintf('Avvio animazione (%d frames totali)...\n', floor(N/step_anim));
    
    for i = 1:step_anim:N
        loop_counter = loop_counter + 1;

        % Posizioni attuali (Sonda)
        sat_x = full_trajectory_au(i, 1);
        sat_y = full_trajectory_au(i, 2);
        
        % Posizioni attuali (Pianeti - sincronizzati dall'indice i)
        earth_x = r_earth_sync(1, i); earth_y = r_earth_sync(2, i);
        mars_x  = r_mars_sync(1, i);  mars_y  = r_mars_sync(2, i);
        sat_x_p = r_saturn_sync(1, i); sat_y_p = r_saturn_sync(2, i);
        
        % UPDATE GRAFICO
        set(h_sat,    'XData', sat_x,   'YData', sat_y);
        set(h_earth,  'XData', earth_x, 'YData', earth_y);
        set(h_mars,   'XData', mars_x,  'YData', mars_y);
        set(h_saturn, 'XData', sat_x_p, 'YData', sat_y_p);

        addpoints(h_trail, sat_x, sat_y);
        
        %UPDATE DATA & EVENTI
        current_jd = full_jd_vector(i);
        curr_date_str = string(datetime(current_jd, 'ConvertFrom', 'juliandate'), 'yyyy-MM-dd');
        
        % Data Corrente
        msg_lines = {['Date: ' char(curr_date_str)]}; 
        % Check fly-by Terra
        if current_jd >= jd_earth_arrival_actual
            msg_lines{end+1, 1} = ['Earth Fly-by: ' char(str_earth_fb)];
        end
        % Check fly-by Marte
        if current_jd >= jd_mars_arrival_actual
            msg_lines{end+1, 1} = ['Mars Fly-by:  ' char(str_mars_fb)];
        end
        % Check fly-by Saturno
        if current_jd >= jd_saturn_arrival_actual
             msg_lines{end+1, 1} = ['Saturn SOI:   ' char(str_sat_arr)];
        end
        % Aggiornamento Testo
        set(h_text_date, 'String', msg_lines);

        % --- UPDATE VELOCITÀ ---
        v_current = full_velocity_mag(i);
        set(h_text_vel, 'String', sprintf('Heliocentric Velocity: %.2f km/s', v_current));

        drawnow limitrate;

        % Video 
        if save_video && mod(loop_counter, video_skip) == 0
            drawnow; 
            writeVideo(v, getframe(gcf));
        end

    end
    
   % Frame finale
    set(h_sat, 'XData', full_trajectory_au(end, 1), 'YData', full_trajectory_au(end, 2));
    drawnow;
    if save_video
        writeVideo(v, getframe(gcf));
        close(v);
        fprintf('Video salvato in: %s\n', pwd);
    end
    fprintf('Animazione completata.\n');


catch ME
    if exist('v', 'var')
        try close(v); catch; end % Proviamo a chiudere, se fallisce amen
    end
    warning('Errore animazione: %s', ME.message); 
end