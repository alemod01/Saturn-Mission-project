%% PLOT 2D: Traiettoria interplanetaria Earth SOI -> Earth SOI (proiezione XY, unità AU)
try
    
    sat_pos_au = state_cruise_earth_earth(:, 1:3); % Nx3 (AU)
    figure('Name', 'Traiettoria Earth->Earth (2D XY)'); clf; hold on
    
    % Sole
    plot(0,0,'y*','MarkerSize',12,'DisplayName','Sole');
    
    % Traiettoria della sonda (eliocentrica)
    plot(sat_pos_au(:,1), sat_pos_au(:,2), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Traiettoria sonda');

    %Traiettoria terra 
    step_m=0.1; 
    jd_vec = jd_start_0:step_m:jd_earth_arrival_actual;
    [~, r_earth_traj, v_earth_traj] = planet_orbit_coplanar(planets_elements.earth, jd_start_0, jd_earth_arrival_actual, jd_vec);
    plot(r_earth_traj(1, :), r_earth_traj(2, :), 'b--', 'LineWidth', 1.2, 'DisplayName', 'Orbita Terra');
    
    
    % Punti chiave
    plot(r_sat_earth_sp(1), r_sat_earth_sp(2), 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Uscita SOI Terra');
    plot(r_sat_interplanetary_earth_earth(1), r_sat_interplanetary_earth_earth(2), 'kx', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Sat end');
    plot(r_earth_sp(1), r_earth_sp(2), 'bo', 'MarkerFaceColor', 'b', 'DisplayName', 'Terra (alla partenza)');
    plot(r_earth_arr(1), r_earth_arr(2), 'co', 'MarkerFaceColor', 'c', 'DisplayName', 'Terra (arrivo)');
    
    % SOI (proiezione 2D) in AU
    theta = linspace(0,2*pi,240)';
    earth_soi_au = soi_earth / au;
    plot(r_earth_sp(1) + earth_soi_au*cos(theta), r_earth_sp(2) + earth_soi_au*sin(theta), 'b--', 'DisplayName', 'SOI Terra (proj)');
    plot(r_earth_arr(1) + earth_soi_au*cos(theta), r_earth_arr(2) + earth_soi_au*sin(theta), 'b--', 'DisplayName', 'SOI Terra (proj)');

    % Draw velocity direction arrows (scaled for visibility)
    try
        arrow_days = 5; % represent ~5 days displacement for visibility
        scale_seconds = 3600 * 24 * arrow_days;
        v_earth_vec = v_earth_sp * scale_seconds; % AU
        v_earth_arr_vec = v_earth_arr * scale_seconds;  % AU
        quiver(r_earth_sp(1), r_earth_sp(2), v_earth_vec(1), v_earth_vec(2), 0, 'b', 'LineWidth', 1.4, 'MaxHeadSize', 0.5);
        quiver(r_earth_arr(1), r_earth_arr(2), v_earth_arr_vec(1), v_earth_arr_vec(2), 0, 'c', 'LineWidth', 1.4, 'MaxHeadSize', 0.5);
        text(r_earth_sp(1) + v_earth_vec(1)*0.06, r_earth_sp(2) + v_earth_vec(2)*0.06, 'v_{Earth,start} (scaled)', 'Color', 'b');
        text(r_earth_arr(1) + v_earth_arr_vec(1)*0.06, r_earth_arr(2) + v_earth_arr_vec(2)*0.06, 'v_{Earth,arr} (scaled)', 'Color', 'c');
    catch
        % non-critico se quiver fallisce
    end

    xlabel('X [AU]'); ylabel('Y [AU]');
    title('Traiettoria interplanetaria: Earth SOI -> Earth SOI (XY projection)');
    axis equal; grid on; legend('Location','best');

    % Print human-readable dates for plotted planet positions and trajectory interval
    try
        fprintf('\nPlotted Earth position at (JD = %.6f): %s\n', jd_earth_sp, datestr(datetime(jd_earth_sp,'convertfrom','juliandate','TimeZone',timezone)));
    catch
        try
            fprintf('\nPlotted Earth position at (approx JD = %.6f): %s\n', t_vec_escape(end)/86400, datestr(datetime(t_vec_escape(end)/86400,'convertfrom','juliandate','TimeZone',timezone)));
        catch
        end
    end
    try
        fprintf('Plotted Earth arrival at (JD = %.6f): %s\n', jd_earth_arrival_actual, datestr(datetime(jd_earth_arrival_actual,'convertfrom','juliandate','TimeZone',timezone)));
    catch
    end
    try
        t0_jd = t_vec_cruise_earth_earth(1)/86400;
        t1_jd = t_vec_cruise_earth_earth(end)/86400;
        fprintf('Trajectory plotted from %s (JD=%.6f) to %s (JD=%.6f)\n', datestr(datetime(t0_1jd,'convertfrom','juliandate','TimeZone',timezone)), t0_jd, datestr(datetime(t1_jd,'convertfrom','juliandate','TimeZone',timezone)), t1_jd);
    catch
    end
catch ME
    warning(ME.noPlot2D , 'Impossibile creare il plot 2D Earth->Earth: %s', ME.message)
end
