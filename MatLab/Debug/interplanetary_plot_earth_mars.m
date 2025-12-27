%% PLOT 2D: Traiettoria interplanetaria Earth SOI -> Mars SOI (proiezione XY, unità AU)
try
    
    sat_pos_au = state_cruise_earth_mars(:, 1:3); % Nx3 (AU)
    figure('Name', 'Traiettoria Earth->Mars (2D XY)'); clf; hold on
    
    % Sole
    plot(0,0,'y*','MarkerSize',12,'DisplayName','Sole');
    
    % Traiettoria della sonda (eliocentrica)
    plot(sat_pos_au(:,1), sat_pos_au(:,2), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Traiettoria sonda');

    %Traiettoria marte 
    step_m=0.1; 
    jd_vec = jd_start:step_m:jd_mars_arrival_actual;
    [~, r_mars_traj, v_mars_traj] = planet_orbit_coplanar(planets_elements.mars, jd_start, jd_mars_arrival_actual, jd_vec);
    plot(r_mars_traj(1, :), r_mars_traj(2, :), 'm--', 'LineWidth', 1.2, 'DisplayName', 'Orbita Marte');
    
    
    % Punti chiave
    plot(r_sat_earthfb_sp(1), r_sat_earthfb_sp(2), 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Uscita SOI Terra');
    plot(r_sat_interplanetary_earth_mars(1), r_sat_interplanetary_earth_mars(2), 'kx', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Sat end');
    plot(r_earthfb_sp(1), r_earthfb_sp(2), 'bo', 'MarkerFaceColor', 'b', 'DisplayName', 'Terra (alla partenza)');
    %plot(r_mars_traj(1, 1), r_mars_traj(2, 1), 'mo', 'MarkerFaceColor', 'w', 'DisplayName', 'Marte (Start)');
    plot(r_mars_arr(1), r_mars_arr(2), 'mo', 'MarkerFaceColor', 'm', 'DisplayName', 'Marte (arrivo)');
    
    % SOI (proiezione 2D) in AU
    theta = linspace(0,2*pi,240)';
    earth_soi_au = soi_earth / au;
    mars_soi_au = soi_mars / au;
    plot(r_earthfb_sp(1) + earth_soi_au*cos(theta), r_earthfb_sp(2) + earth_soi_au*sin(theta), 'b--', 'DisplayName', 'SOI Terra (proj)');
    plot(r_mars_arr(1) + mars_soi_au*cos(theta), r_mars_arr(2) + mars_soi_au*sin(theta), 'm--', 'DisplayName', 'SOI Marte (proj)');

    % Draw velocity direction arrows (scaled for visibility)
    try
        arrow_days = 5; % represent ~5 days displacement for visibility
        scale_seconds = 3600 * 24 * arrow_days;
        v_earth_vec = v_earthfb_sp * scale_seconds; % AU
        v_mars_vec = v_mars_arr * scale_seconds;  % AU
        quiver(r_earthfb_sp(1), r_earthfb_sp(2), v_earth_vec(1), v_earth_vec(2), 0, 'b', 'LineWidth', 1.4, 'MaxHeadSize', 0.5);
        quiver(r_mars_arr(1), r_mars_arr(2), v_mars_vec(1), v_mars_vec(2), 0, 'm', 'LineWidth', 1.4, 'MaxHeadSize', 0.5);
        text(r_earthfb_sp(1) + v_earth_vec(1)*0.06, r_earthfb_sp(2) + v_earth_vec(2)*0.06, 'v_{Earth} (scaled)', 'Color', 'b');
        text(r_mars_arr(1) + v_mars_vec(1)*0.06, r_mars_arr(2) + v_mars_vec(2)*0.06, 'v_{Mars} (scaled)', 'Color', 'm');
    catch
        % non-critico se quiver fallisce
    end

    xlabel('X [AU]'); ylabel('Y [AU]');
    title('Traiettoria interplanetaria: Earth SOI -> Mars SOI (XY projection)');
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
        fprintf('Plotted Mars position at arrival (JD = %.6f): %s\n', jd_mars_arrival_actual, datestr(datetime(jd_mars_arrival_actual,'convertfrom','juliandate','TimeZone',timezone)));
    catch
    end
    try
        t0_jd = t_vec_cruise_earth_mars(1)/86400;
        t1_jd = t_vec_cruise_earth_mars(end)/86400;
        fprintf('Trajectory plotted from %s (JD=%.6f) to %s (JD=%.6f)\n', datestr(datetime(t0_jd,'convertfrom','juliandate','TimeZone',timezone)), t0_jd, datestr(datetime(t1_jd,'convertfrom','juliandate','TimeZone',timezone)), t1_jd);
    catch
    end
catch ME
    warning(ME.noPlot2D , 'Impossibile creare il plot 2D Earth->Mars: %s', ME.message)
end