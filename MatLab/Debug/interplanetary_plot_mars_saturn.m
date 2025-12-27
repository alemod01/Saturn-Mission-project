%% PLOT 2D: Traiettoria interplanetaria Mars SOI -> Saturn SoI (proiezione XY, unità AU)
try
    
    sat_pos_au = state_cruise_mars_saturn(:, 1:3); % Nx3 (AU)
    figure('Name', 'Traiettoria Mars -> Saturn (2D XY)'); clf; hold on 
    
    % Sole
    plot(0,0,'y*','MarkerSize',12,'DisplayName','Sole');
    
    
    % Traiettoria della sonda (eliocentrica)
    plot(sat_pos_au(:,1), sat_pos_au(:,2), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Traiettoria Sonda');

    %Traiettoria marte 
    step_m=0.1; 
    jd_vec = jd_start:step_m:jd_saturn_arrival_actual;
    [~, r_saturn_traj, v_saturn_traj] = planet_orbit_coplanar(planets_elements.saturn, jd_start, jd_saturn_arrival_actual, jd_vec);
    plot(r_saturn_traj(1, :), r_saturn_traj(2, :), '--', 'Color', [1, 0.5, 0],'LineWidth', 1.2, 'DisplayName', 'Orbita Saturno');
    
    
    % Punti chiave
    plot(r_sat_marsfb_sp(1), r_sat_marsfb_sp(2), 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Uscita SOI Marte');
    text(r_sat_marsfb_sp(1) + r_sat_marsfb_sp(1)*0.01, r_sat_marsfb_sp(2) + r_sat_marsfb_sp(2)*0.01, 'Uscita SoI Marte', 'Color', 'r');
    plot(r_sat_interplanetary_mars_saturn(1), r_sat_interplanetary_mars_saturn(2), 'kx', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Sat end');
    plot(r_mars_sp(1), r_mars_sp(2), 'mo', 'MarkerFaceColor', 'm', 'DisplayName', 'Marte (alla partenza)');
    plot(r_saturn_arr(1), r_saturn_arr(2), 'o', 'Color', [1, 0.5, 0], 'MarkerFaceColor', [1, 0.5, 0],'DisplayName', 'Saturno (arrivo)');
    
    % SOI (proiezione 2D) in AU
    theta = linspace(0,2*pi,240)';
    saturn_soi_au = soi_saturn / au;
    mars_soi_au = soi_mars / au;
    plot(r_mars_sp(1) + mars_soi_au*cos(theta), r_mars_sp(2) + mars_soi_au*sin(theta), 'm--', 'DisplayName', 'SOI Marte (proj)');
    plot(r_saturn_arr(1) + saturn_soi_au*cos(theta), r_saturn_arr(2) + saturn_soi_au*sin(theta), '--', 'Color', [1, 0.5, 0], 'DisplayName', 'SOI Saturno (proj)');


    % Draw velocity direction arrows (scaled for visibility)
    try
        arrow_days = 5; % represent ~5 days displacement for visibility
        scale_seconds = 3600 * 24 * arrow_days;
        v_mars_vec = v_mars_sp * scale_seconds;  % AU
        v_saturn_vec = v_saturn_arr * scale_seconds; % AU
        quiver(r_mars_sp(1), r_mars_sp(2), v_mars_vec(1), v_mars_vec(2), 0,'m', 'LineWidth', 3.4, 'MaxHeadSize', 0.5);
        quiver(r_saturn_arr(1), r_saturn_arr(2), v_saturn_vec(1), v_saturn_vec(2), 0,  '-', 'Color', [1, 0.5, 0] , 'LineWidth', 1.4, 'MaxHeadSize', 0.5);
        text(r_mars_sp(1) + v_mars_vec(1)*1.06, r_mars_sp(2) + v_mars_vec(2)*1.06, 'v_{Mars}', 'Color', 'm');
        text(r_saturn_arr(1) + v_saturn_vec(1)*1.06, r_saturn_arr(2) + v_saturn_vec(2)*1.06, 'v_{Saturn}' , 'Color', [1, 0.5, 0]);
    catch
        % non-critico se quiver fallisce
    end

    xlabel('X [AU]'); ylabel('Y [AU]');
    title('Traiettoria interplanetaria: Mars SOI -> Saturn SOI (XY projection)');
    axis equal; grid on; 
    %legend('Sole' ,'Traiettoria Sonda' , 'Orbita Saturno' ,  'Uscita SOI Marte' ,'Ingresso SoI Saturno' );
    %legend('Location','best');
    legend('Sole','Traiettoria Sonda' , 'Orbita Saturno' , 'Uscita SOI Marte' , 'Sat end'  )

    % Print human-readable dates for plotted planet positions and trajectory interval
    try
        fprintf('\nPlotted Mars position at (JD = %.6f): %s\n', jd_mars_sp, datestr(datetime(jd_mars_sp,'convertfrom','juliandate','TimeZone',timezone)));
    catch
        try
            fprintf('\nPlotted Mars position at (approx JD = %.6f): %s\n', t_vec_mars_escape(end)/86400, datestr(datetime(t_vec_mars_escape(end)/86400,'convertfrom','juliandate','TimeZone',timezone)));
        catch
        end
    end
    try
        fprintf('Plotted Saturn position at arrival (JD = %.6f): %s\n', jd_saturn_arrival_actual, datestr(datetime(jd_saturn_arrival_actual,'convertfrom','juliandate','TimeZone',timezone)));
    catch
    end
    try
        t0_jd = t_vec_cruise_mars_saturn(1)/86400;
        t1_jd = t_vec_cruise_mars_saturn(end)/86400;
        fprintf('Trajectory plotted from %s (JD=%.6f) to %s (JD=%.6f)\n', datestr(datetime(t0_jd,'convertfrom','juliandate','TimeZone',timezone)), t0_jd, datestr(datetime(t1_jd,'convertfrom','juliandate','TimeZone',timezone)), t1_jd);
    catch
    end
catch ME
    warning(ME.noPlot2D , 'Impossibile creare il plot 2D Earth->Mars: %s', ME.message)
end