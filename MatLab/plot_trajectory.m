function plot_trajectory(planets_elements , jd_start,jd_end , step)
    % Compute planets states from start to end date with a step of a day or
    % some hours

    jd_vec = jd_start:step:jd_end;
    [~, r, ~] = planet_orbit_coplanar(planets_elements, jd_start, jd_end, jd_vec );

    figure('Name', 'Displacement Check', 'Color', 'w');
    hold on; grid on; axis equal;
    xlabel('X [AU]'); ylabel('Y [AU]');
    title({'Traiettoria '});

    % 1. Disegna il Sole (Riferimento centrale)
    plot(0, 0, 'y.', 'MarkerSize', 40, 'DisplayName', 'Sole');

    % Disegna la scia di pallini ("pallini" ogni step ora)
    plot(r(1, :), r(2, :), 'r.', 'MarkerSize', 6, 'DisplayName', 'Posizione (ogni 24h)');

    % 3. Evidenzia Inizio e Fine per capire il verso
    plot(r(1, 1), r(2, 1), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 8, 'DisplayName', 'Start');
    plot(r(1, end), r(2, end), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8, 'DisplayName', 'Arrivo');

    legend('Location', 'best');

    % Aggiungi una freccia o testo per chiarire la direzione (Opzionale)
    text(r(1,1), r(2,1)+0.1, '  Start', 'Color', 'b');
    text(r(1,end), r(2,end)+0.1, '  End', 'Color', 'k');
    end