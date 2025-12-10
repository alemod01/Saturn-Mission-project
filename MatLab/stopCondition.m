function [value, isterminal, direction] = stopCondition(t, y, R_soi)
    % STOPCONDITION Event function for ode45 to stop integration at SOI border
    %
    % Inputs:
    %   t     - Time
    %   y     - State vector
    %   R_soi - Radius of the Sphere of Influence (limit)
    
    % 1. Calcola la distanza attuale dal centro del corpo
    r = y(1:3);
    current_dist = norm(r);

    % 2. Definisci l'evento (VALUE)
    % MATLAB cerca quando questa variabile attraversa lo Zero.
    % Se dist = R_soi, allora value = 0 -> EVENTO TROVATO!
    value = current_dist - R_soi;

    % 3. ISTERMINAL = 1 significa "Ferma l'integrazione quando trovi l'evento"
    isterminal = 1;

    % 4. DIRECTION = 1 significa "Ferma solo se la funzione sta crescendo"
    % (Cioè se stiamo USCENDO dalla sfera, non se stiamo entrando)
    direction = 1;
end