function [value, isterminal, direction] = stopCondition_interplanetary(t, y,jd_start , planet_elements, R_soi)
    % STOPCONDITION Event function for ode45 to stop integration at SOI
    % border interplanetary 
    % 
    % Inputs:
    %   t                   - Time
    %   y                   - State vector
    %   jd_start            - jd start date ode integration 
    %   planet_elements     - planet 
    %   R_soi               - Radius of the Sphere of Influence (limit)
    
    % 1. Calcola la distanza attuale dal centro del corpo
    r = y(1:3);
    

    %2. Individua la posizione del pianeta all'istante considerato frame
    %[J200]
    jd_end=jd_start + t/60/60/24 ; 
    [~, r_mars_actual , ~] = planet_orbit_coplanar(planet_elements, jd_start, jd_end , [jd_start, jd_end]);
    vec_diff= abs(r_mars_actual(: , end) - r ); 

    % 2. Definisci l'evento (VALUE)
    % MATLAB cerca quando questa variabile attraversa lo Zero.
    % Se dist = R_soi, allora value = 0 -> EVENTO TROVATO!
    value =  norm (vec_diff) - R_soi;

    % 3. ISTERMINAL = 1 significa "Ferma l'integrazione quando trovi l'evento"
    isterminal = 1;

    % 4. DIRECTION = 1 significa "Ferma solo se la funzione sta crescendo"
    % (Cioè se stiamo USCENDO dalla sfera, non se stiamo entrando)
    direction = 1;
end