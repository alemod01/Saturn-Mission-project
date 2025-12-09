function dydt = satellite_ode(t, y, mu)
    % SATELLITE_ODE ODE system for the 2-body problem
    %
    % Inputs:
    %   t  - Time (not used in 2-body dynamics, but required by ode45)
    %   y  - State vector [rx; ry; rz; vx; vy; vz] (6x1)
    %   mu - Gravitational parameter of the central body
    %
    % Output:
    %   dydt - Derivative of the state [vx; vy; vz; ax; ay; az] (6x1)

    % 1.position (r) e velocity (v) from y
    r = y(1:3); % Primi 3 elementi: Posizione
    v = y(4:6); % Ultimi 3 elementi: Velocità

    % 2. Calcola la distanza dal centro (norma del raggio)
    r_norm = norm(r);

    % 3. Calcola l'accelerazione dovuta alla gravità (Legge di Newton)
    % a = -mu / r^3 * r_vector
    a = -mu / (r_norm^3) * r;

    % 4. Costruisci il vettore derivata (dydt)
    % La derivata della posizione è la velocità (v)
    % La derivata della velocità è l'accelerazione (a)
    dydt = [v; a]; 
end