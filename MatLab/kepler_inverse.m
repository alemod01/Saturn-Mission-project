function E_deg = kepler_inverse(e, M_deg, E0, tol)
% This function solves the inverse kepler problem, given the orbital
% elements and the mean anomaly in radiants.

M_rad = deg2rad(M_deg);

% Elliptic case:
if e < 1
    % Compute the initial guess E0
    E0_vec = M_rad + e*sin(M_rad);
    E_n = E0_vec;
    dE = 1;
    
    % Iterate till tolerance condition is satisfied
    while 1 - prod(abs(dE) <= tol)
        % Modify the Kepler equation to handle non-zero E0
        dM = M_rad - (E_n - e * sin(E_n)) + (E0 - e * sin(E0));
        dE = dM / (1 - e * cos(E_n));

        % Update E for the next iteration
        E_n = E_n + dE;  
    end
    E_rad = E_n;

    E_deg = rad2deg(E_rad);

% Hyperbolic case:
elseif e > 1
        % Initial guess for hyperbolic case using F0 (hyperbolic eccentric anomaly)
    F0_vec = M_rad + e * sinh(M_rad);
    F_n = F0_vec;
    dF = 1;
    
    % Iterate till tolerance condition is satisfied
    while 1 - prod(abs(dF) <= tol)
        % Modify the hyperbolic Kepler equation to handle non-zero F0
        dM = M_rad - (e * sinh(F_n) - F_n) + (e * sinh(E0) - E0);
        dF = dM / (e * cosh(F_n) - 1);
        F_n = F_n + dF;  % Update F for next iteration
    end
    E_rad = F_n; 

    E_deg = rad2deg(E_rad);
end
