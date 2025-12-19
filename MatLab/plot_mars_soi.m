function plot_mars_soi(t_vec, r_sat, soi_mars, v_mars, R_planet)
% PLOT_MARS_SOI Plot della propagazione del satellite dentro la SOI di Marte
%   plot_mars_soi(t_vec, r_sat, soi_mars, v_mars)
%   - t_vec: vettore tempi (s)
%   - r_sat: posizione del satellite (3xN o Nx3) in km, relativo al centro di Marte
%   - soi_mars: raggio SOI di Marte in km (opzionale, default 577400 km)
%   - v_mars: vettore velocita' di Marte in km/s (3x1) opzionale; se fornito
%             viene disegnata una freccia che indica la direzione del moto
%
% Esempio d'uso (da `main_rosetta.m` dopo la propagazione):
%   plot_mars_soi(t_vec_mars_escape, r_sat_mars_escape, soi_mars, v_mars_sp)

if nargin < 3 || isempty(soi_mars)
    soi_mars = 577400; % km
end

% Normalizza forma r_sat a 3xN
if size(r_sat,1) == 3
    R = r_sat;
else
    R = r_sat';
end

% 3D plot: Marte (solido), SOI (wireframe), traiettoria
fig1 = figure('Name','Traiettoria 3D dentro SOI Marte');
hold on; grid on; axis equal;
% Marte
[xs, ys, zs] = sphere(48);
surf(R_planet*xs, R_planet*ys, R_planet*zs, 'FaceColor',[1 0.6 0.4],'EdgeColor','none','FaceAlpha',0.9);
% SOI wireframe
[XS,YS,ZS] = sphere(36);
mesh(soi_mars*XS, soi_mars*YS, soi_mars*ZS, 'EdgeColor',[0.5 0.5 0.5],'FaceAlpha',0,'LineStyle','--');

% Traiettoria
plot3(R(1,:), R(2,:), R(3,:), '-b', 'LineWidth', 1.5);
plot3(R(1,1), R(2,1), R(3,1), 'og', 'MarkerFaceColor','g', 'DisplayName','Ingresso');
plot3(R(1,end), R(2,end), R(3,end), 'or', 'MarkerFaceColor','r', 'DisplayName','Uscita');

xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
title('Traiettoria sonda dentro SOI di Marte');
legend({'Marte','SOI','Traiettoria','Ingresso','Uscita'}, 'Location','best');
view(45,25);
% Plot della direzione di moto di Marte se fornita
if nargin >= 4 && ~isempty(v_mars)
    v = v_mars(:);
    if norm(v) > 0
        v_dir = v / norm(v);
        arrow_len = soi_mars * 0.18; % lunghezza della freccia per visibilita'
        arrow = v_dir * arrow_len;
        quiver3(0,0,0, arrow(1), arrow(2), arrow(3), 0, 'k', 'LineWidth', 2, 'MaxHeadSize', 1);
        text(arrow(1)*1.05, arrow(2)*1.05, arrow(3)*1.05, sprintf('v_{Mars}=%.3f km/s', norm(v)), 'FontSize', 9);
    end
end

% 2D XY projection
fig2 = figure('Name','Proiezione XY della Traiettoria');
plot(R(1,:), R(2,:), '-b', 'LineWidth', 1.2); hold on; grid on; axis equal;
theta = linspace(0,2*pi,720);
plot(R_planet*cos(theta), R_planet*sin(theta), 'Color',[1 0.6 0.4], 'LineWidth', 1.5);
plot(soi_mars*cos(theta), soi_mars*sin(theta), '--', 'Color',[0.5 0.5 0.5]);
plot(R(1,1), R(2,1), 'og', 'MarkerFaceColor','g');
plot(R(1,end), R(2,end), 'or', 'MarkerFaceColor','r');
xlabel('X (km)'); ylabel('Y (km)');
title('Proiezione XY della traiettoria dentro SOI di Marte');
legend('Traiettoria','Marte (r)','SOI','Ingresso','Uscita','Location','best');
% 2D arrow for Mars velocity
if nargin >= 4 && ~isempty(v_mars)
    v = v_mars(:);
    if norm(v) > 0
        v_dir = v / norm(v);
        arrow_len = soi_mars * 0.18;
        arrow2 = v_dir(1:2) * arrow_len;
        quiver(0,0, arrow2(1), arrow2(2), 0, 'k', 'LineWidth', 2, 'MaxHeadSize', 1);
    end
end
end
