function M_drag = MD_momento_drag(F_drag_B, CP)
% MD_MOMENTO_DRAG_FROM_FORCE
% Entradas:
%   F_drag_B      : 3x1 fuerza de drag en BODY [N]
%   CP            : 3x1 vector del CM al centro de presión en BODY [m]
% Salida:
%   M_drag           : 3x1 momento (torque) en BODY [N*m]

M_drag = (cross(CP, F_drag_B))';
end
