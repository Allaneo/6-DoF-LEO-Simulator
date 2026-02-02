function M_drag = MD_momento_drag(F_drag_B, CP,l_minimo_drag)
% --- MD_momento_drag ---
% Calcula el torque aerodinámico por drag en la terna Body a partir de la 
% fuerza de drag y el brazo de palanca correspondiente al centro de presión
% equivalente. El brazo se proyecta perpendicularmente a la dirección de la
% fuerza para obtener el componente efectivo que genera momento, y se 
% impone un brazo mínimo l_minimo_drag cuando existe componente perpendicular
% pequeña, con el objetivo de evitar momentos nulos o demasiado bajos como
% "buena práctica".

% Entradas:
%   F_drag_B      : 3x1 fuerza de drag en BODY [N]
%   CP            : 3x1 vector del CM al centro de presión en BODY [m]
% l_minimo_drag: valor escalar mínimo para el brazo de palanca [m]

% Salida:
%   M_drag           : 3x1 momento (torque) en BODY [N*m]

Fn = norm(F_drag_B);

% 1. PROTECCIÓN INMEDIATA: Antes de dividir, chequeamos si hay fuerza.
if Fn == 0
    M_drag = [0;0;0];
    return
end

Fhat = F_drag_B / Fn;

% Cálculo del brazo perpendicular
CP_perp = CP - (CP.' * Fhat) * Fhat;
rp = norm(CP_perp);

% 2. Lógica del brazo de palanca
if rp>0 && rp < l_minimo_drag
    % Si existe un desalineamiento pero es muy pequeño, lo amplificamos
    % para cumplir el requisito de "peor caso" o robustez.
    CP_perp = CP_perp .* (l_minimo_drag / rp);
elseif rp ==0
    CP_perp = [l_minimo_drag;0;0];
end

% 3. Cálculo final
M_drag = cross(CP_perp, F_drag_B);

end