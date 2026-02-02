function F_drag_I = MD_drag(v_rel_I, Aproj, rho, CD)
% --- MD_DRAG ---
% Calcula fuerza aerodinámica utilizando la velocidad relativa considerando atmosfera corrotante, 
% el área proyectada en la dirección de la velocidad, la densidad atmosférica y el coeficiente de arrastre.

% Inputs:
%  r_I, v_I         : 3x1 ECI [m],[m/s]
%  q                : 4x1 quaternion scalar-first (inertial->body)
%  Cd               : coeficiente de arrastre (scalar)
%  rho              : densidad atmosférica en punto (scalar) [kg/m^3]
%  A_proj           : área proyectada en la dirección de la velocidad
% Outputs:
%  F_B              : 3x1 fuerza drag en BODY [N]

% 1) velocidad relativa en ECI (atmósfera co-rotante)
v = norm(v_rel_I);
vhat_I = v_rel_I / v;

% 2) magnitud y vector de la fuerza en ECI
F_mag = 0.5 * rho * CD * Aproj * v^2;
F_drag_I = - F_mag * vhat_I;    % fuerza total en ECI (componentes en ECI)

end
