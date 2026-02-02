function M_srp_B = MD_momento_srp(F_srp_B, CP_srp_B)
% MD_momento_srp
% Calcula el torque por presión solar (SRP) en la terna Body. Reemplazar
% por bloque cross en simulink?
%
% Entradas:
%   F_srp_B   : 3x1 fuerza SRP en Body [N]
%   CP_srp_B  : 3x1 vector desde el origen Body (CG) al CP equivalente SRP [m]
%
% Salida:
%   M_srp_B   : 3x1 torque SRP en Body [N*m]

M_srp_B = cross(CP_srp_B, F_srp_B);

end
