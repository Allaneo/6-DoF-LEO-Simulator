function Aproj = MT_area_proyectada(v_rel_I, q, A_transversal)
% MD_AREA_PROJECTED  Calcula area proyectada (satélite modelado como paralelepípedo/cubo)
% Inputs:
%   v_rel_I : 3x1 relativa en ECI (m/s)
%   q       : 4x1 quaternion scalar-first
%   A       : 3x1 área transversal

% Output:
%   Aproj   : scalar projected area [m^2]

% relative velocity in BODY
v_rel_B = MT_ECI2BODY(q,v_rel_I);
v_rel = norm(v_rel_B);
vhat = v_rel_B / v_rel;

% projected area (orthographic approximation)
Aproj = A_transversal(1)*abs(vhat(1)) + A_transversal(2)*abs(vhat(2)) + A_transversal(3)*abs(vhat(3));
end
