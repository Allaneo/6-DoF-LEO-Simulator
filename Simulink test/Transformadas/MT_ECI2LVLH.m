function v_lvlh = MT_ECI2LVLH(r_eci, v_eci, vec_in_eci)
% MT_ECI2LVLH
% Convierte un vector expresado en ECI a componentes en la terna LVLH.
% Definición LVLH:
%   z = -r/|r|
%   y = -(r x v)/|r x v|
%   x = y x z
%
% Inputs:
%   r_eci     [3x1] posición en ECI
%   v_eci     [3x1] velocidad en ECI
%   v_in_eci  [3x1] vector a transformar (en ECI)
%
% Output:
%   v_lvlh    [3x1] componentes del vector en LVLH

r = r_eci(:);
v = v_eci(:);
u = vec_in_eci(:);

z = -r / norm(r);

h = cross(r, v);
y = -h / norm(h);

x = cross(y, z);
x = x / norm(x);

% Proyección a la base LVLH (componentes)
v_lvlh = [dot(x,u); dot(y,u); dot(z,u)];
end
