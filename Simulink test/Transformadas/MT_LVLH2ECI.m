function vec_eci = MT_LVLH2ECI(r_eci, v_eci, vec_in_lvlh)
% MT_LVLH2ECI
% Convierte componentes LVLH a vector expresado en ECI.
% Definición LVLH:
%   z = -r/|r|
%   y = -(r x v)/|r x v|
%   x = y x z
%
% Inputs:
%   r_eci       [3x1] posición en ECI
%   v_eci_orb   [3x1] velocidad orbital en ECI (para definir LVLH)
%   v_in_lvlh   [3x1] componentes en LVLH
%
% Output:
%   v_eci       [3x1] vector en ECI

r = r_eci(:);
v = v_eci(:);
uL = vec_in_lvlh(:);

z = -r / norm(r);

h = cross(r, v);
y = -h / norm(h);

x = cross(y, z);
x = x / norm(x);

C = [x y z];

% Reconstrucción del vector en ECI
vec_eci = C*vec_in_lvlh;
end
