function [r, lambda, sinphi, cosphi, s, rhoP,phi] = MT_ECEF2GEOCENTRIC(r_ECEF)
% Convierte r_ECEF -> r, phi, lambda y variables auxiliares.
%   r_ECEF: vector posicion 3x1 [m]
%   r, phi, lambda : vector posicion en coordenadas geocentricas [m, rad]
%   sinphi, cosphi : seno y coseno de phi
%   s               : sqrt(1 - sin(phi)^2) = abs(cosphi)
%   rhoP            : (a_ref/r)^(0:N+2) vector de 4x1 (se utiliza N=3)

% --- configuración fija 
N = 3;
a_ref = 6378136.3;  % [m] - radio referencia

x = r_ECEF(1); y = r_ECEF(2); z = r_ECEF(3);
r = sqrt(x*x + y*y + z*z);

% coordenadas geocéntricas (numéricamente estables)
rho = sqrt(x*x + y*y);          % proyección en XY
phi = atan2(z, rho);            % latitud geocéntrica
lambda = atan2(y, x);           % longitud

sinphi = sin(phi);
cosphi = cos(phi);

s = sqrt(max(0, 1 - sinphi^2));

% potencias radiales
rhoP = (a_ref / r) .^ (0:(N));
end
