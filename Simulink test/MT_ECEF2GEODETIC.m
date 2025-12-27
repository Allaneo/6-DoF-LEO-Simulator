function [phi, lambda, h] = MT_ECEF2GEODETIC(r_ecef,a_world,b_world,e2_world,es2_world)
% MT_ECEF2GEODETIC
% Convierte ECEF -> geodésicas WGS84 usando método de Bowring.
%
% Entrada:
%   r_ecef : vector 3x1 [m]
% Salidas:
%   phi    : latitud geodésica [deg]
%   lambda : longitud [deg]
%   h      : altura sobre el elipsoide [m]

x = r_ecef(1); y = r_ecef(2); z = r_ecef(3);

% Longitud
lambda = atan2(y, x)*180/pi;

% Distancia al eje Z
p = sqrt(x*x + y*y);

% Bowring: ángulo auxiliar
theta = atan2(z*a_world, p*b_world)*180/pi;

st = sind(theta); ct = cosd(theta);

% Latitud geodésica
phi = atan2(z + es2_world*b_world*st^3, p - e2_world*a_world*ct^3)*180/pi;

sphi = sind(phi); cphi = cosd(phi);

% Radio de curvatura N y altura
N = a_world / sqrt(1 - e2_world*sphi*sphi);
h = p/cphi - N;

end
