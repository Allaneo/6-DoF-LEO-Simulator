function r_ecef = MT_GEOCENTRIC2ECEF(r_geo, phi, lambda)
% Entradas:
%   r_geo  : posicion en GEO 3x1 
%   phi    : latitud geocéntrica (rad)
%   lambda : longitud (rad)
% Salida:
%   r_ecef : posicion en ECEF 3x1

cphi = cos(phi); sphi = sin(phi);
cl = cos(lambda); sl = sin(lambda);

% matriz de transformacion
T = [ cphi*cl,  -sphi*cl,  -sl;
      cphi*sl,  -sphi*sl,   cl;
      sphi,      cphi,      0 ];

% multiplicacion matricial
r_ecef = T * r_geo;
end
