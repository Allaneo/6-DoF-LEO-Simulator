function r_ecef = MT_ECI2ECEF(r_eci, era)
%   Convierte posición de ECI a ECEF
%   r_ecef: vector de posicion en ECEF (3x1)
%   r_eci: vector de posición en ECI (3x1)
%   era: earth rotation angle [0,2pi)

% Matriz de rotación ECEF->ECI
c = cos(era);
s = sin(era);
R = [c, s, 0;
     -s,  c, 0;
     0,  0, 1];

% Para ECI->ECEF usamos la transpuesta
r_ecef = R * r_eci;
end