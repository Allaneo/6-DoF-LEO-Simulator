% DATA_earth.m
% Script de data sobre el modelo terrestre para el cálculo gravitacional.
% Permite elegir (en INIT_parametros) entre:
%   - "N4"     : campo hasta grado/orden 4 (Cbar,Sbar provistos)
%   - "J2"     : solo término C20 (equivalente a J2)
%   - "SIMPLE" : punto-masa (solo GM)

K = [1.000000000000000	0	0	0
1.732050807568877	1.732050807568877	0	0
2.236067977499790	1.290994448735806	0.645497224367903	0
2.645751311064591	1.080123449734643	0.341565025531987	0.139443337755679];

a_world   = 6378137.2;                                   % semieje mayor [m]
f_world   = 1/298.257223563;                             % aplanamiento
b_world   = a_world*(1-f_world);                         % semieje menor [m]
e2_world  = f_world*(2 - f_world);                       % excentricidad^2
es2_world = (a_world^2 - b_world^2)/(b_world^2);         % segunda excentricidad^2
GM = 3.986004418e14;

Cbar = zeros(4);
Sbar = zeros(4);

switch GRAV_MODEL
    case 1
        Cbar = [0	0	0	0
        0	0	0	0
        -4.841651437908150e-04	-2.066155090741760e-10	2.439383573283130e-06	0
        9.571612070934730e-07	2.030462010478640e-06	9.047878948095281e-07	7.213217571215680e-07];
        
        Sbar = [0	0	0	0
        0	0	0	0
        0	1.384413891379790e-09	-1.400273703859340e-06	0
        0	2.482004158568720e-07	-6.190054751776180e-07	1.414349261929410e-06];
        disp("Modelo gravitacional: N4")

    case 2
        % C20 normalizado (relación con J2: C20 = -J2/sqrt(5))
        Cbar(3,1) = -4.841651437908150e-04;
        disp("Modelo gravitacional: J2")
    case 3
        disp("Modelo gravitacional: simple")
end