% parámetros iniciales para iniciar la simulación

clc; clear;
%% Parámetros del Satélite
m = 1000;                       % masa del satélite
I = [100 0 0
    0 200 0
    0 0 150];                   % matriz de inercia del satélite
A_transversal = [1;0.1;0.1];    % área transversal
m_res = [1;1;1]*1e-2;
CD = 2.2;
Cr = 2;
CP = [0,0,0.01];                % posición del centro de presiones respecto del centro de gravedad

%% Condiciones Iniciales
posicion_inicial = [6988;0;0]*1e3; % metros (en total, km adentro del corchete)
velocidad_inicial = [0;7600;0];    % m/s

CI_lineal = [posicion_inicial;velocidad_inicial];
CI_angular = [1;0;0;0;0;0;0];

%% Parámetros Gravitacionales (no modificar)
Cbar = [0	0	0	0
0	0	0	0
-4.841651437908150e-04	-2.066155090741760e-10	2.439383573283130e-06	0
9.571612070934730e-07	2.030462010478640e-06	9.047878948095281e-07	7.213217571215680e-07];

Sbar = [0	0	0	0
0	0	0	0
0	1.384413891379790e-09	-1.400273703859340e-06	0
0	2.482004158568720e-07	-6.190054751776180e-07	1.414349261929410e-06];

K = [1.000000000000000	0	0	0
1.732050807568877	1.732050807568877	0	0
2.236067977499790	1.290994448735806	0.645497224367903	0
2.645751311064591	1.080123449734643	0.341565025531987	0.139443337755679];

a_world   = 6378137.2;                 % semieje mayor [m]
f_world   = 1/298.257223563;           % aplanamiento
b_world   = a_world*(1-f_world);       % semieje menor [m]
e2_world  = f_world*(2 - f_world);     % excentricidad^2
es2_world = (a_world^2 - b_world^2)/(b_world^2);         % segunda excentricidad^2

% Cbar = zeros(4);
% Sbar = zeros(4);
% K = zeros(4);
% K(1,1) = 1;

GM = 3.986004418e14;

% Parametros gravitacionales "modificables"
delta =0.1; % distancia para calcular la derivada espacial para el momento gravitacional

%% Parámetros para presión solar (no modificar)
solar_omegas = 282.94; %Omega+omega
solar_epsilon = 23.43929111;
solar_cosep = cosd(solar_epsilon);
solar_sinep = sind(solar_epsilon);