%% INIT_parametros.m
% Script de inicialización principal. Contiene los parámetros iniciales
% modificables para iniciar la simulación y llama a los scripts DATA
% necesarios para los distintos modelos MD.

clc; clear;

%% ============================
% Parámetros del Satélite
% ============================
m = 12;                                              % masa del satélite [kg]
L             = [0.3; 0.2; 0.1];                     % longitudes del satélite [m]
A_transversal = [L(2)*L(3); L(1)*L(3); L(1)*L(2)];   % área frontal por eje x,y,z [m^2]

%   Elección de comportamiento del área proyectada. Cannonball para área
%   constante (toma primer valor de A_transversal) ó prisma (calcula área
%   proyectada cuando lo requiera)

SAT_SHAPE = Simulink.Parameter(int32(1));            % 1-Cannonball / 2-Prisma
SAT_SHAPE.CoderInfo.StorageClass = 'SimulinkGlobal';

fprintf('Modelo de Satélite: %s (SAT_SHAPE=%d)\n', ...
    subsref({'Cannonball','Prisma'}, substruct('{}',{double(SAT_SHAPE.Value)})), ...
    double(SAT_SHAPE.Value));

I = [L(2)^2+L(3)^2 0 0
     0 L(1)^2+L(3)^2 0
     0 0 L(2)^2+L(1)^2]*m/12;                         % inercia (prisma homogéneo)

m_res = [1;1;1]*1e-2;           % dipolo residual del satélite
CD = 3;                         % Coeficiente de drag del satélite. Valor estándar 3
Cr = 2;                         % Coefificiente de reflectividad del satélite (entre 1 y 2, valor estándar 1.3)
r_CG_in_G = [0; 0; 0];          % Vector posición del CG respecto centro geométrico
l_minimo_drag = 0.05;           % Brazo de palanca mínimo para momento por drag

%% ============================
% Tiempo inicial (UTC) + Transformación ECI-ECEF
% ============================
utc0 = [2026,2,10,0,0,0];      % [año, mes, día, hora]

% deltaAT = TAI-UTC (leap seconds). Desde 2017 es 37 s.
deltaAT = 37;  % [s]

% Selector de transformación:
% 1: ERA-only
% 2: FIXED_CORR  (ECI->ECEF = R3(-ERA)*C_FIX_val)
XFORM_MODE = Simulink.Parameter(int32(2));
XFORM_MODE.CoderInfo.StorageClass = 'SimulinkGlobal';

xm = double(XFORM_MODE.Value);
if ~(xm==1 || xm==2)
    error("Wrong Input for XFORM_MODE. Use 1 (ERA-only) or 2 (FIXED_CORR).");
end

% --- SIEMPRE calculo EOP y jd_ut_inicial desde la MISMA función ---
% Esto evita que en modo 1 te quede jd_ut_inicial=jd_utc0 (UTC) y ERA mal.
[Cfix_tmp, infoCfix] = MD_build_Cfix(utc0, deltaAT);

% Exporto SIEMPRE lo coherente para el resto del proyecto
deltaUT1  = infoCfix.deltaUT1;   % [s]
xp_rad        = infoCfix.xp_rad;     % [rad]
yp_rad        = infoCfix.yp_rad;     % [rad]
jd_utc0       = infoCfix.jd_utc0;    % JD(UTC) t0
jd_ut_inicial = infoCfix.jd_ut10;    % JD(UT1) t0  <<< CLAVE para ERA

% Defino C_FIX_val según modo (pero SIEMPRE como Simulink.Parameter global)
if xm == 2
    C_FIX_num = Cfix_tmp;        % FIXED_CORR
else
    C_FIX_num = eye(3);          % ERA-only
end

C_FIX_val = Simulink.Parameter(C_FIX_num);
C_FIX_val.CoderInfo.StorageClass = 'SimulinkGlobal';

% Print transformación (arreglado)
labels = {'ERA-only','FIXED_CORR'};
fprintf('Tipo de transformación ECI-ECEF: %s (XFORM_MODE=%d)\n', labels{xm}, xm);

% Print EOP / JD
fprintf("EOP t0: deltaUT1=%.6f s | xp=%.3e rad | yp=%.3e rad\n", deltaUT1, xp_rad, yp_rad);
fprintf("JD: utc0=%.8f | ut10=%.8f | diff=%.3e days\n", jd_utc0, jd_ut_inicial, jd_ut_inicial-jd_utc0);

% --- CHECK DURO en t0 (solo informativo; en modo 2 debería ser ~0 numérico) ---
if xm == 2
    % ERA0 usando EXACTAMENTE el mismo pipeline que usa la simulación:
    [ERA0_sim, ~] = MD_era(0, jd_ut_inicial);

    c = cos(ERA0_sim); s = sin(ERA0_sim);
    R0 = [ c  s  0;
          -s  c  0;
           0  0  1];  % R3(-ERA0)

    C_full0 = dcmeci2ecef('IAU-2000/2006', utc0, deltaAT, deltaUT1, [xp_rad yp_rad]);

    Cfix_num = double(C_FIX_val.Value);
    fprintf("CHECK t0: ||R0*Cfix - Cfull0||_F = %.3e\n", norm(R0*Cfix_num - C_full0,'fro'));
end

%% ============================
% Condiciones Iniciales de órbita / actitud
% ============================
% Selector de caso orbital:
% 0 = CUSTOM
% 1 = 700 km circular ecuatorial
% 2 = 700 km circular polar
% 3 = 700 km circular i=98 deg, RAAN=40 deg, nu=90 deg
ORBIT_CASE = int32(3);

% --- CUSTOM CASE (SE USA SOLAMENTE PARA  ORBIT_CASE == 0) ---
posicion_inicial  = [633;-755;7010]*1e3;   % [m] ECI
velocidad_inicial = [-5749;-4824;0];       % [m/s] ECI

% Actitud inicial: [q0;q1;q2;q3; wx;wy;wz]
CI_angular = [0;0;0;0; 0;0;0];

% --- Presets ---
switch ORBIT_CASE
    case 0
        case_name = 'CUSTOM';

    case 1
        case_name = 'ECI_700km_CIRC_EQUATORIAL';
        posicion_inicial  = [7078.137; 0; 0]*1e3;
        velocidad_inicial = [0; 7504.28649; 0];

    case 2
        case_name = 'ECI_700km_CIRC_POLAR';
        posicion_inicial  = [5004.99867; 0; 5004.99867]*1e3;
        velocidad_inicial = [-5306.33187; 0; 5306.33187];

    case 3
        case_name = 'ECI_700km_CIRC_i98_Om40_nu90';
        posicion_inicial  = [633.20125; -754.61987; 7009.25306]*1e3;
        velocidad_inicial = [-5748.61697; -4823.66238; 0];

    otherwise
        error('ORBIT_CASE inválido. Use 0..3.');
end

CI_lineal = [posicion_inicial; velocidad_inicial];
fprintf('Caso orbital: %s (ORBIT_CASE=%d)\n', case_name, ORBIT_CASE);

%% ============================
% Parámetros de la Tierra / modelos
% ============================

% --- Selección de Modelo Gravitatorio ---
GRAV_MODEL = 3;                                  % 1: N4 / 2: J2 / 3: SIMPLE

% --- Selección (y switch) de Modelo Atmosférico para Drag ---
ATM_SEL = Simulink.Parameter(int32(3));          % 1: US / 2: EXP / 3: NONE
ATM_SEL.CoderInfo.StorageClass = 'SimulinkGlobal';

DATA_earth;
fprintf('Modelo atmosférico: %s (ATM_SEL=%d)\n', ...
    subsref({'US','EXP','NONE'}, substruct('{}',{double(ATM_SEL.Value)})), ...
    double(ATM_SEL.Value));

% --- Switch ON-OFF de modelo geomagnético (torque geomagnético) ---
GEOMAG_SEL = Simulink.Parameter(int32(1));       % 1: ON / 2: OFF
GEOMAG_SEL.CoderInfo.StorageClass = 'SimulinkGlobal';

fprintf('Modelo GEOMAGNETICO: %s (GEOMAG_SEL=%d)\n', ...
    subsref({'ON','OFF'}, substruct('{}',{double(GEOMAG_SEL.Value)})), ...
    double(GEOMAG_SEL.Value));

%% ============================
% Parámetros para posición solar / SRP
% ============================
DATA_sun;

% --- Switch ON-OFF para SRP ---
SRP_SEL = Simulink.Parameter(int32(2));          % 1: ON / 2: OFF
SRP_SEL.CoderInfo.StorageClass = 'SimulinkGlobal';

fprintf('Modelo SRP: %s (SRP_SEL=%d)\n', ...
    subsref({'ON','OFF'}, substruct('{}',{double(SRP_SEL.Value)})), ...
    double(SRP_SEL.Value));
