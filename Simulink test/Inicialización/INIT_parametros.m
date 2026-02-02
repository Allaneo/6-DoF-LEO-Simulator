%% INIT_parametros.m
% Script de inicialización principal. Contiene los parámetros iniciales
% modificables para iniciar la simulación y llama a los scripts DATA
% necesarios para los distintos modelos MD. Variables modificables:
% - Parámetros del Satélite: m, L, A_transversal, I, m_res, CD, Cr, r_CG_in_G, l_minimo_drag
% - Condiciones iniciales: posicion_inicial, velocidad_inicial, jd_ut_inicial
% - Parámetros de la Tierra: GRAV_MODEL, ATM_SEL


clc; clear;
%% Parámetros del Satélite
m = 12;                                              % masa del satélite
L             = [0.3; 0.2; 0.1];                     % longitudes del satélite
A_transversal = [L(2)*L(3); L(1)*L(3); L(1)*L(2)];   % área frontal por eje x,y,z

I = [L(2)^2+L(3)^2 0 0
    0 L(1)^2+L(3)^2 0
    0 0 L(2)^2+L(1)^2]*m/12;                         % matriz de inercia del satélite asumiendo prisma rectangular homogéneo

m_res = [1;1;1]*1e-2;           % nombre?      
CD = 3;                         % Coeficiente de resistencia (drag). 3 para placa plana según DSMC.
Cr = 2;                         % Coeficiente de reflectividad (SRP). 1.3 como valor estándar; 2 como valor máximo.
r_CG_in_G = [0; 0; 0];          % posición del centro de gravedad RESPECTO del centro geométrico
l_minimo_drag = 0.05;           % brazo de palanca MÍNIMO para el cálculo de momento por drag.

%% Condiciones Iniciales
utc0 = [2026,1,20,0,0,0];
jd_utc0 = juliandate(utc0(1),utc0(2),utc0(3),utc0(4),utc0(5),utc0(6));

%% === EOP + C_FIX (para FIXED_CORR) ===
XFORM_MODE = Simulink.Parameter(int32(2));              % 1: ERA-only, 2: FIXED_CORR 
XFORM_MODE.CoderInfo.StorageClass = 'SimulinkGlobal';

if XFORM_MODE.Value == 2
% deltaAT = TAI-UTC (leap seconds). Desde 2017 viene siendo 37 s.
deltaAT = 37;  % [s]  (TAI-UTC)

% --- leer EOP (UT1-UTC, xp, yp) desde finals2000A.data ---
[deltaUT1, xp_rad, yp_rad] = MD_getEOP_finals2000A(utc0);

% JD(UT1) en t0 para tu ERA
jd_ut_inicial = jd_utc0 + deltaUT1/86400.0;

% DCM "full" solo en t0 (Aerospace Toolbox)
polarmotion = [xp_rad yp_rad];  % [rad rad]
C_full0 = dcmeci2ecef('IAU-2000/2006', utc0, deltaAT, deltaUT1, polarmotion);

% ERA0 con tu función
[ERA0,~] = MD_era(0, jd_ut_inicial);
c0 = cos(ERA0); s0 = sin(ERA0);
R3_0 = [ c0  s0  0;
        -s0  c0  0;
         0   0   1];

% C_FIX tal que: C_full0 ≈ C_FIX * R3_0
C_FIX_val = C_full0 * (R3_0.');

elseif XFORM_MODE.Value == 1
C_FIX_val = eye(3);

else
    error("Wrong Input for XFORM_MODE. Read the script.")
end
fprintf('Tipo de transformación ECI-ECEF: %s (XFORM_MODE=%d)\n', XFORM_MODE.Value);


%% Condiciones Iniciales
% Selector de caso orbital:
% 0 = CUSTOM (usa los vectores que definas abajo)
% 1 = 700 km circular ecuatorial (control)
% 2 = 700 km circular polar (max latitud)
% 3 = 700 km circular i=98 deg (tipo SSO), con RAAN=40 deg y nu=90 deg
ORBIT_CASE = int32(2);

% --- CUSTOM (se usa solo si ORBIT_CASE==0) ---
posicion_inicial  = [633;-755;7010]*1e3;   % [m] en terna ECI (metros en total, km adentro del corchete). semi eje mayor de referencia a = 6378 km.
velocidad_inicial = [-5749;-4824;0];       % [m/s] ECI

% Actitud inicial: [q0;q1;q2;q3; wx;wy;wz]
% Recomendación para "encender" el gravity-gradient: offset de 45 deg sobre eje Y (BODY)
q_offset_45deg_y = [0.9238795325; 0; 0.3826834324; 0]; % cos(22.5), 0, sin(22.5), 0
CI_angular = [q_offset_45deg_y; 0;0;0];                % Actitud y velocidad angular inicial del satélite.

% --- Presets (sobrescriben posicion/velocidad si ORBIT_CASE ~= 0) ---
switch ORBIT_CASE
    case 0
        case_name = 'CUSTOM';

    case 1  % 700 km circular ecuatorial (r=[R,0,0], v=[0,V,0])
        case_name = 'ECI_700km_CIRC_EQUATORIAL';
        posicion_inicial  = [7078.137; 0; 0]*1e3;        % [m]
        velocidad_inicial = [0; 7504.28649; 0];          % [m/s]

    case 2  % 700 km circular polar (nu=45 deg, Omega=0, i=90)
        case_name = 'ECI_700km_CIRC_POLAR';
        posicion_inicial  = [5004.99867; 0; 5004.99867]*1e3;
        velocidad_inicial = [-5306.33187; 0; 5306.33187];

    case 3  % 700 km circular i=98 deg, Omega=40 deg, nu=90 deg
        case_name = 'ECI_700km_CIRC_i98_Om40_nu90';
        posicion_inicial  = [633.20125; -754.61987; 7009.25306]*1e3;
        velocidad_inicial = [-5748.61697; -4823.66238; 0];

    otherwise
        error('ORBIT_CASE inválido. Use 0..3.');
end

CI_lineal = [posicion_inicial; velocidad_inicial];      % Posición inicial y velocidad inicial del satélite.

fprintf('Caso orbital: %s (ORBIT_CASE=%d)\n', case_name, ORBIT_CASE);
%% Parámetros de la Tierra 
GRAV_MODEL = 1;                                  % Modelos gravitatorios disponibles: "1: N4" / "2: J2" / "3: SIMPLE"
ATM_SEL = Simulink.Parameter(int32(1));          % Modelo atmosférico para cálculo de densidad. 1: US1962/1976 / 2: modelo exponencial
ATM_SEL.CoderInfo.StorageClass = 'SimulinkGlobal';
DATA_earth;
fprintf('Modelo atmosférico: %s (ATM_SEL=%d)\n', subsref({'US','EXP'}, substruct('{}',{double(ATM_SEL.Value)})), ATM_SEL.Value);

%% Parámetros para posición solar
DATA_sun;



%% === Local function: parse finals2000A.data ===
function [dut1_sec, xp_rad, yp_rad] = MD_getEOP_finals2000A(utc_vec)
% Lee UT1-UTC y polar motion (xp, yp) para la fecha utc_vec desde:
% https://maia.usno.navy.mil/ser7/finals2000A.data
% Formato de columnas según MathWorks (aeroReadIERSData). :contentReference[oaicite:3]{index=3}

    url = 'https://maia.usno.navy.mil/ser7/finals2000A.data';
    local = fullfile(tempdir, 'finals2000A.data');

    if ~isfile(local)
        try
            websave(local, url);
        catch
            error('No pude descargar finals2000A.data. Bajalo manualmente de USNO y ponelo en: %s', local);
        end
    end

    y = utc_vec(1); m = utc_vec(2); d = utc_vec(3);

    fid = fopen(local,'r');
    if fid < 0, error('No pude abrir %s', local); end

    dut1_sec = 0; xp_rad = 0; yp_rad = 0;
    found = false;

    % finals2000A.data es fixed-width. Usamos las columnas oficiales:
    % Year(1-2), Month(3-4), Day(5-6)
    % PM-x (19-27) arcsec, PM-y (38-46) arcsec, UT1-UTC (59-68) sec :contentReference[oaicite:4]{index=4}
    while ~feof(fid)
        line = fgetl(fid);
        if ~ischar(line) || numel(line) < 68, continue; end

        yy = str2double(line(1:2));
        mm = str2double(line(3:4));
        dd = str2double(line(5:6));

        if isnan(yy) || isnan(mm) || isnan(dd), continue; end

        % Regla del archivo: año real depende de MJD; para fechas modernas es 2000+yy.
        % Para 2026 esto es correcto.
        yyyy = 2000 + yy;

        if yyyy==y && mm==m && dd==d
            pmx_asec = str2double(line(19:27));
            pmy_asec = str2double(line(38:46));
            dut1_sec = str2double(line(59:68)); % UT1-UTC [s]

            asec2rad = (pi/180.0)/3600.0;
            xp_rad = pmx_asec * asec2rad;
            yp_rad = pmy_asec * asec2rad;

            found = true;
            break;
        end
    end
    fclose(fid);

    if ~found
        error('No encontré EOP para %04d-%02d-%02d en finals2000A.data', y,m,d);
    end
end