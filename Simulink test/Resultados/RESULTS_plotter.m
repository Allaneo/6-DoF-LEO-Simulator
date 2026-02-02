%% RESULTS_plotter_GMAT.m  (SCRIPT)
% Plotter + validación contra GMAT (configurable, sin inferir formatos).
%
% Fuente de datos:
%  - logsout (Dataset) si existe, si no -> SDI (último run)
%
% Señales esperadas (según tu captura / naming actual):
%  x_lineal (6x1), x_angular (7x1), v_rel (3x1), t, shadow, rho,
%  r_solar (3x1), r_cp (3x1), phi, lambda,
%  Fuerza Drag/Solar/Gravedad/Total ECI (3x1),
%  Momento Drag/Geomagnético TOTAL/Gravedad/Solar/Total BODY (3x1),
%  era, Aproj_drag
%
% NOTA: La comparación GMAT NO se ejecuta si no configurás GMAT_FILE + GMAT_COLMAP.

clc;
nl = sprintf('\n');

%% =========================
% 0) Parámetros / convenciones
% =========================
omegaE = 7.2921150e-5; % rad/s

% Quaternion:
%   quat_mode = 'CIB' -> q define C_ib (body->inertial): v_I = C_ib * v_B
%   quat_mode = 'CBI' -> q define C_bi (inertial->body): v_B = C_bi * v_I
quat_mode = 'CIB';

% Eje del cuerpo a comparar contra nadir (en body frame)
body_axis_for_nadir = [0 0 1];

% Diagnósticos orbitales
PLOT_ORBIT_DIAGNOSTICS = true;
muEarth = 3.986004418e14; % [m^3/s^2]

% Ground track: si usás phi/lambda logueados, definí unidades explícitas (NO infiero)
% Opciones: 'rad' o 'deg'
PHI_LAMBDA_UNITS = 'rad';

%% =========================
% GMAT: CONFIGURACIÓN OBLIGATORIA (si querés comparación)
% =========================
% 1) Path al ReportFile/CSV/TXT exportado de GMAT:
GMAT_FILE = '';  % <-- EJ: 'C:\...\gmat_report.csv'

% 2) Mapeo de columnas (NOMBRES EXACTOS tal como aparecen en el header de tu archivo GMAT)
%    Si alguno queda vacío, NO se compara (y NO se infiere).
GMAT_COLMAP = struct();
GMAT_COLMAP.t  = '';   % tiempo en segundos desde t0 (recomendado) o algo que vos conviertas a segundos
GMAT_COLMAP.rx = '';   % posición ECI x [m]
GMAT_COLMAP.ry = '';   % posición ECI y [m]
GMAT_COLMAP.rz = '';   % posición ECI z [m]
GMAT_COLMAP.vx = '';   % velocidad ECI x [m/s]
GMAT_COLMAP.vy = '';   % velocidad ECI y [m/s]
GMAT_COLMAP.vz = '';   % velocidad ECI z [m/s]

% (Opcional) densidad si la exportaste desde GMAT:
GMAT_COLMAP.rho = '';  % [kg/m^3]

%% =========================
% 1) Señales (NOMBRES EXACTOS según tu captura)
% =========================
SIG.x_lineal  = 'x_lineal';
SIG.x_angular = 'x_angular';

SIG.v_rel   = 'v_rel';
SIG.t       = 't';
SIG.shadow  = 'shadow';
SIG.rho     = 'rho';

SIG.r_solar = 'r_solar';
SIG.r_cp    = 'r_cp';

SIG.phi     = 'phi';
SIG.lambda  = 'lambda';

SIG.M_total = 'Momento Total BODY';
SIG.M_grav  = 'Momento Gravedad';
SIG.M_geo   = 'Momento Geomagnético TOTAL';
SIG.M_drag  = 'Momento Drag';
SIG.M_solar = 'Momento Solar';

SIG.F_total = 'Fuerza Total ECI';
SIG.F_solar = 'Fuerza Solar';
SIG.F_grav  = 'Fuerza Gravedad';
SIG.F_drag  = 'Fuerza Drag';

SIG.era     = 'era';
SIG.Aproj   = 'Aproj_drag';

sigKeys = { ...
    'v_rel','t','shadow','rho','r_solar','r_cp','phi','lambda', ...
    'M_total','M_grav','M_geo','M_drag','M_solar', ...
    'F_total','F_solar','F_grav','F_drag', ...
    'era','Aproj' ...
};

sigNames = { ...
    SIG.v_rel, SIG.t, SIG.shadow, SIG.rho, SIG.r_solar, SIG.r_cp, SIG.phi, SIG.lambda, ...
    SIG.M_total, SIG.M_grav, SIG.M_geo, SIG.M_drag, SIG.M_solar, ...
    SIG.F_total, SIG.F_solar, SIG.F_grav, SIG.F_drag, ...
    SIG.era, SIG.Aproj ...
};

%% =========================
% 2) Fuente: logsout o SDI
% =========================
SRC_MODE = '';
logs = [];

sdiNames = {};
sdiTS    = {};
sdiRunID = [];

if exist('logsout','var') == 1
    try
        if isa(logsout,'Simulink.SimulationData.Dataset') && logsout.numElements > 0
            logs = logsout;
            SRC_MODE = 'workspace.logsout (Dataset)';
        end
    catch
    end
end

if isempty(SRC_MODE)
    runIDs = Simulink.sdi.getAllRunIDs;
    if isempty(runIDs)
        error(['No encontré logsout en workspace y tampoco hay runs en SDI.' nl ...
               'No hay datos registrados para graficar.']);
    end

    sdiRunID = runIDs(end);
    r = Simulink.sdi.getRun(sdiRunID);

    try
        ids = r.getSignalIDs;
    catch
        error(['Hay runs en SDI pero no puedo listar señales (getSignalIDs no disponible).' nl ...
               'Abrí el SDI y verificá que las señales estén registradas.']);
    end

    nS = numel(ids);
    sdiNames = cell(nS,1);
    sdiTS    = cell(nS,1);

    for i = 1:nS
        sigObj = [];
        try
            sigObj = Simulink.sdi.getSignal(ids(i));
        catch
            sigObj = [];
        end

        if isempty(sigObj)
            sdiNames{i} = '';
            sdiTS{i} = [];
            continue;
        end

        try
            sdiNames{i} = sigObj.Name;
        catch
            sdiNames{i} = '';
        end

        sdiTS{i} = local_get_sdi_timeseries(sigObj);
    end

    SRC_MODE = sprintf('SDI (último run ID %d)', sdiRunID);
end

fprintf('Fuente de datos: %s\n', SRC_MODE);

%% =========================
% 3) Cargar señales en struct data (strcmp/strcmpi)
% =========================
data = struct();

for k = 1:numel(sigKeys)
    key  = sigKeys{k};
    name = sigNames{k};

    data.(key).ok   = false;
    data.(key).name = name;
    data.(key).t    = [];
    data.(key).X    = [];

    ts = [];

    if ~isempty(logs)
        ts = local_get_dataset_timeseries_exact(logs, name);
    else
        ts = local_get_sdi_timeseries_by_name(sdiNames, sdiTS, name);
    end

    if ~isempty(ts)
        try
            [t, X] = local_read_ts(ts);
            data.(key).ok = true;
            data.(key).t  = t;
            data.(key).X  = X;
        catch
            data.(key).ok = false;
        end
    end
end

%% =========================
% 4) Leer x_lineal (OBLIGATORIO): r,v en ECI
% =========================
xTS = [];
if ~isempty(logs)
    xTS = local_get_dataset_timeseries_exact(logs, SIG.x_lineal);
else
    xTS = local_get_sdi_timeseries_by_name(sdiNames, sdiTS, SIG.x_lineal);
end

if isempty(xTS)
    error(['No encontré la señal "' SIG.x_lineal '" en logsout/SDI.' nl ...
           'Asegurate de que esté registrada exactamente con ese nombre.']);
end

[tx, Xstate] = local_read_ts(xTS);
if size(Xstate,2) < 6
    error([SIG.x_lineal ' debe contener al menos 6 estados: r(1:3) y v(4:6).']);
end

rECI = Xstate(:,1:3);
vECI = Xstate(:,4:6);

xNameUsed = SIG.x_lineal;

%% =========================
% 5) Leer x_angular (OPCIONAL): [q0 q1 q2 q3 wx wy wz]
% =========================
xAngTS = [];
if ~isempty(logs)
    xAngTS = local_get_dataset_timeseries_exact(logs, SIG.x_angular);
else
    xAngTS = local_get_sdi_timeseries_by_name(sdiNames, sdiTS, SIG.x_angular);
end

if isempty(xAngTS)
    fprintf(['Aviso: no encontré "' SIG.x_angular '". Se saltean figuras de actitud.' nl]);
    haveAtt = false;
else
    [tA, Xang] = local_read_ts(xAngTS);
    haveAtt = true;

    if size(Xang,2) < 7
        fprintf(['Aviso: "' SIG.x_angular '" tiene menos de 7 componentes. Se saltean figuras de actitud.' nl]);
        haveAtt = false;
    end
end

if haveAtt
    q_norm = Xang(:,1:4);
    wB    = Xang(:,5:7);

    qn_raw = sqrt(sum(q_norm.^2,2));
    q_norm_err = qn_raw - 1;

    omega_mag = sqrt(sum(wB.^2,2));

    q_plot = local_quat_continuity(q_norm);

    n_safe = qn_raw; n_safe(n_safe < 1e-12) = 1;
    q_used = q_plot ./ [n_safe n_safe n_safe n_safe];

    [c11,c12,c13,c21,c22,c23,c31,c32,c33] = local_quat_to_dcm_cib(q_used, quat_mode);

    [yaw, pitch, roll] = local_dcm_to_euler_zyx(c11,c21,c31,c32,c33);
    yaw_deg   = local_unwrap_deg(yaw   * 180/pi);
    pitch_deg = local_unwrap_deg(pitch * 180/pi);
    roll_deg  = local_unwrap_deg(roll  * 180/pi);

    axb = body_axis_for_nadir(:);
    axb = axb / max(norm(axb), eps);

    axI_x = c11*axb(1) + c12*axb(2) + c13*axb(3);
    axI_y = c21*axb(1) + c22*axb(2) + c23*axb(3);
    axI_z = c31*axb(1) + c32*axb(2) + c33*axb(3);
    axI = [axI_x axI_y axI_z];

    rx = interp1(tx, rECI(:,1), tA, 'linear', 'extrap');
    ry = interp1(tx, rECI(:,2), tA, 'linear', 'extrap');
    rz = interp1(tx, rECI(:,3), tA, 'linear', 'extrap');
    rI = [rx ry rz];
    rnorm = sqrt(sum(rI.^2,2));
    rhat = rI ./ [rnorm rnorm rnorm];

    nadir = -rhat;

    cosang = sum(axI .* nadir, 2);
    cosang = max(min(cosang, 1), -1);
    nadir_err = acos(cosang) * 180/pi;
end

%% =========================
% 6) Derivada: Gravedad NO radial total (perp a r)
% =========================
data.F_grav_nonrad.ok   = false;
data.F_grav_nonrad.name = 'Fuerza Gravedad NO radial (ECI)';
data.F_grav_nonrad.t    = [];
data.F_grav_nonrad.X    = [];

if data.F_grav.ok
    tFg = data.F_grav.t(:);
    Fg  = data.F_grav.X;

    % Fg podría venir con dims extra -> asegurar NxM y luego quedarnos con 3 primeras si existen
    Fg2 = squeeze(Fg);
    if size(Fg2,1) ~= numel(tFg) && size(Fg2,2) == numel(tFg)
        Fg2 = permute(Fg2, [2 1 3:ndims(Fg2)]);
    end
    if size(Fg2,1) ~= numel(tFg)
        try
            Fg2 = reshape(Fg2, numel(tFg), []);
        catch
            Fg2 = Fg2;
        end
    else
        Fg2 = reshape(Fg2, numel(tFg), []);
    end

    if size(Fg2,2) >= 3
        Fg3 = Fg2(:,1:3);
    else
        % si no hay 3 componentes no puedo derivar no radial
        Fg3 = [];
    end

    if ~isempty(Fg3)
        rx = interp1(tx, rECI(:,1), tFg, 'linear', 'extrap');
        ry = interp1(tx, rECI(:,2), tFg, 'linear', 'extrap');
        rz = interp1(tx, rECI(:,3), tFg, 'linear', 'extrap');
        rI = [rx ry rz];

        rnorm = sqrt(sum(rI.^2, 2));
        bad = (rnorm <= 0);
        rnorm(bad) = NaN;

        rhat = rI ./ [rnorm rnorm rnorm];

        Frad_mag = sum(Fg3 .* rhat, 2);
        F_rad    = rhat .* [Frad_mag Frad_mag Frad_mag];

        F_perp = Fg3 - F_rad;

        data.F_grav_nonrad.ok = true;
        data.F_grav_nonrad.t  = tFg;
        data.F_grav_nonrad.X  = F_perp;
    end
end

%% =========================
% 7) Sombra: segmentos [tStart tEnd]
% =========================
shadowSegments = [];
if data.shadow.ok
    tt = data.shadow.t(:);
    s  = data.shadow.X(:) > 0.5;

    if numel(tt) == numel(s) && numel(tt) > 1
        ds = diff([false; s; false]);
        i1 = find(ds == 1);
        i2 = find(ds == -1) - 1;
        shadowSegments = [tt(i1) tt(i2)];
    end
end

%% =========================
% 8) FIGURA 1: Fuerzas perturbadoras (magnitud)
% =========================
forceKeys   = {'F_drag','F_solar','F_grav_nonrad'};
forceLabels = {'Drag','Solar','Gravedad NO radial'};

figure('Name','Perturbaciones - Fuerzas (Magnitud)','Color','w');
ax1 = axes; hold(ax1,'on'); grid(ax1,'on');

hF = []; lblF = {};

for i = 1:numel(forceKeys)
    k = forceKeys{i};

    if isfield(data,k) && data.(k).ok
        t = data.(k).t(:);
        X = data.(k).X;

        mag = local_mag_timeseries_samplewise(X, t); % Nx1
        mag = max(mag, realmin);                     % evita log(0)

        mag = local_mag_timeseries_samplewise(X, t); % Nx1

if strcmp(k,'F_solar') && data.shadow.ok
    sh = interp1(data.shadow.t(:), double(data.shadow.X(:)>0.5), t, 'nearest', 'extrap') > 0.5;
    mag(sh) = NaN;              % cortar en sombra
end

mag(mag<=0) = NaN;              % log-safe general
h = semilogy(ax1, t, mag, 'LineWidth', 1.4);

        h = semilogy(ax1, t, mag, 'LineWidth', 1.4);
        hF(end+1) = h; %#ok<SAGROW>
        lblF{end+1} = forceLabels{i}; %#ok<SAGROW>
    end
end

set(ax1,'YScale','log');   % fuerza log sí o sí
grid(ax1,'on');            % opcional: minor grid ayuda a "ver" el log
set(ax1,'YMinorGrid','on');

if isempty(hF)
    text(0.02,0.5,'No hay señales suficientes para fuerzas perturbadoras.','Units','normalized');
else
    xlabel(ax1,'Tiempo [s]');
    ylabel(ax1,'|F| (norma)');
    title(ax1,'Fuerzas perturbadoras (sin componente radial de gravedad)');
    local_apply_shadow(ax1, shadowSegments);
    legend(ax1, hF, lblF, 'Location','best');
end

%% =========================
% 8b) FIGURA 2: IDEM sin fuerza gravitacional
% =========================
figure('Name','Perturbaciones - Fuerzas (Magnitud)','Color','w');
ax1 = axes; hold(ax1,'on'); grid(ax1,'on');

hF = []; lblF = {};

for i = 1:numel(forceKeys)-1
    k = forceKeys{i};

    if isfield(data,k) && data.(k).ok
        t = data.(k).t;
        X = data.(k).X;

        mag = local_mag_timeseries_samplewise(X, t); % Nx1 garantizado

        h = semilogy(ax1, t, mag, 'LineWidth', 1.4);
        hF(end+1) = h; %#ok<SAGROW>
        lblF{end+1} = forceLabels{i}; %#ok<SAGROW>
    end
end

if isempty(hF)
    text(0.02,0.5,'No hay señales suficientes para fuerzas perturbadoras.','Units','normalized');
else
    xlabel(ax1,'Tiempo [s]');
    ylabel(ax1,'|F| (norma)');
    title(ax1,'Fuerzas perturbadoras (sin componente radial de gravedad)');
    local_apply_shadow(ax1, shadowSegments);
    legend(ax1, hF, lblF, 'Location','best');
end


%% =========================
% 10) Ground track
%  Preferencia:
%   A) Si existen phi/lambda logueados -> los uso
%   B) Si no existen -> calculo desde rECI usando ERA (+ C_FIX/C_FIX_val si corresponde)
% =========================
useLoggedPhiLambda = data.phi.ok && data.lambda.ok;

if useLoggedPhiLambda
    tpos = data.phi.t(:);
    phi  = data.phi.X(:);

    % lambda puede tener su propio timebase -> interpolar a tpos (determinista)
    lam = interp1(data.lambda.t(:), data.lambda.X(:), tpos, 'linear', 'extrap');

    if strcmpi(PHI_LAMBDA_UNITS,'rad')
        lat_deg = phi * 180/pi;
        lon_deg = lam * 180/pi;
    elseif strcmpi(PHI_LAMBDA_UNITS,'deg')
        lat_deg = phi;
        lon_deg = lam;
    else
        error('PHI_LAMBDA_UNITS debe ser ''rad'' o ''deg''.');
    end

else
    tpos = tx(:);

    % --- ERA(t) ---
    haveEra = isfield(data,'era') && data.era.ok;
    if haveEra
        tEra = data.era.t(:);
        era  = data.era.X(:);
        theta = interp1(tEra, era, tpos, 'linear', 'extrap'); % [rad]
    else
        if exist('jd_ut_inicial','var') == 1
            if exist('MD_era','file') == 2 || exist('MD_era','file') == 6
                theta = zeros(size(tpos));
                for i = 1:numel(tpos)
                    theta(i) = MD_era(tpos(i), jd_ut_inicial);
                end
            else
                theta = local_era_from_jd_ut1(jd_ut_inicial, tpos);
                fprintf('Aviso: no existe MD_era en path. Uso fórmula IERS/IAU para ERA.\n');
            end
        else
            theta = omegaE * tpos;
            fprintf('Aviso: no hay señal era ni jd_ut_inicial. Uso theta=omegaE*t (aprox).\n');
        end
    end

    % --- R3(theta): ECEF_tmp = R3(theta)*ECI ---
    ct = cos(theta); st = sin(theta);

    xeci = rECI(:,1);
    yeci = rECI(:,2);
    zeci = rECI(:,3);

    xecef = ct.*xeci + st.*yeci;
    yecef = -st.*xeci + ct.*yeci;
    zecef = zeci;

    % --- FIXED_CORR: aplicar C_FIX o C_FIX_val si existen ---
    useCFIX = false;
    CF = [];

    if exist('XFORM_MODE','var') == 1
        try
            xm = XFORM_MODE;
            if isa(xm,'Simulink.Parameter')
                xm_val = double(xm.Value);
            else
                xm_val = double(xm);
            end

            if xm_val == 2
                if exist('C_FIX','var') == 1
                    if isa(C_FIX,'Simulink.Parameter'), CF = C_FIX.Value; else, CF = C_FIX; end
                    useCFIX = true;
                elseif exist('C_FIX_val','var') == 1
                    CF = C_FIX_val;
                    useCFIX = true;
                end
            end
        catch
            useCFIX = false;
        end
    end

    if useCFIX
        try
            rtmp = [xecef yecef zecef];   % Nx3
            rfix = rtmp * (CF.');         % Nx3
            xecef = rfix(:,1);
            yecef = rfix(:,2);
            zecef = rfix(:,3);
        catch
            fprintf('Aviso: no pude aplicar C_FIX/C_FIX_val. Continúo con ERA-only.\n');
        end
    end

    [lat_deg, lon_deg] = local_ecef2lla_wgs84_deg(xecef, yecef, zecef);
end

lon_deg = mod(lon_deg + 180, 360) - 180;
jump = [false; abs(diff(lon_deg)) > 180];
lat_deg(jump) = NaN;
lon_deg(jump) = NaN;

shadowOn = false(size(tpos));
if data.shadow.ok
    tsh = data.shadow.t(:);
    sh  = data.shadow.X(:) > 0.5;
    if numel(tsh) > 1
        shadowOn = interp1(tsh, double(sh), tpos, 'nearest', 'extrap') > 0.5;
    end
end

lat_sun = lat_deg; lon_sun = lon_deg;
lat_sh  = lat_deg; lon_sh  = lon_deg;
lat_sun(shadowOn) = NaN; lon_sun(shadowOn) = NaN;
lat_sh(~shadowOn) = NaN; lon_sh(~shadowOn) = NaN;

figure('Name','Trayectoria (Ground Track) - Sombra vs Iluminado','Color','w');
ax4 = axes; hold(ax4,'on'); grid(ax4,'on');

didMap = false;
try
    load coastlines
    plot(ax4, coastlon, coastlat, 'LineWidth', 0.8, 'HandleVisibility','off');
    didMap = true;
catch
    try
        load coast
        plot(ax4, long, lat, 'LineWidth', 0.8, 'HandleVisibility','off');
        didMap = true;
    catch
        didMap = false;
    end
end

hSun = plot(ax4, lon_sun, lat_sun, 'LineWidth', 1.6);
hSh  = plot(ax4, lon_sh,  lat_sh,  'LineWidth', 1.6, 'LineStyle','--');

xlabel(ax4,'Longitud [deg]');
ylabel(ax4,'Latitud [deg]');

if didMap
    title(ax4, sprintf('Ground track - pos de %s(:,1:3) en ECI', xNameUsed));
else
    title(ax4, sprintf('Ground track (sin mapa base) - pos de %s(:,1:3) en ECI', xNameUsed));
end

xlim(ax4,[-180 180]);
ylim(ax4,[-90 90]);
legend(ax4, [hSun hSh], {'Iluminado','Sombra'}, 'Location','best');

%% =========================
% 11) Figuras de actitud (si existe x_angular)
% =========================
if haveAtt
    figure('Name','Actitud - Euler ZYX (deg)','Color','w');
    axE = axes; hold(axE,'on'); grid(axE,'on');
    he1 = plot(axE, tA, yaw_deg,   'LineWidth', 1.3);
    he2 = plot(axE, tA, pitch_deg, 'LineWidth', 1.3);
    he3 = plot(axE, tA, roll_deg,  'LineWidth', 1.3);
    xlabel(axE,'Tiempo [s]');
    ylabel(axE,'Ángulo [deg]');
    title(axE, sprintf('Euler ZYX desde q (mode=%s) - q normalizado solo para DCM', quat_mode));
    local_apply_shadow(axE, shadowSegments);
    legend(axE, [he1 he2 he3], {'Yaw \psi','Pitch \theta','Roll \phi'}, 'Location','best');

    figure('Name','Actitud - Velocidades angulares (rad/s)','Color','w');
    axW = axes; hold(axW,'on'); grid(axW,'on');
    hw1 = plot(axW, tA, wB(:,1), 'LineWidth', 1.2);
    hw2 = plot(axW, tA, wB(:,2), 'LineWidth', 1.2);
    hw3 = plot(axW, tA, wB(:,3), 'LineWidth', 1.2);
    hw4 = plot(axW, tA, omega_mag, 'LineWidth', 1.6);
    xlabel(axW,'Tiempo [s]');
    ylabel(axW,'\omega [rad/s]');
    title(axW,'Rates (body) + norma');
    local_apply_shadow(axW, shadowSegments);
    legend(axW, [hw1 hw2 hw3 hw4], {'\omega_x','\omega_y','\omega_z','|\omega|'}, 'Location','best');

    figure('Name','Actitud - Error apuntamiento a nadir (deg)','Color','w');
    axP = axes; hold(axP,'on'); grid(axP,'on');
    plot(axP, tA, nadir_err, 'LineWidth', 1.6);
    xlabel(axP,'Tiempo [s]');
    ylabel(axP,'Error [deg]');
    title(axP, sprintf('Error apuntamiento a nadir (axis=[%.0f %.0f %.0f]_B, mode=%s)', ...
        body_axis_for_nadir(1), body_axis_for_nadir(2), body_axis_for_nadir(3), quat_mode));
    local_apply_shadow(axP, shadowSegments);
end

%% =========================
% 12) Diagnósticos orbitales (a,e,i,RAAN)
% =========================
if PLOT_ORBIT_DIAGNOSTICS
    nT = size(rECI,1);
    a = nan(nT,1); e = nan(nT,1); inc = nan(nT,1);
    raan = nan(nT,1); argp = nan(nT,1); nu = nan(nT,1);

    for ii = 1:nT
        [a(ii), e(ii), inc(ii), raan(ii), argp(ii), nu(ii)] = local_coe_from_rv(rECI(ii,:).', vECI(ii,:).', muEarth);
    end

    figure('Name','Órbita - Semieje mayor','Color','w');
    axO = axes; hold(axO,'on'); grid(axO,'on');
    plot(axO, tx, a/1e3, 'LineWidth', 1.2);
    xlabel(axO,'Tiempo [s]'); ylabel(axO,'a [km]');
    title(axO,'Semieje mayor');
    local_apply_shadow(axO, shadowSegments);

    figure('Name','Órbita - Excentricidad','Color','w');
    axEcc = axes; hold(axEcc,'on'); grid(axEcc,'on');
    plot(axEcc, tx, e, 'LineWidth', 1.2);
    xlabel(axEcc,'Tiempo [s]'); ylabel(axEcc,'e [-]');
    title(axEcc,'Excentricidad');
    local_apply_shadow(axEcc, shadowSegments);

    figure('Name','Órbita - Inclinación','Color','w');
    axI = axes; hold(axI,'on'); grid(axI,'on');
    plot(axI, tx, inc*180/pi, 'LineWidth', 1.2);
    xlabel(axI,'Tiempo [s]'); ylabel(axI,'i [deg]');
    title(axI,'Inclinación');
    local_apply_shadow(axI, shadowSegments);

    figure('Name','Órbita - RAAN','Color','w');
    axR = axes; hold(axR,'on'); grid(axR,'on');
    plot(axR, tx, local_unwrap_deg(raan*180/pi), 'LineWidth', 1.2);
    xlabel(axR,'Tiempo [s]'); ylabel(axR,'\Omega [deg]');
    title(axR,'RAAN');
    local_apply_shadow(axR, shadowSegments);
end

%% =========================
% 13) Resumen de disponibilidad
% =========================
fprintf('\n=== RESUMEN DISPONIBILIDAD DE SEÑALES ===\n');
for k = 1:numel(sigKeys)
    key = sigKeys{k};
    if data.(key).ok
        fprintf('[OK]   %s\n', data.(key).name);
    else
        fprintf('[MISS] %s\n', data.(key).name);
    end
end

if data.F_grav_nonrad.ok
    fprintf('[OK]   %s\n', data.F_grav_nonrad.name);
else
    fprintf('[MISS] %s (no se pudo derivar)\n', data.F_grav_nonrad.name);
end

fprintf('[OK]   Estado usado (órbita):  %s\n', SIG.x_lineal);

if haveAtt
    fprintf('[OK]   Estado usado (actitud): %s\n', SIG.x_angular);
    fprintf('[INFO] Max |omega| = %.6g rad/s\n', max(omega_mag));
    fprintf('[INFO] Max (||q_raw||-1) = %.6g\n', max(abs(q_norm_err)));
    fprintf('[INFO] Max nadir err = %.6g deg\n', max(nadir_err));
else
    fprintf('[MISS] Estado de actitud: %s\n', SIG.x_angular);
end

%% =========================
% 14) COMPARACIÓN CONTRA GMAT (si está configurado)
% =========================
doGMAT = ~isempty(GMAT_FILE) && exist(GMAT_FILE,'file') == 2;

need = {'t','rx','ry','rz','vx','vy','vz'};
missingMap = false;
for i = 1:numel(need)
    if ~isfield(GMAT_COLMAP, need{i}) || isempty(GMAT_COLMAP.(need{i}))
        missingMap = true;
    end
end

if doGMAT && ~missingMap
    fprintf('\n=== GMAT: leyendo %s ===\n', GMAT_FILE);
    T = local_read_table_file(GMAT_FILE);

    % Validar columnas requeridas
    local_require_cols(T, GMAT_COLMAP);

    tG = T.(GMAT_COLMAP.t)(:);
    rG = [T.(GMAT_COLMAP.rx)(:), T.(GMAT_COLMAP.ry)(:), T.(GMAT_COLMAP.rz)(:)];
    vG = [T.(GMAT_COLMAP.vx)(:), T.(GMAT_COLMAP.vy)(:), T.(GMAT_COLMAP.vz)(:)];

    % Interpolar simulador a tiempos GMAT
    rS = [ ...
        interp1(tx, rECI(:,1), tG, 'linear', 'extrap'), ...
        interp1(tx, rECI(:,2), tG, 'linear', 'extrap'), ...
        interp1(tx, rECI(:,3), tG, 'linear', 'extrap')];

    vS = [ ...
        interp1(tx, vECI(:,1), tG, 'linear', 'extrap'), ...
        interp1(tx, vECI(:,2), tG, 'linear', 'extrap'), ...
        interp1(tx, vECI(:,3), tG, 'linear', 'extrap')];

    dr = rS - rG;
    dv = vS - vG;

    drn = sqrt(sum(dr.^2,2));
    dvn = sqrt(sum(dv.^2,2));

    figure('Name','GMAT vs Sim - ||Δr||','Color','w');
    ax = axes; hold(ax,'on'); grid(ax,'on');
    plot(ax, tG, drn, 'LineWidth', 1.4);
    xlabel(ax,'Tiempo GMAT [s]');
    ylabel(ax,'||r_{sim} - r_{GMAT}|| [m]');
    title(ax,'Error posición (norma)');

    figure('Name','GMAT vs Sim - ||Δv||','Color','w');
    ax = axes; hold(ax,'on'); grid(ax,'on');
    plot(ax, tG, dvn, 'LineWidth', 1.4);
    xlabel(ax,'Tiempo GMAT [s]');
    ylabel(ax,'||v_{sim} - v_{GMAT}|| [m/s]');
    title(ax,'Error velocidad (norma)');

    fprintf('[GMAT] max||Δr|| = %.6g m | RMS = %.6g m\n', max(drn), sqrt(mean(drn.^2)));
    fprintf('[GMAT] max||Δv|| = %.6g m/s | RMS = %.6g m/s\n', max(dvn), sqrt(mean(dvn.^2)));

    if isfield(GMAT_COLMAP,'rho') && ~isempty(GMAT_COLMAP.rho) && data.rho.ok
        if ismember(GMAT_COLMAP.rho, T.Properties.VariableNames)
            rhoG = T.(GMAT_COLMAP.rho)(:);
            rhoS = interp1(data.rho.t(:), data.rho.X(:), tG, 'linear', 'extrap');

            figure('Name','GMAT vs Sim - rho','Color','w');
            ax = axes; hold(ax,'on'); grid(ax,'on');
            plot(ax, tG, rhoS, 'LineWidth', 1.2);
            plot(ax, tG, rhoG, 'LineWidth', 1.2);
            xlabel(ax,'Tiempo GMAT [s]');
            ylabel(ax,'\rho [kg/m^3]');
            title(ax,'Densidad atmosférica');
            legend(ax, {'Sim','GMAT'}, 'Location','best');
            set(ax,'YScale','log');
        end
    end

elseif doGMAT && missingMap
    fprintf(['\n[GMAT] Archivo definido pero falta mapear columnas en GMAT_COLMAP.' nl ...
             'No ejecuto comparación (sin inferencias).' nl]);
elseif ~doGMAT
    fprintf('\n[GMAT] No hay comparación (GMAT_FILE vacío o inexistente).\n');
end

fprintf('\nListo.\n');

%% =========================
% Local functions
% =========================
function ts = local_get_dataset_timeseries_exact(logs, name)
    ts = [];
    if isempty(logs), return; end

    try
        elem = logs.get(name);
        if ~isempty(elem)
            ts = elem.Values;
            return;
        end
    catch
    end

    try
        nEl = logs.numElements;
        for i = 1:nEl
            nm = '';
            try, nm = logs{i}.Name; catch, end
            if ~isempty(nm) && strcmpi(nm, name)
                try
                    ts = logs{i}.Values;
                    return;
                catch
                    ts = [];
                    return;
                end
            end
        end
    catch
    end
end

function ts = local_get_sdi_timeseries(sigObj)
    ts = [];
    try, ts = sigObj.Values; return; catch, end
    try, ts = sigObj.getValues; return; catch, end
end

function ts = local_get_sdi_timeseries_by_name(sdiNames, sdiTS, name)
    ts = [];
    if isempty(sdiNames), return; end

    idx = find(strcmp(sdiNames, name), 1, 'first');
    if ~isempty(idx) && ~isempty(sdiTS{idx})
        ts = sdiTS{idx}; return;
    end

    idx = find(strcmpi(sdiNames, name), 1, 'first');
    if ~isempty(idx) && ~isempty(sdiTS{idx})
        ts = sdiTS{idx}; return;
    end
end

function [t, X] = local_read_ts(ts)
    t = ts.Time(:);
    Xraw = ts.Data;
    Xraw = squeeze(Xraw);

    if size(Xraw,1) == numel(t)
        X = Xraw;
    elseif size(Xraw,2) == numel(t)
        X = Xraw.';
    else
        X = Xraw;
    end

    if isrow(X), X = X(:); end
end

function mag = local_mag_timeseries_samplewise(X, t)
% Devuelve magnitud Nx1 por muestra, incluso si Data tiene dims extra.
% Evita que plot() genere múltiples líneas (y que explote h(end+1)=...).
    n = numel(t);

    Xs = squeeze(X);

    % Alinear dimensión temporal a filas
    if size(Xs,1) ~= n && size(Xs,2) == n
        Xs = permute(Xs, [2 1 3:ndims(Xs)]);
    end

    % Forzar N x M
    if size(Xs,1) ~= n
        try
            Xs = reshape(Xs, n, []);
        catch
            error('No puedo alinear ts.Data con t: size(Data)=%s, numel(t)=%d', mat2str(size(X)), n);
        end
    else
        Xs = reshape(Xs, n, []);
    end

    mag = sqrt(sum(Xs.^2, 2));
end

function local_apply_shadow(ax, segs)
    if isempty(segs), return; end
    yl = get(ax,'YLim');

    for r = 1:size(segs,1)
        patch(ax, ...
            [segs(r,1) segs(r,2) segs(r,2) segs(r,1)], ...
            [yl(1) yl(1) yl(2) yl(2)], ...
            [0.85 0.85 0.85], ...
            'EdgeColor','none', 'FaceAlpha',0.25, 'HandleVisibility','off');
    end

    uistack(findall(ax,'Type','line'),'top');
end

function [lat_deg, lon_deg] = local_ecef2lla_wgs84_deg(x, y, z)
    a  = 6378137.0;
    f  = 1/298.257223563;
    e2 = f*(2-f);
    b  = a*sqrt(1-e2);
    ep2 = (a^2 - b^2)/b^2;

    lon = atan2(y, x);
    p   = sqrt(x.^2 + y.^2);

    theta = atan2(z*a, p*b);
    st = sin(theta);
    ct = cos(theta);

    lat = atan2(z + ep2*b*st.^3, p - e2*a*ct.^3);

    lat_deg = lat * 180/pi;
    lon_deg = lon * 180/pi;
end

function q = local_quat_continuity(q)
    for i = 2:size(q,1)
        if sum(q(i,:).*q(i-1,:)) < 0
            q(i,:) = -q(i,:);
        end
    end
end

function [c11,c12,c13,c21,c22,c23,c31,c32,c33] = local_quat_to_dcm_cib(q, mode)
    q0 = q(:,1); q1 = q(:,2); q2 = q(:,3); q3 = q(:,4);

    C11 = q0.^2 + q1.^2 - q2.^2 - q3.^2;
    C12 = 2*(q1.*q2 - q0.*q3);
    C13 = 2*(q1.*q3 + q0.*q2);

    C21 = 2*(q1.*q2 + q0.*q3);
    C22 = q0.^2 - q1.^2 + q2.^2 - q3.^2;
    C23 = 2*(q2.*q3 - q0.*q1);

    C31 = 2*(q1.*q3 - q0.*q2);
    C32 = 2*(q2.*q3 + q0.*q1);
    C33 = q0.^2 - q1.^2 - q2.^2 + q3.^2;

    if strcmpi(mode,'CIB')
        c11=C11; c12=C12; c13=C13;
        c21=C21; c22=C22; c23=C23;
        c31=C31; c32=C32; c33=C33;
    else
        c11=C11; c12=C21; c13=C31;
        c21=C12; c22=C22; c23=C32;
        c31=C13; c32=C23; c33=C33;
    end
end

function [yaw, pitch, roll] = local_dcm_to_euler_zyx(c11,c21,c31,c32,c33)
    s = -c31;
    s = max(min(s, 1), -1);
    pitch = asin(s);
    roll  = atan2(c32, c33);
    yaw   = atan2(c21, c11);
end

function a_deg = local_unwrap_deg(a_deg)
    for i = 2:numel(a_deg)
        da = a_deg(i) - a_deg(i-1);
        if da > 180
            a_deg(i:end) = a_deg(i:end) - 360;
        elseif da < -180
            a_deg(i:end) = a_deg(i:end) + 360;
        end
    end
end

function [a, e, inc, raan, argp, nu] = local_coe_from_rv(rI, vI, mu)
% COE (osculating) desde r,v inerciales.
% Devuelve ángulos en rad en [0,2*pi) cuando aplica.
% Convenciones en degenerados:
%  - Ecuatorial (inc ~ 0): RAAN=0 y argp = longitude of periapsis (varpi)
%  - Circular (e ~ 0):     argp=0 y nu = argument of latitude (u) si no ecuatorial;
%                          si además ecuatorial: nu = true longitude (ell)

    tol_n = 1e-12;
    tol_e = 1e-12;

    r = rI(:); v = vI(:);
    R = norm(r);
    V = norm(v);

    % Momento angular específico
    h = cross(r, v);
    hnorm = norm(h);

    % Vector nodo
    k = [0;0;1];
    n = cross(k, h);
    nnorm = norm(n);

    % Vector excentricidad
    evec = (1/mu)*((V^2 - mu/R)*r - (dot(r,v))*v);
    e = norm(evec);

    % Energía específica (2-body equivalente)
    energy = V^2/2 - mu/R;
    if abs(energy) > 1e-16
        a = -mu/(2*energy);
    else
        a = Inf;
    end

    % Inclinación (con clamp)
    if hnorm < tol_n
        inc = 0;
    else
        c = h(3)/hnorm;
        c = max(min(c,1),-1);
        inc = acos(c);
    end

    % RAAN
    if nnorm < tol_n
        raan = 0;
    else
        raan = atan2(n(2), n(1));
        if raan < 0, raan = raan + 2*pi; end
    end

    % Argumento del periapsis y anomalía verdadera
    if e < tol_e
        % Circular: periapsis indefinido
        argp = 0;

        if nnorm < tol_n
            % Circular ecuatorial: usar true longitude ell
            nu = atan2(r(2), r(1));
            if nu < 0, nu = nu + 2*pi; end
        else
            % Circular inclinada: usar argument of latitude u (desde nodo a r)
            % u = atan2( ( (n x r) · h ) / (|n||r||h|), (n · r)/(|n||r|) )
            sin_u = dot(cross(n, r), h) / (nnorm*R*hnorm);
            cos_u = dot(n, r) / (nnorm*R);
            nu = atan2(sin_u, cos_u);
            if nu < 0, nu = nu + 2*pi; end
        end

    else
        % No circular: nu con atan2 robusto
        % nu = atan2( ( (e x r) · h )/(|e||r||h|), (e · r)/(|e||r|) )
        sin_nu = dot(cross(evec, r), h) / (e*R*hnorm);
        cos_nu = dot(evec, r) / (e*R);
        nu = atan2(sin_nu, cos_nu);
        if nu < 0, nu = nu + 2*pi; end

        if nnorm < tol_n
            % Ecuatorial no circular: RAAN indefinido, usar longitude of periapsis varpi
            % varpi = atan2(e_y, e_x)
            argp = atan2(evec(2), evec(1));
            if argp < 0, argp = argp + 2*pi; end
            raan = 0;
        else
            % argp con atan2 robusto
            % argp = atan2( ( (n x e) · h )/(|n||e||h|), (n · e)/(|n||e|) )
            sin_w = dot(cross(n, evec), h) / (nnorm*e*hnorm);
            cos_w = dot(n, evec) / (nnorm*e);
            argp = atan2(sin_w, cos_w);
            if argp < 0, argp = argp + 2*pi; end
        end
    end
end

function theta = local_era_from_jd_ut1(jd_ut1_0, tsec)
    jd_ut1 = jd_ut1_0 + tsec(:)/86400;
    fraction = jd_ut1 - 2451545.0;
    rev = 0.7790572732640 + 1.00273781191135448 * fraction;
    rev_frac = rev - floor(rev);
    theta = 2*pi*rev_frac;
end

function T = local_read_table_file(fname)
% Lee CSV/TXT a tabla con header. No infiere mapeos.
    try
        % Si tu MATLAB soporta preservar nombres, esto minimiza sorpresas
        T = readtable(fname, 'VariableNamingRule','preserve');
    catch
        try
            T = readtable(fname);
        catch
            try
                T = readtable(fname, 'Delimiter', ',');
            catch ME
                error('No pude leer GMAT_FILE con readtable: %s', ME.message);
            end
        end
    end
end

function local_require_cols(T, MAP)
    nl = sprintf('\n');
    vars = T.Properties.VariableNames;

    req = {'t','rx','ry','rz','vx','vy','vz'};
    for i = 1:numel(req)
        nm = MAP.(req{i});
        if ~ismember(nm, vars)
            error(['GMAT: no encuentro la columna "%s".' nl ...
                   'Columnas disponibles:' nl '%s'], nm, strjoin(vars, ', '));
        end
    end

    if isfield(MAP,'rho') && ~isempty(MAP.rho)
        if ~ismember(MAP.rho, vars)
            error(['GMAT: pediste rho="%s" pero no existe en el archivo.' nl ...
                   'Columnas disponibles:' nl '%s'], MAP.rho, strjoin(vars, ', '));
        end
    end
end
