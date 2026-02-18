%% TEST_validate_vs_AeroTB.m
clear; clc;

INIT_parametros;
DATA_earth;

%% Config
cfg.Ncases     = 50;
cfg.seed       = 1;
cfg.h_ref_m    = 700e3;
cfg.h_span_m   = 100e3;

cfg.plot_xform = true;        % false para no plottear
cfg.t_span_sec = 30*24*3600;  % 1 mes

% tolerancias
tol.vec_abs = 1e-10;
tol.pos_m   = 1e-2;
tol.geo_mix = 1e-6;
tol.g_mps2  = 1e-8;

rng(cfg.seed);

fprintf('\n=== VALIDACION AUTOMATICA vs AERO TB ===\n');
fprintf('Casos: %d\n', cfg.Ncases);

%% Elipsoide
world.a   = a_world;
world.b   = b_world;
world.e2  = e2_world;
world.es2 = es2_world;
world.f   = (world.a - world.b)/world.a;

%% Casos LLA / quats
N = cfg.Ncases;

lat_deg = -80 + 160*rand(N,1);
lon_deg = -180 + 360*rand(N,1);
h_m     = cfg.h_ref_m + cfg.h_span_m*(2*rand(N,1)-1);
lla_deg_m = [lat_deg, lon_deg, h_m];

% quats random (scalar-first) y normalizados
Q = randn(4,N);
Q = Q ./ vecnorm(Q);

results = struct('test',{},'status',{},'maxErr',{},'rmsErr',{},'note',{});

%% ============================
% A: MT_ECI2BODY vs quat2dcm (AeroTB)
% ============================
try
    name = "A: MT_ECI2BODY vs quat2dcm (AeroTB)";
    errs = zeros(N,1);

    for k=1:N
        q  = Q(:,k);
        vI = randn(3,1);

        % AeroTB: quat2dcm devuelve C_BI (I->B) en conv aero estándar
        C_BI = quat2dcm(q.');   % q row [q0 q1 q2 q3]
        vB_ref = C_BI * vI;

        vB_u = MT_ECI2BODY(q, vI);
        errs(k) = norm(vB_u - vB_ref);
    end

    results(end+1) = packResult(name, errs, tol.vec_abs, "");
catch ME
    results(end+1) = packError(name, ME);
end

%% ============================
% B: MT_BODY2ECI vs quat2dcm (AeroTB)
% ============================
try
    name = "B: MT_BODY2ECI vs quat2dcm (AeroTB)";
    errs = zeros(N,1);

    for k=1:N
        q  = Q(:,k);
        vB = randn(3,1);

        C_BI = quat2dcm(q.');
        C_IB = C_BI.';          % B->I
        vI_ref = C_IB * vB;

        vI_u = MT_BODY2ECI(q, vB);
        errs(k) = norm(vI_u - vI_ref);
    end

    results(end+1) = packResult(name, errs, tol.vec_abs, "");
catch ME
    results(end+1) = packError(name, ME);
end

%% ============================
% C: MT_ECEF2GEODETIC vs ecef2lla
% Métrica: sqrt(dN^2+dE^2+dh^2) [m]
% ============================
try
    name = "C: MT_ECEF2GEODETIC vs ecef2lla (mismo elipsoide)";
    errs = zeros(N,1);

    for k=1:N
        r_ecef = lla2ecef(lla_deg_m(k,:), world.f, world.a).';
        r_ecef = r_ecef(:);

        [lat_u_deg, lon_u_deg, h_u_m] = MT_ECEF2GEODETIC(r_ecef, world.a, world.b, world.e2, world.es2);

        lla_ref = ecef2lla(r_ecef.', world.f, world.a);
        lat_ref = lla_ref(1);
        lon_ref = lla_ref(2);
        h_ref   = lla_ref(3);

        dlat_rad = deg2rad(lat_u_deg - lat_ref);
        dlon_rad = deg2rad(wrapTo180Local(lon_u_deg - lon_ref));
        dh       = h_u_m - h_ref;

        phi  = deg2rad(lat_ref);
        sinp = sin(phi);
        denom = sqrt(1 - world.e2*sinp*sinp);
        RN = world.a / denom;
        RM = world.a*(1-world.e2) / (denom^3);

        dN = dlat_rad * (RM + h_ref);
        dE = dlon_rad * (RN + h_ref) * cos(phi);

        errs(k) = sqrt(dN*dN + dE*dE + dh*dh);
    end

    results(end+1) = packResult(name, errs, tol.pos_m, "Métrica: sqrt(dN^2+dE^2+dh^2) [m].");
catch ME
    results(end+1) = packError(name, ME);
end

%% ============================
% D: MT_NED2ECEF vs AeroTB (DCM ECEF<->NED)
% ============================
try
    name = "D: MT_NED2ECEF vs AeroTB dcmecef2ned (nT->T)";
    errs = zeros(N,1);

    for k = 1:N
        phi_deg = lla_deg_m(k,1);
        lam_deg = lla_deg_m(k,2);

        % Vector en NED [nT]
        B_ned_nT = 50000*randn(3,1);

        % Tu implementación (entrada [nT], salida [T])
        B_ecef_u_T = MT_NED2ECEF(B_ned_nT, phi_deg, lam_deg);

        % Referencia AeroTB: DCM ECEF->NED, entonces NED->ECEF = DCM'
        phi = deg2rad(phi_deg);
        lam = deg2rad(lam_deg);

        C_NE = dcmecef2ned(phi_deg, lam_deg);      % ECEF -> NED
        C_EN = C_NE.';                     % NED  -> ECEF

        B_ecef_ref_T = C_EN * (B_ned_nT * 1e-9);

        errs(k) = norm(B_ecef_u_T - B_ecef_ref_T);
    end

    results(end+1) = packResult(name, errs, 1e-12, "Ref: dcmecef2ned. Entrada [nT] salida [T].");
catch ME
    results(end+1) = packError(name, ME);
end


%% ============================
% E: MT_ECEF2BODY_matrix vs referencia (quat2dcm + tu XFORM)
% ============================
try
    name = "E: MT_ECEF2BODY_matrix vs referencia (quat2dcm + XFORM)";
    errs = zeros(N,1);

    % XFORM_MODE robusto
    xm = 1;
    if exist('XFORM_MODE','var')
        if isa(XFORM_MODE,'Simulink.Parameter'), xm = double(XFORM_MODE.Value);
        else, xm = double(XFORM_MODE);
        end
    end

    % C_FIX_val robusto
    if exist('C_FIX_val','var')==0
        C_FIX_num = eye(3);
    elseif isa(C_FIX_val,'Simulink.Parameter')
        C_FIX_num = double(C_FIX_val.Value);
    else
        C_FIX_num = double(C_FIX_val);
    end

    for k=1:N
        q   = Q(:,k);
        era = 2*pi*rand();

        A = randn(3,3);
        J_ecef = A'*A + 1e-3*eye(3);

        % función a validar (tu implementación)
        J_body_u = MT_ECEF2BODY_matrix(q, era, J_ecef, C_FIX_num);

        % ---- referencia ----
        % C_BI desde AeroTB
        C_BI = quat2dcm(q.');

        % ECEF->ECI coherente con tu MT_ECEF2ECI
        c = cos(era); s = sin(era);
        R3p = [ c, -s, 0;
                s,  c, 0;
                0,  0, 1];

        if xm == 2
            % tu convención actual: ECI->ECEF = R3(-ERA) * Cfix
            % => ECEF->ECI = (ECI->ECEF)' = Cfix' * R3(+ERA) = C_FIX' * R3p
            C_IE = (C_FIX_num.') * R3p;
        else
            C_IE = R3p;
        end

        % ECEF->BODY = (ECI->BODY)*(ECEF->ECI)
        C_BE = C_BI * C_IE;

        J_body_ref = C_BE * J_ecef * C_BE.';
        errs(k) = norm(J_body_u - J_body_ref, 'fro');
    end

    results(end+1) = packResult(name, errs, 1e-9, "Valida ECEF->BODY consistente con tu XFORM.");
catch ME
    results(end+1) = packError(name, ME);
end

%% ============================
% F: MT_ECEF2GEOCENTRIC vs cart2sph (equivalente en m)
% ============================
try
    name = "F: MT_ECEF2GEOCENTRIC vs cart2sph (equivalente en m)";
    errs = zeros(N,1);

    for k=1:N
        r_ecef = 1e6*randn(3,1);

        [r_u, lam_u, ~, ~, ~, ~, phi_u] = MT_ECEF2GEOCENTRIC(r_ecef);
        [az, el, r_ref] = cart2sph(r_ecef(1), r_ecef(2), r_ecef(3));

        dr   = (r_u - r_ref);
        dlam = wrapToPiLocal(lam_u - az);
        dphi = (phi_u - el);

        dang_m = r_ref * sqrt( (dphi)^2 + (cos(el)*dlam)^2 );  % el ~= phi geocéntrica
        errs(k) = sqrt(dr^2 + dang_m^2);
    end

    results(end+1) = packResult(name, errs, tol.geo_mix, "sqrt(dr^2 + (r*dang)^2) [m].");
catch ME
    results(end+1) = packError(name, ME);
end
%% ============================
% G: MT_GEOCENTRIC2ECEF vs referencia AeroTB (angle2dcm)
% Ref: angle2dcm (AeroTB) + transpose (convención pasiva) + permutación de ejes
% Base de tu función: v_geo = [v_r; v_phi; v_lam] con:
%   rhat (radial/out), phihat (north), lamhat (east)
% ============================
try
    name = "G: MT_GEOCENTRIC2ECEF vs angle2dcm (AeroTB)";
    errs = zeros(N,1);

    for k=1:N
        phi = (-pi/2) + pi*rand();      % [rad] geocéntrica
        lam = (-pi)   + 2*pi*rand();    % [rad]
        v_geo = randn(3,1);             % [vr; vphi; vlam]

        v_u = MT_GEOCENTRIC2ECEF(v_geo, phi, lam);

        % AeroTB: angle2dcm(psi,theta,phi,'ZYX') ~ DCM pasiva (inertial->body).
        % Con psi=lam, theta=-phi, roll=0:
        % - La matriz ACTIVA geoX->ECEF sería Rz(lam)*Ry(-phi), cuyas columnas son [rhat lamhat phihat].
        % - Si angle2dcm devuelve la PASIVA, esa es la transpuesta de la activa. Por eso invertimos con (').
        C_geoX_ecef = angle2dcm(lam, -phi, 0, 'ZYX');  % pasiva (ECEF->geoX)
        C_ecef_geoX = C_geoX_ecef.';                   % activa (geoX->ECEF) = [rhat lamhat phihat]

        % Tu base geocéntrica es [rhat phihat lamhat] (no [rhat lamhat phihat])
        C_ecef_geo = C_ecef_geoX(:, [1 3 2]);          % [rhat phihat lamhat]

        v_ref = C_ecef_geo * v_geo;

        errs(k) = norm(v_u - v_ref);
    end

    results(end+1) = packResult(name, errs, 1e-12, "Ref: angle2dcm (AeroTB) usando transpose (pasiva->activa) + permutación.");
catch ME
    results(end+1) = packError(name, ME);
end



%% ============================
% H: MD_aceleracion_geopotencial vs gravitysphericalharmonic (Custom, deg=3)
% ============================
try
    name = "H: MD_aceleracion_geopotencial vs gravitysphericalharmonic (Custom, deg=3)";

    deg = 3;
    Re  = 6378136.3;
    mu  = GM;

    C = Cbar; S = Sbar;
    C(1,1) = 1;
    degree = deg;
    GM = mu; 

    tmpFile = fullfile(tempdir, "__grav_custom_deg3.mat");
    save(tmpFile, "Re", "GM", "degree", "C", "S");

    errs = zeros(N,1);
    for k=1:N
        r_ecef = lla2ecef(lla_deg_m(k,:), world.f, world.a).';
        r_ecef = r_ecef(:);

        [r, lambda, sinphi, cosphi, s, rhoP, phi] = MT_ECEF2GEOCENTRIC(r_ecef);
        a_u = MD_aceleracion_geopotencial(r, lambda, sinphi, cosphi, s, rhoP, C, S, K, mu, phi);

        [gx,gy,gz] = gravitysphericalharmonic(r_ecef.', "Custom", deg, {tmpFile, @load}, "None");
        a_tb = [gx;gy;gz];

        errs(k) = norm(a_u(:) - a_tb(:));
    end

    results(end+1) = packResult(name, errs, tol.g_mps2, "Deg=3 (4x4). Custom MAT con GM.");
catch ME
    results(end+1) = packError(name, ME);
end

%% ============================
% I. XFORM: ECI->ECEF vs AeroTB (1 mes, multi-r, EOP t0 fijo)
% - tiempos ORDENADOS
% - varios vectores r_eci por instante (set fijo para evitar serrucho por muestreo)
% ============================
try
    name = "I. XFORM: ECI->ECEF vs AeroTB (1 mes, multi-r, EOP t0)";

    if exist('utc0','var')==0,          error("No existe utc0."); end
    if exist('jd_ut_inicial','var')==0, error("No existe jd_ut_inicial."); end
    if exist('deltaAT','var')==0, deltaAT = 37; end
    if exist('deltaUT1','var')==0 || exist('xp_rad','var')==0 || exist('yp_rad','var')==0
        error("Faltan deltaUT1/xp_rad/yp_rad en workspace (INIT_parametros).");
    end

    % C_FIX_val robusto
    if exist('C_FIX_val','var')==0
        C_FIX_num = eye(3);
    elseif isa(C_FIX_val,'Simulink.Parameter')
        C_FIX_num = double(C_FIX_val.Value);
    else
        C_FIX_num = double(C_FIX_val);
    end

    utc0_dt = utc_to_datetime_robust_XFORM(utc0);

    % EOP fijo t0
    deltaUT1_0 = double(deltaUT1);
    polarmotion0 = [double(xp_rad) double(yp_rad)];

    % experimento
    Ntime = N;
    Nvec  = 200;                 % más grande = max más estable
    t_sec = linspace(0, cfg.t_span_sec, Ntime).';  % ordenado

    rmag = world.a + cfg.h_ref_m;

    % set fijo de direcciones (evita serrucho por muestreo)
    U = randn(3,Nvec);
    U = U ./ vecnorm(U);

    err_trad_mean = zeros(Ntime,1);
    err_fix_mean  = zeros(Ntime,1);
    err_trad_max  = zeros(Ntime,1);
    err_fix_max   = zeros(Ntime,1);

    for k = 1:Ntime
        tk = t_sec(k);

        utc_t = utc0_dt + seconds(tk);
        C_full = dcmeci2ecef('IAU-2000/2006', utc_t, deltaAT, deltaUT1_0, polarmotion0);

        [ERA_t, ~] = MD_era(tk, jd_ut_inicial);

        c = cos(ERA_t); s = sin(ERA_t);
        R3m = [ c,  s, 0;
               -s,  c, 0;
                0,  0, 1];   % ECI->ECEF (ERA-only)

        C_trad = R3m;

        % tu convención actual: ECI->ECEF = R3(-ERA) * Cfix
        C_fix  = R3m * C_FIX_num;

        errs_trad = zeros(Nvec,1);
        errs_fix  = zeros(Nvec,1);

        for j = 1:Nvec
            r_eci = rmag * U(:,j);

            r_tb   = C_full * r_eci;
            r_trad = C_trad * r_eci;
            r_fix  = C_fix  * r_eci;

            errs_trad(j) = norm(r_trad - r_tb);
            errs_fix(j)  = norm(r_fix  - r_tb);
        end

        err_trad_mean(k) = mean(errs_trad);
        err_fix_mean(k)  = mean(errs_fix);
        err_trad_max(k)  = max(errs_trad);
        err_fix_max(k)   = max(errs_fix);
    end

    mxTrad  = max(err_trad_max);
    mxFix   = max(err_fix_max);
    rmsTrad = sqrt(mean(err_trad_mean.^2));
    rmsFix  = sqrt(mean(err_fix_mean.^2));

    note = sprintf("Trad: max=%.3e m rms(mean)=%.3e | Fix: max=%.3e m rms(mean)=%.3e", ...
                   Nvec, mxTrad, rmsTrad, mxFix, rmsFix);

    results(end+1) = struct('test',string(name),'status',"INFO", ...
                            'maxErr',mxFix,'rmsErr',rmsFix,'note',string(note));

    if cfg.plot_xform
        t_days = t_sec/86400;
        figure('Name','XFORM error vs time (1 mes, multi-r fijo)');
        plot(t_days, err_trad_mean); hold on;
        plot(t_days, err_fix_mean);
        plot(t_days, err_trad_max);
        plot(t_days, err_fix_max);
        grid on;
        xlabel('t [days]');
        ylabel('|| r_{ecef}^{(model)} - r_{ecef}^{(AeroTB)} || [m]');
        title('Error ECI->ECEF vs tiempo (EOP t0 fijo, multi-r fijo)');
        legend('ERA-only mean','FIXED\_CORR mean','ERA-only max','FIXED\_CORR max', ...
               'Location','best');
    end

catch ME
    results(end+1) = packError(name, ME);
end

%% ============================
% Reporte
% ============================
fprintf('\n--- REPORTE ---\n');
for i=1:numel(results)
    r = results(i);
    fprintf('%-60s  %-6s', r.test, r.status);
    if isfinite(r.maxErr)
        fprintf('  max=% .3e  rms=% .3e', r.maxErr, r.rmsErr);
    end
    if strlength(r.note) > 0
        fprintf('  | %s', r.note);
    end
    fprintf('\n');
end
fprintf('=== FIN ===\n');

%% ============================
% Helpers
% ============================

function out = packResult(name, errs, tol, note)
mx = max(errs);
rmsv = sqrt(mean(errs.^2));
status = "PASS";
if mx > tol
    status = "FAIL";
end
out = struct('test',string(name),'status',string(status),'maxErr',mx,'rmsErr',rmsv,'note',string(note));
end

function out = packError(name, ME)
out = struct('test',string(name),'status',"ERROR",'maxErr',nan,'rmsErr',nan,'note',string(ME.message));
end

function ang = wrapTo180Local(angdeg)
ang = mod(angdeg + 180, 360) - 180;
end

function ang = wrapToPiLocal(angrad)
ang = mod(angrad + pi, 2*pi) - pi;
end

function dt = utc_to_datetime_robust_XFORM(x)
% Acepta:
% - datetime
% - string/char con fecha
% - datevec [Y M D h m s]
% - datenum (≈ 7e5 para años modernos)
% - juliandate (≈ 2.4e6)
% - MJD (≈ 6e4) -> JD = MJD + 2400000.5

    if isa(x,'datetime')
        dt = x;
        if isempty(dt.TimeZone), dt.TimeZone = 'UTC'; end
        return
    end

    if isstring(x) || ischar(x)
        dt = datetime(x,'TimeZone','UTC');
        return
    end

    if isnumeric(x)
        x = double(x);

        if isvector(x) && numel(x)==6
            dt = datetime(x(1),x(2),x(3),x(4),x(5),x(6),'TimeZone','UTC');
            return
        end

        if isscalar(x)
            if x > 1e9
                dt = datetime(x,'ConvertFrom','posixtime','TimeZone','UTC');
            elseif x > 2e6
                dt = datetime(x,'ConvertFrom','juliandate','TimeZone','UTC');
            elseif x > 5e5
                dt = datetime(x,'ConvertFrom','datenum','TimeZone','UTC');
            elseif x > 3e4
                dt = datetime(x + 2400000.5,'ConvertFrom','juliandate','TimeZone','UTC');
            else
                error("utc0 numérico demasiado chico/ambiguo (%.3g). Pasalo como datetime o datenum real.", x);
            end

            if year(dt) < 1
                error("utc0 convertido a año inválido (%d). utc0 está mal formateado.", year(dt));
            end
            return
        end
    end

    error("utc0 debe ser datetime, string, datevec[6], datenum, juliandate o MJD.");
end
