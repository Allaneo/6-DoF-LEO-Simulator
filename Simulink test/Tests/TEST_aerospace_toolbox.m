%% TEST_validate_vs_AeroTB.m
% Valida funciones propias vs Aerospace Toolbox 

clear; clc;
INIT_parametros;

%% Config
cfg.Ncases   = 50;
cfg.seed     = 1;
cfg.h_ref_m  = 700e3;
cfg.h_span_m = 100e3;

% Tolerancias
tol.vec_abs  = 1e-10;
tol.dcm_fro  = 1e-12;
tol.pos_m    = 1e-2;     % 1 cm
tol.ang_deg  = 1e-7;     % ~1e-7 deg
tol.geo_mix  = 1e-9;     % para test geocéntrico mezclado
tol.g_mps2   = 1e-8;     % gravedad [m/s^2] (custom vs AeroTB)
tol.Ja_fro   = 1e-10;    % Jacobiano de gravedad (Frobenius)
tol.rho_rel  = 5e-3;     % 0.5% relativo (densidad: continuidad/knots)

rng(cfg.seed);

fprintf('\n=== VALIDACION AUTOMATICA vs AERO TB ===\n');
fprintf('Casos: %d\n', cfg.Ncases);

%% Parámetros del elipsoide (existen y punto)
world.a   = a_world;
world.b   = b_world;
world.e2  = e2_world;
world.es2 = es2_world;
world.f   = (world.a - world.b)/world.a;

%% Descubrir funciones
files = dir('*.m');
myNames = string(erase({files.name}, '.m'));
hasMy = @(name) any(myNames == string(name));
hasTB = @(name) exist(name,'file')==2 || exist(name,'file')==5;

%% Casos LLA y quaterniones
N = cfg.Ncases;
lat_deg = -80 + 160*rand(N,1);
lon_deg = -180 + 360*rand(N,1);
h_m     = cfg.h_ref_m + cfg.h_span_m*(2*rand(N,1)-1);
lla_deg_m = [lat_deg, lon_deg, h_m];

% Quaterniones aleatorios (scalar-first q=[q0;q1;q2;q3])
Q = randn(4,N);
Q = Q ./ vecnorm(Q);

results = struct('test',{},'status',{},'maxErr',{},'rmsErr',{},'note',{});

%% ============================
% TEST A: MT_ECI2BODY vs quat2dcm
% ============================
% Aerospace TB: quat2dcm(q) devuelve DCM pasiva I->B: vB = C_BI*vI
if hasMy('MT_ECI2BODY') && hasTB('quat2dcm')
    name = "MT_ECI2BODY vs quat2dcm";
    errs = zeros(N,1);

    for k=1:N
        q  = Q(:,k);
        vI = randn(3,1);

        C_BI = quat2dcm(q.');     % I -> B (pasiva)
        vB_ref = C_BI*vI;

        vB_u = MT_ECI2BODY(q, vI);
        errs(k) = norm(vB_u - vB_ref);
    end

    results(end+1) = packResult(name, errs, tol.vec_abs, "");
else
    results(end+1) = packSkip("MT_ECI2BODY vs quat2dcm", ...
        missingWhy(hasMy('MT_ECI2BODY'), hasTB('quat2dcm'), "MT_ECI2BODY","quat2dcm"));
end

%% ============================
% TEST B: MT_BODY2ECI vs quat2dcm
% ============================
% vI = C_IB*vB con C_IB = C_BI'
if hasMy('MT_BODY2ECI') && hasTB('quat2dcm')
    name = "MT_BODY2ECI vs quat2dcm";
    errs = zeros(N,1);

    for k=1:N
        q  = Q(:,k);
        vB = randn(3,1);

        C_BI = quat2dcm(q.');     % I -> B
        C_IB = C_BI.';            % B -> I
        vI_ref = C_IB*vB;

        vI_u = MT_BODY2ECI(q, vB);
        errs(k) = norm(vI_u - vI_ref);
    end

    results(end+1) = packResult(name, errs, tol.vec_abs, "");
else
    results(end+1) = packSkip("MT_BODY2ECI vs quat2dcm", ...
        missingWhy(hasMy('MT_BODY2ECI'), hasTB('quat2dcm'), "MT_BODY2ECI","quat2dcm"));
end

%% ============================
% TEST C: DCM (via MT_ECI2BODY) vs quat2dcm (fro)
% ============================
if hasMy('MT_ECI2BODY') && hasTB('quat2dcm')
    name = "DCM (via MT_ECI2BODY) vs quat2dcm (fro)";
    errs = zeros(N,1);

    for k=1:N
        q = Q(:,k);
        C_BI_ref = quat2dcm(q.');

        e1=[1;0;0]; e2=[0;1;0]; e3=[0;0;1];
        C_BI_u = [MT_ECI2BODY(q,e1), MT_ECI2BODY(q,e2), MT_ECI2BODY(q,e3)];

        errs(k) = norm(C_BI_u - C_BI_ref, 'fro');
    end

    results(end+1) = packResult(name, errs, tol.dcm_fro, "");
end

%% ============================
% TEST D: MT_ECEF2GEODETIC vs ecef2lla (mismo elipsoide)
% ============================
if hasMy('MT_ECEF2GEODETIC') && hasTB('ecef2lla') && hasTB('lla2ecef')
    name = "MT_ECEF2GEODETIC vs ecef2lla (mismo elipsoide)";
    errs = zeros(N,1);

    for k=1:N
        r_ecef = lla2ecef(lla_deg_m(k,:), world.f, world.a);
        r_ecef = r_ecef(:);

        [lat_u_deg, lon_u_deg, h_u_m] = MT_ECEF2GEODETIC(r_ecef, world.a, world.b, world.e2, world.es2);

        lla_ref = ecef2lla(r_ecef.', world.f, world.a);
        lat_ref = lla_ref(1);
        lon_ref = lla_ref(2);
        h_ref   = lla_ref(3);

        dlat = lat_u_deg - lat_ref;
        dlon = wrapTo180Local(lon_u_deg - lon_ref);
        dh   = h_u_m - h_ref;

        dh_deg_equiv = dh / 111320;
        errs(k) = norm([dlat; dlon; dh_deg_equiv]);
    end

    results(end+1) = packResult(name, errs, max(tol.ang_deg, tol.pos_m/111320), ...
        "Error en deg-equivalentes (alt escalada).");
else
    results(end+1) = packSkip("MT_ECEF2GEODETIC vs ecef2lla", ...
        missingWhy(hasMy('MT_ECEF2GEODETIC'), hasTB('ecef2lla') && hasTB('lla2ecef'), ...
        "MT_ECEF2GEODETIC","ecef2lla/lla2ecef"));
end

%% ============================
% TEST E: MT_NED2ECEF (vector) vs fórmula cerrada
% ============================
if hasMy('MT_NED2ECEF')
    name = "MT_NED2ECEF (vector) vs fórmula NED->ECEF (+ nT->T)";
    errs = zeros(N,1);

    for k=1:N
        phi_deg = lla_deg_m(k,1);
        lam_deg = lla_deg_m(k,2);

        B_ned_nT = 50000*randn(3,1);

        B_ecef_u_T = MT_NED2ECEF(B_ned_nT, phi_deg, lam_deg);

        phi = deg2rad(phi_deg);
        lam = deg2rad(lam_deg);

        C_ned_ecef = [ ...
            -sin(phi)*cos(lam),  -sin(phi)*sin(lam),   cos(phi); ...
            -sin(lam),            cos(lam),            0;        ...
            -cos(phi)*cos(lam),  -cos(phi)*sin(lam),  -sin(phi)];

        B_ecef_ref_T = (C_ned_ecef.') * (B_ned_nT * 1e-9);

        errs(k) = norm(B_ecef_u_T - B_ecef_ref_T);
    end

    results(end+1) = packResult(name, errs, 1e-12, "Rotación + conversión de unidades.");
else
    results(end+1) = packSkip("MT_NED2ECEF vs fórmula", "no existe MT_NED2ECEF.");
end

%% ============================
% TEST F: MT_ECI2ECEF + MT_ECEF2ECI (inversa interna)
% ============================
if hasMy('MT_ECI2ECEF') && hasMy('MT_ECEF2ECI')
    name = "MT_ECI2ECEF <-> MT_ECEF2ECI (consistencia interna)";
    errs = zeros(N,1);

    for k=1:N
        r_eci = 1e6*randn(3,1);
        era   = 2*pi*rand();

        r_ecef = MT_ECI2ECEF(r_eci, era);
        r_eci2 = MT_ECEF2ECI(r_ecef, era);

        errs(k) = norm(r_eci2 - r_eci);
    end

    results(end+1) = packResult(name, errs, 1e-9, "Tol numérica (m).");
else
    results(end+1) = packSkip("MT_ECI2ECEF <-> MT_ECEF2ECI", "faltan funciones.");
end

%% ============================
% TEST G: MT_ECEF2BODY_matrix vs referencia
% ============================
if hasMy('MT_ECEF2BODY_matrix') && hasTB('quat2dcm')
    name = "MT_ECEF2BODY_matrix vs referencia";
    errs = zeros(N,1);

    for k=1:N
        q   = Q(:,k);
        era = 2*pi*rand();

        A = randn(3,3);
        J_ecef = A'*A + 1e-3*eye(3);

        J_body_u = MT_ECEF2BODY_matrix(q, era, J_ecef);

        c = cos(era); s = sin(era);
        C_IE = [c, -s, 0;
                s,  c, 0;
                0,  0, 1];     % ECEF -> ECI

        C_BI = quat2dcm(q.');  % ECI -> Body
        C_BE = C_BI * C_IE;    % ECEF -> Body

        J_body_ref = C_BE * J_ecef * C_BE.';
        errs(k) = norm(J_body_u - J_body_ref, 'fro');
    end

    results(end+1) = packResult(name, errs, 1e-12, "Frobenius.");
else
    results(end+1) = packSkip("MT_ECEF2BODY_matrix vs referencia", "faltan funciones.");
end

%% ============================
% TEST H: MT_ECEF2GEOCENTRIC vs cart2sph
% ============================
if hasMy('MT_ECEF2GEOCENTRIC')
    name = "MT_ECEF2GEOCENTRIC vs cart2sph";
    errs = zeros(N,1);

    for k=1:N
        r_ecef = 1e6*randn(3,1);

        [r_u, lam_u, ~, ~, ~, ~, phi_u] = MT_ECEF2GEOCENTRIC(r_ecef);

        [az, el, r_ref] = cart2sph(r_ecef(1), r_ecef(2), r_ecef(3));

        errs(k) = norm([r_u - r_ref;
                        wrapToPiLocal(lam_u - az);
                        (phi_u - el)]);
    end

    results(end+1) = packResult(name, errs, tol.geo_mix, "Tol relajada (mezcla m/rad).");
else
    results(end+1) = packSkip("MT_ECEF2GEOCENTRIC vs cart2sph", "no existe MT_ECEF2GEOCENTRIC.");
end

%% ============================
% TEST I: MT_GEOCENTRIC2ECEF sanity
% ============================
if hasMy('MT_GEOCENTRIC2ECEF')
    name = "MT_GEOCENTRIC2ECEF sanity (base rhat/phihat/lambdahat)";
    errs = zeros(N,1);

    for k=1:N
        phi = (-pi/2) + pi*rand();
        lam = (-pi) + 2*pi*rand();
        v_geo = randn(3,1);

        v_u = MT_GEOCENTRIC2ECEF(v_geo, phi, lam);

        cphi = cos(phi); sphi = sin(phi);
        cl = cos(lam);   sl = sin(lam);

        rhat   = [ cphi*cl;  cphi*sl;  sphi ];
        phihat = [-sphi*cl; -sphi*sl;  cphi ];
        lamhat = [   -sl;       cl;      0  ];

        T = [rhat, phihat, lamhat];
        v_ref = T*v_geo;

        errs(k) = norm(v_u - v_ref);
    end

    results(end+1) = packResult(name, errs, 1e-12, "");
else
    results(end+1) = packSkip("MT_GEOCENTRIC2ECEF sanity", "no existe MT_GEOCENTRIC2ECEF.");
end

%% ============================
% TEST J: MT_area_proyectada vs referencia (quat2dcm correcto)
% ============================
if hasMy('MT_area_proyectada') && hasTB('quat2dcm')
    name = "MT_area_proyectada vs referencia (quat2dcm)";
    errs = zeros(N,1);

    for k=1:N
        q = Q(:,k);
        vI = randn(3,1);
        A_trans = abs(randn(3,1)) + 0.1;

        A_u = MT_area_proyectada(vI, q, A_trans);

        C_BI = quat2dcm(q.');   % I -> B
        vB = C_BI*vI;
        vhat = vB / norm(vB);

        A_ref = A_trans(1)*abs(vhat(1)) + A_trans(2)*abs(vhat(2)) + A_trans(3)*abs(vhat(3));
        errs(k) = abs(A_u - A_ref);
    end

    results(end+1) = packResult(name, errs, 1e-12, "");
else
    results(end+1) = packSkip("MT_area_proyectada vs referencia", "faltan funciones/toolbox.");
end

%% ============================
% TEST K (grado alineado): MD_aceleracion_geopotencial vs gravitysphericalharmonic('Custom')
% ============================
if hasMy('MD_aceleracion_geopotencial') && hasTB('gravitysphericalharmonic') && hasTB('lla2ecef')

    name = "MD_aceleracion_geopotencial vs gravitysphericalharmonic (Custom, degree=3)";

    deg = 3;
    C3  = Cbar(1:deg+1, 1:deg+1);
    S3  = Sbar(1:deg+1, 1:deg+1);
    K3  = K   (1:deg+1, 1:deg+1);

    Re  = 6378136.3;         % consistente con tu MT_ECEF2GEOCENTRIC
    GMc = GM;

    tmpFile = fullfile(tempdir, "__grav_custom_deg3.mat");
    degree = deg; 
    C = C3; S = S3; GM = GMc;
    save(tmpFile, "Re", "GM", "degree", "C", "S");

    errs = zeros(N,1);

    for k=1:N
        r_ecef = lla2ecef(lla_deg_m(k,:), world.f, world.a);
        r_ecef = r_ecef(:);

        [r, lambda, sinphi, cosphi, s, rhoP, phi] = MT_ECEF2GEOCENTRIC(r_ecef);
        a_u = MD_aceleracion_geopotencial(r, lambda, sinphi, cosphi, s, rhoP, C3, S3, K3, GMc, phi);

        [gx,gy,gz] = gravitysphericalharmonic(r_ecef.', "Custom", deg, {tmpFile, @load}, "None");
        a_tb = [gx;gy;gz];

        errs(k) = norm(a_u - a_tb);
    end

    results(end+1) = packResult(name, errs, 1e-8, "Degree alineado (3) y C/S truncados.");
else
    results(end+1) = packSkip("MD_aceleracion_geopotencial vs gravitysphericalharmonic (degree=3)", ...
        "faltan funciones/toolbox.");
end

%% ============================
% TEST L: MD_Ja_forward vs diferencias centradas (J_a)
% ============================
if hasMy('MD_Ja_forward') && hasMy('MD_aceleracion_geopotencial') && hasMy('MT_ECEF2GEOCENTRIC') && hasTB('lla2ecef')

    name = "MD_Ja_forward vs central-diff (Fro)";

    deltaJ = 1.0; % [m]
    errs = zeros(N,1);

    for k=1:N
        R = lla2ecef(lla_deg_m(k,:), world.f, world.a);
        R = R(:);

        [r0, lam0, sphi0, cphi0, s0, rhoP0, phi0] = MT_ECEF2GEOCENTRIC(R);
        a0 = MD_aceleracion_geopotencial(r0, lam0, sphi0, cphi0, s0, rhoP0, Cbar, Sbar, K, GM, phi0);

        J_fwd = MD_Ja_forward(R, Sbar, Cbar, K, deltaJ, GM, a0);

        ex=[1;0;0]; ey=[0;1;0]; ez=[0;0;1];
        J_cd = zeros(3,3);
        E = [ex ey ez];

        for i=1:3
            Rp = R + deltaJ*E(:,i);
            Rm = R - deltaJ*E(:,i);

            [rp, lamp, sphp, cphp, sp, rhoPp, phip] = MT_ECEF2GEOCENTRIC(Rp);
            ap = MD_aceleracion_geopotencial(rp, lamp, sphp, cphp, sp, rhoPp, Cbar, Sbar, K, GM, phip);

            [rm, lamm, sphm, cphm, sm, rhoPm, phim] = MT_ECEF2GEOCENTRIC(Rm);
            am = MD_aceleracion_geopotencial(rm, lamm, sphm, cphm, sm, rhoPm, Cbar, Sbar, K, GM, phim);

            J_cd(:,i) = (ap - am)/(2*deltaJ);
        end

        errs(k) = norm(J_fwd - J_cd, "fro");
    end

    results(end+1) = packResult(name, errs, tol.Ja_fro, "Forward vs central-diff.");
else
    results(end+1) = packSkip("MD_Ja_forward vs central-diff", "faltan funciones/toolbox.");
end

%% ============================
% TEST G0/G1: Diagnóstico de gravedad (central y J2)
% ============================
mu = GM;
Re = 6378136.3;

%% ----------------------------
% G0b) TU MD central-only vs analítico
% ----------------------------
if hasMy('MD_aceleracion_geopotencial') && hasMy('MT_ECEF2GEOCENTRIC') && hasTB('lla2ecef')
    name = "G0b: MD_aceleracion_geopotencial central-only vs analítico";

    C0   = zeros(4,4);
    S0   = zeros(4,4);
    Ksub = K(1:4,1:4);

    errs = zeros(N,1);
    for k=1:N
        r_ecef = lla2ecef(lla_deg_m(k,:), world.f, world.a);
        r_ecef = r_ecef(:);

        [r, lambda, sinphi, cosphi, s, rhoP, phi] = MT_ECEF2GEOCENTRIC(r_ecef);

        a_u  = MD_aceleracion_geopotencial(r, lambda, sinphi, cosphi, s, rhoP, C0, S0, Ksub, mu, phi);
        a_an = -mu * r_ecef / norm(r_ecef)^3;

        errs(k) = norm(a_u - a_an);
    end

    results(end+1) = packResult(name, errs, 1e-10, "");
else
    results(end+1) = packSkip("G0b: TU central-only vs analítico", "faltan funciones.");
end

%% ----------------------------
% G1) J2-only (solo C20): TU vs TB vs Analítico
% ----------------------------
if hasMy('MD_aceleracion_geopotencial') && hasMy('MT_ECEF2GEOCENTRIC') && hasTB('gravitysphericalharmonic') && hasTB('lla2ecef')
    baseName = "G1: J2-only (C20) triángulo";

    deg = 3;

    Cj2 = zeros(4,4);
    Sj2 = zeros(4,4);
    Cj2(1,1) = 1;
    Cj2(3,1) = Cbar(3,1);

    tmpfile = fullfile(tempdir,'__CustomJ2.mat');
    degree = deg; 
    C = Cj2; S = Sj2; GM = mu;
    save(tmpfile,'Re','GM','degree','C','S');

    Ksub = K(1:4,1:4);

    errs_u_tb  = zeros(N,1);
    errs_u_an  = zeros(N,1);
    errs_tb_an = zeros(N,1);

    for k=1:N
        r_ecef = lla2ecef(lla_deg_m(k,:), world.f, world.a);
        r_ecef = r_ecef(:);

        [r, lambda, sinphi, cosphi, s, rhoP, phi] = MT_ECEF2GEOCENTRIC(r_ecef);

        a_u = MD_aceleracion_geopotencial(r, lambda, sinphi, cosphi, s, rhoP, Cj2, Sj2, Ksub, mu, phi);
        a_u = a_u(:);

        [gx,gy,gz] = gravitysphericalharmonic(r_ecef.', 'Custom', deg, {tmpfile @load}, 'None');
        a_tb = [gx;gy;gz];

        a_an = accelJ2_fromC20norm(r_ecef, mu, Re, Cj2(3,1));

        errs_u_tb(k)  = norm(a_u  - a_tb);
        errs_u_an(k)  = norm(a_u  - a_an);
        errs_tb_an(k) = norm(a_tb - a_an);
    end

    results(end+1) = packResult(baseName+" | TU vs TB", errs_u_tb,  1e-8,  "");
    results(end+1) = packResult(baseName+" | TU vs AN", errs_u_an,  1e-8,  "");
    results(end+1) = packResult(baseName+" | TB vs AN", errs_tb_an, 1e-10, "");
else
    results(end+1) = packSkip("G1: J2-only triángulo", "faltan funciones/toolbox.");
end


%% ============================
% TEST M: MD_densidad (knots + continuidad + monotonicidad)
% ============================
if hasMy('MD_densidad')

    h0 = [0 25 30 40 50 60 70 80 90 100 110 120 130 140 150 180 200 250 300 350 400 450 500 600 700 800 900 1000];
    rho0 = [1.225e+0 3.899e-2 1.774e-2 3.972e-3 1.057e-3 3.206e-4 8.770e-5 1.905e-5 3.396e-6 5.297e-7 ...
            9.661e-8 2.438e-8 8.484e-9 3.845e-9 2.070e-9 5.464e-10 2.789e-10 7.248e-11 2.418e-11 9.518e-12 ...
            3.725e-12 1.585e-12 6.967e-13 1.454e-13 3.614e-14 1.170e-14 5.245e-15 3.019e-15];

    name = "MD_densidad: knots exactos";
    errs = zeros(numel(h0),1);
    for i=1:numel(h0)
        r = (6378 + h0(i))*1000;
        rho = MD_densidad([r;0;0]);
        errs(i) = abs(rho - rho0(i)) / rho0(i);
    end
    results(end+1) = packResult(name, errs, 1e-12, "Error relativo en puntos tabulados.");

    name = "MD_densidad: continuidad en fronteras";
    eps_km = 1e-3;
    errs = zeros(numel(h0)-1,1);
    for i=1:(numel(h0)-1)
        hb = h0(i+1);

        r1 = (6378 + (hb - eps_km))*1000;
        r2 = (6378 + (hb + eps_km))*1000;

        rho1 = MD_densidad([r1;0;0]);
        rho2 = MD_densidad([r2;0;0]);

        errs(i) = abs(rho2 - rho1) / max(rho1, rho2);
    end
    results(end+1) = packResult(name, errs, tol.rho_rel, "Salto relativo en fronteras.");

    name = "MD_densidad: monotonicidad";
    hs = linspace(0, 1000, 500);
    rho_s = zeros(size(hs));
    for i=1:numel(hs)
        r = (6378 + hs(i))*1000;
        rho_s(i) = MD_densidad([r;0;0]);
    end
    viol = max(rho_s(2:end) - rho_s(1:end-1), 0);
    errs = viol(:) ./ max(rho_s(1:end-1).', 1e-300);

    results(end+1) = packResult(name, errs, 0, "0 = sin violaciones.");
else
    results(end+1) = packSkip("MD_densidad tests", "no existe MD_densidad.");
end

%% ============================
% TEST N: MD_solar_force (shadow gating + magnitud/dirección)
% ============================
if hasMy('MD_solar_force')

    name = "MD_solar_force: shadow + dirección";

    errs = zeros(N,1);
    for k=1:N
        P  = 4.5e-6;                 % Pa
        Cr = 1.2;
        A  = 2.0;                    % m^2
        u  = randn(3,1); u = u/norm(u);

        F0 = MD_solar_force(P, A, Cr, u, 0);
        F1 = MD_solar_force(P, A, Cr, u, 1);

        F0_ref = -P*Cr*A*u;
        F1_ref = [0;0;0];

        errs(k) = norm(F0 - F0_ref) + norm(F1 - F1_ref);
    end

    results(end+1) = packResult(name, errs, 1e-15, "");
else
    results(end+1) = packSkip("MD_solar_force tests", "no existe MD_solar_force.");
end

%% ============================
% TEST O: MD_sun_position (periodicidad)
% ============================
if hasMy('MD_sun_position')

    name = "MD_sun_position: periodicidad 360deg";
    errs = zeros(N,1);

    for k=1:N
        M = 360*rand();
        w = 360*rand();

        L1 = MD_sun_position(M, w);
        L2 = MD_sun_position(M+360, w);

        errs(k) = abs((L2 - L1) - 360);
    end

    results(end+1) = packResult(name, errs, 1e-12, "");
else
    results(end+1) = packSkip("MD_sun_position tests", "no existe MD_sun_position.");
end
%% ============================
% TEST G3: Aislar término por término (n,m) que rompe vs Toolbox (degree=3)
% ============================
if hasMy('MD_aceleracion_geopotencial') && hasMy('MT_ECEF2GEOCENTRIC') && ...
   hasTB('gravitysphericalharmonic') && hasTB('lla2ecef')

    deg  = 3;
    Re   = 6378136.3;     % consistente con MT_ECEF2GEOCENTRIC
    mu   = GM;            % en tu script ya existe GM (lo venís usando arriba)
    Ksub = K(1:4,1:4);    % porque tu MD está hardcodeada hasta grado 3 (4x4)

    nm_list = [ ...
        1 0;
        1 1;
        2 0;
        2 1;
        2 2;
        3 0;
        3 1;
        3 2;
        3 3];

    % para ranking del peor término
    termName = strings(size(nm_list,1),1);
    termMax  = zeros(size(nm_list,1),1);

    for ii = 1:size(nm_list,1)
        n = nm_list(ii,1);
        m = nm_list(ii,2);

        % Construyo modelo con SOLO C00 y un (n,m)
        C = zeros(4,4);
        S = zeros(4,4);
        C(1,1) = 1;

        C(n+1,m+1) = Cbar(n+1,m+1);
        S(n+1,m+1) = Sbar(n+1,m+1);

        tmpfile = fullfile(tempdir, sprintf('__Custom_n%d_m%d.mat', n, m));
        degree = deg; GM = mu;
        save(tmpfile, 'Re', 'GM', 'degree', 'C', 'S');

        errs = zeros(N,1);

        for k = 1:N
            r_ecef = lla2ecef(lla_deg_m(k,:), world.f, world.a);
            r_ecef = r_ecef(:);

            [r, lambda, sinphi, cosphi, s, rhoP, phi] = MT_ECEF2GEOCENTRIC(r_ecef);

            a_u = MD_aceleracion_geopotencial(r, lambda, sinphi, cosphi, s, rhoP, C, S, Ksub, mu, phi);
            a_u = a_u(:);

            [gx,gy,gz] = gravitysphericalharmonic(r_ecef.', "Custom", deg, {tmpfile, @load}, "None");
            a_tb = [gx;gy;gz];

            errs(k) = norm(a_u - a_tb);
        end

        testName = sprintf("G3 term-by-term: n=%d m=%d", n, m);
        results(end+1) = packResult(testName, errs, 1e-10, "");

        termName(ii) = string(testName);
        termMax(ii)  = max(errs);
    end

    % Ranking rápido: cuál término es el culpable
    [termMaxSorted, idx] = sort(termMax, 'descend');
    fprintf('\n--- G3 ranking (peor termino primero) ---\n');
    for jj = 1:numel(idx)
        fprintf('%-35s  maxErr=% .3e\n', termName(idx(jj)), termMaxSorted(jj));
    end
    fprintf('--- FIN G3 ranking ---\n');

end

%% ============================
% Reporte
% ============================
fprintf('\n--- REPORTE ---\n');
for i=1:numel(results)
    r = results(i);
    fprintf('%-55s  %-6s', r.test, r.status);
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

function out = packSkip(name, note)
out = struct('test',string(name),'status',"SKIP",'maxErr',nan,'rmsErr',nan,'note',string(note));
end

function msg = missingWhy(hasMy, hasTb, myName, tbName)
msg = "";
if ~hasMy
    msg = msg + "no existe " + myName + ". ";
end
if ~hasTb
    msg = msg + "no existe toolbox " + tbName + ". ";
end
msg = strtrim(msg);
end

function ang = wrapTo180Local(angdeg)
ang = mod(angdeg + 180, 360) - 180;
end

function ang = wrapToPiLocal(angrad)
ang = mod(angrad + pi, 2*pi) - pi;
end

function a = accelJ2_fromC20norm(r_ecef, mu, Re, C20_norm)
% Para coeficientes fully-normalized tipo EGM: C20_norm = -J2/sqrt(5)
J2 = -C20_norm*sqrt(5);

x = r_ecef(1); y = r_ecef(2); z = r_ecef(3);
r2 = x*x + y*y + z*z;
r  = sqrt(r2);

zx = (z*z)/r2;
k  = 1.5*J2*(Re*Re/r2);

fxy = 1 - k*(5*zx - 1);
fz  = 1 - k*(5*zx - 3);

a = [-mu*x/r^3 * fxy;
     -mu*y/r^3 * fxy;
     -mu*z/r^3 * fz];
end
