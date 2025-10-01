% test_geopotential_sch.m  — test del geopotentialAccel_ECEF + J2

clear; clc;

% 1) Parámetros de carga
egmfile  = 'EGM2008_to2190_TideFree';  % nombre de tu fichero de coefes
Nmax     = 3;
useCache = true;                       % caché en tempdir/EGMcache

[Cbar,Sbar,~,meta] = loadEGMcoeffs(egmfile, Nmax, useCache);

% 2) Parámetros físicos
params.GM = 3.986004418e14;
params.a  = 6378136.3;
params.N  = Nmax;

% 3) Posición en el ecuador, h = 700 km
Re     = params.a;
alt    = 700e3;
r_ecef = [Re+alt; 0; 0];
rnorm  = norm(r_ecef);

% 4) Aceleración geopotencial completa
a_geo = geopotentialAccel_ECEF(r_ecef, Cbar, Sbar, params);
fprintf('a_geo (ECEF)          = [%.6e  %.6e  %.6e] m/s^2\n', a_geo);

% 5) Two-body (central)
a_tb = - params.GM / rnorm^3 * r_ecef;
fprintf('two-body accel        = [%.6e  %.6e  %.6e] m/s^2\n', a_tb);
fprintf('diff total - two-body = %.3e m/s^2\n\n', norm(a_geo - a_tb));

% 6) Perturbación J2 sintética (solo grado 2)
Ctest      = zeros(size(Cbar));
Ctest(3,1) = Cbar(3,1);     % fully-normalized C20
Stest      = zeros(size(Sbar));
a_full2    = geopotentialAccel_ECEF(r_ecef, Ctest, Stest, params);
a_J2_synth = a_full2 - a_tb;

% 7) Closed-form J2: usar exactamente -Cbar(3,1) para coeficiente
J2norm = - Cbar(3,1);  % directamente del fully-normalized C20
mu     = params.GM;
x = r_ecef(1);
z = r_ecef(3);

ax_j2 = (3*mu*Re^2*J2norm)/(2*rnorm^5) * x * (5*z^2/rnorm^2 - 1);
ay_j2 = 0;
az_j2 = (3*mu*Re^2*J2norm)/(2*rnorm^5) * z * (5*z^2/rnorm^2 - 3);
a_j2_cf = [ax_j2; ay_j2; az_j2];

fprintf('J2 synth norm = %.6e m/s^2\n', norm(a_J2_synth));
fprintf('J2 cf    norm = %.6e m/s^2\n', norm(a_j2_cf));
fprintf('diff (synth–cf) = %.3e m/s^2\n', norm(a_J2_synth - a_j2_cf));
