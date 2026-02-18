%% RESULTS_plotter_main.m (SCRIPT)
% Plotter modular.
% Convención fija: q=[q0 q1 q2 q3] scalar-first define C_IB (Body->Inertial).

clc;
nl = sprintf('\n');

%% =========================
% Config
% =========================
cfg = struct();
cfg.omegaE = 7.2921150e-5;          % [rad/s]
cfg.muEarth = 3.986004418e14;       % [m^3/s^2]
cfg.body_axis_for_nadir = [0 0 1];  % eje body a comparar contra nadir
cfg.PHI_LAMBDA_UNITS = 'rad';       % 'rad' o 'deg'
cfg.PLOT_ORBIT_DIAGNOSTICS = true;

% Si querés "mean elements" además de osculantes:
cfg.PLOT_MEAN_ELEMENTS = false;     % true => agrega curvas suavizadas (ventana ~1 órbita)

%% =========================
% Señales (nombres exactos)
% =========================
SIG = struct();
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
% Variables opcionales del workspace (para ground track)
% =========================
cfg.jd_ut_inicial = [];
if exist('jd_ut_inicial','var') == 1
    cfg.jd_ut_inicial = jd_ut_inicial;
end

cfg.XFORM_MODE_val = [];
if exist('XFORM_MODE','var') == 1
    cfg.XFORM_MODE_val = RP_helpers().get_param_value(XFORM_MODE);
end

cfg.C_FIX = [];
if exist('C_FIX','var') == 1
    cfg.C_FIX = RP_helpers().get_param_value(C_FIX);
elseif exist('C_FIX_val','var') == 1
    cfg.C_FIX = C_FIX_val;
end

cfg.has_MD_era = (exist('MD_era','file') == 2) || (exist('MD_era','file') == 6);

%% =========================
% Helpers
% =========================
H = RP_helpers();

%% =========================
% Cargar logs (logsout o SDI)
% =========================
logsout_in = [];
if exist('logsout','var') == 1
    logsout_in = logsout;
end

D = RP_load_data(logsout_in, SIG, sigKeys, sigNames, H);

fprintf('Fuente de datos: %s\n', D.src_mode);

%% =========================
% Derivados (shadowSegments, grav_nonrad, actitud derivada)
% =========================
D = RP_add_derived(D, cfg, H);

%% =========================
% Plots
% =========================
RP_plot_forces(D, cfg, H);
RP_plot_moments(D, cfg, H);
RP_plot_groundtrack(D, cfg, H);
RP_plot_attitude(D, cfg, H);

if cfg.PLOT_ORBIT_DIAGNOSTICS
    RP_plot_orbit_diagnostics(D, cfg, H);
end

RP_print_summary(D);

fprintf('\nListo.\n');
