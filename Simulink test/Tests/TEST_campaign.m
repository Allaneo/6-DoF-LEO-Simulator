function test_campaign(opts)
% TEST_CAMPAIGN  Validación simple/teórica desde logsout.
% Tests: SANITY_STATE, J2_RAAN, DRAG, DRAG_DA, RV_SIGN, GG_TORQUE, MAG_STEP
%
% Requiere: logsout con "x_lineal" = [r_eci; v_eci]
% Usa si existen: rho, v_rel, Fuerza Drag, x_angular, Momento Gravedad, Momento Geomagnético TOTAL

clc;
if nargin < 1, opts = struct(); end

% ---------------------
% Constantes
% ---------------------
mu  = 3.986004418e14;      % [m^3/s^2]
Re  = 6378137.0;           % [m]
J2  = 1.08262668e-3;       % [-]

% ---------------------
% Umbrales
% ---------------------
TOL.nan_frac_max        = getOpt(opts,'nan_frac_max', 0.0);
TOL.rraan_rel_pct_max   = getOpt(opts,'rraan_rel_pct_max', 5.0);     % [%]
TOL.rraan_abs_deg_day   = getOpt(opts,'rraan_abs_deg_day', 1e-6);    % [deg/day]
TOL.drag_CdA_CV_max     = getOpt(opts,'drag_CdA_CV_max', 0.001);
TOL.drag_da_rel_pct_max = getOpt(opts,'drag_da_rel_pct_max', 10);    % [%]
TOL.rv_sign_max_frac    = getOpt(opts,'rv_sign_max_frac', 0.001);
TOL.gg_rel_rms_max      = getOpt(opts,'gg_rel_rms_max', 0.02);

fit_window = getOpt(opts,'fit_window',[0 Inf]);
GG_MAXPTS  = getOpt(opts,'gg_maxpts',50000);

ZERO_TOL_ABS = getOpt(opts,'zero_tol_abs',1e-18);
STEP_TOL_ABS = getOpt(opts,'step_tol_abs',1e-18);

% ---------------------
% Señales (nombres exactos)
% ---------------------
SIG.x_lineal  = 'x_lineal';
SIG.x_angular = 'x_angular';
SIG.rho       = 'rho';
SIG.v_rel     = 'v_rel';
SIG.F_drag    = 'Fuerza Drag';
SIG.M_grav    = 'Momento Gravedad';
SIG.M_geo     = 'Momento Geomagnético TOTAL';

% =========================
% 0) logsout + selectores init (si existen)
% =========================
base = 'base';
if evalin(base,"exist('logsout','var')") ~= 1
    error('No encontre logsout en el workspace.');
end
logsout = evalin(base,'logsout');

GRAV_MODEL  = try_eval(base,'GRAV_MODEL', NaN);
ATM_SEL     = try_eval_param(base,'ATM_SEL', NaN);
GEOMAG_SEL  = try_eval_param(base,'GEOMAG_SEL', NaN);
SRP_SEL     = try_eval_param(base,'SRP_SEL', NaN);
SAT_SHAPE   = try_eval_param(base,'SAT_SHAPE', NaN);

fprintf('\n=== test_campaign ===\n');
fprintf('Config: GRAV_MODEL=%s | ATM_SEL=%s | GEOMAG_SEL=%s | SRP_SEL=%s | SAT_SHAPE=%s\n', ...
    grav_label(GRAV_MODEL), sel_label(ATM_SEL, {'US','EXP','NONE'}), sel_label(GEOMAG_SEL, {'ON','OFF'}), ...
    sel_label(SRP_SEL, {'ON','OFF'}), sel_label(SAT_SHAPE, {'Cannonball','Prisma'}));

% =========================
% 1) Leer estado obligatorio
% =========================
[t, Xlin] = getSig(logsout, SIG.x_lineal, []);
t = t(:);

if size(Xlin,2) < 6
    error('"%s" debe tener >=6 estados (r(1:3), v(4:6)).', SIG.x_lineal);
end
r_eci = Xlin(:,1:3);
v_eci = Xlin(:,4:6);

oe = oe_from_rv_vec(r_eci, v_eci, mu);

% =========================
% 2) Tabla resultados
% =========================
rows = repmat(blank_row(),0,1);

add = @(test,metric,sim,exp,units,crit,status,note) add_row(rows,test,metric,sim,exp,units,crit,status,note);

% =====================================================================
% TEST: SANITY_STATE
% =====================================================================
print_header('SANITY_STATE','No NaNs en estado y dt razonable.');
nan_frac = mean(~all(isfinite(Xlin),2));
dt = diff(t);

rows = add_row(rows,'SANITY_STATE','nan_frac', nan_frac, NaN,'(-)', ...
    sprintf('Se espera nan_frac <= %.3g',TOL.nan_frac_max), pass(nan_frac<=TOL.nan_frac_max), '');

rows = add_row(rows,'SANITY_STATE','dt_min', min(dt), NaN,'(s)', ...
    'Se espera dt_min > 0 (tiempo estrictamente creciente)', pass(min(dt)>0), '');

rows = add_row(rows,'SANITY_STATE','dt_max', max(dt), NaN,'(s)', ...
    'INFO (depende solver). Si es enorme, puede ocultar dinamica.', 'INFO', '');

% =====================================================================
% TEST: J2_RAAN
% =====================================================================
print_header('J2_RAAN','Precesion RAAN: sim vs J2 (o ~0 si gravedad central).');
y = unwrap(oe.RAAN);
[slope_sim, R2] = fit_slope(t, y, fit_window);
sim_deg_day = slope_sim*(180/pi)*86400;

a_ref = median_finite(oe.a);
e_ref = median_finite(oe.e);
i_ref = median_finite(oe.i);

if ~isfinite(sim_deg_day) || ~isfinite(a_ref) || ~isfinite(e_ref) || ~isfinite(i_ref)
    rows = add_row(rows,'J2_RAAN','Omega_dot',sim_deg_day,NaN,'(deg/day)', ...
        'SKIP: valores no finitos', 'SKIP','');
else
    if GRAV_MODEL == 3
        exp_deg_day = 0;
        err_abs = abs(sim_deg_day);
        crit = sprintf('Esperado ~0. OK si |sim| <= %.3g deg/day',TOL.rraan_abs_deg_day);
        status = pass(err_abs <= TOL.rraan_abs_deg_day);
    else
        exp_rad_s = raan_rate_J2(a_ref,e_ref,i_ref,mu,Re,J2);
        exp_deg_day = exp_rad_s*(180/pi)*86400;

        if abs(exp_deg_day) < 1e-12
            err_abs = abs(sim_deg_day - exp_deg_day);
            crit = sprintf('Teoria ~0. OK si |sim-exp| <= %.3g deg/day',TOL.rraan_abs_deg_day);
            status = pass(err_abs <= TOL.rraan_abs_deg_day);
        else
            err_pct = 100*(sim_deg_day-exp_deg_day)/exp_deg_day;
            crit = sprintf('OK si |err| <= %.2f %% (ref J2)',TOL.rraan_rel_pct_max);
            status = pass(abs(err_pct) <= TOL.rraan_rel_pct_max);
        end
    end

    note = sprintf('fit_R2=%.4f | a=%.3f km e=%.5f i=%.2f deg',R2,a_ref/1e3,e_ref,i_ref*180/pi);
    rows = add_row(rows,'J2_RAAN','Omega_dot',sim_deg_day,exp_deg_day,'(deg/day)',crit,status,note);
end

rows = add_row(rows,'J2_RAAN','fit_R2',R2,NaN,'(-)', ...
    'INFO: calidad del ajuste lineal', 'INFO','');
% =====================================================================
% TEST: DRAG (CdA)
% =====================================================================
print_header('DRAG','Reconstruye CdA(t)=2|F|/(rho|v|^2) y chequea constancia.');

doDrag = hasSig(logsout,SIG.rho) && hasSig(logsout,SIG.v_rel) && hasSig(logsout,SIG.F_drag);
CdA_med = NaN;

if ~doDrag
    rows = add_row(rows,'DRAG','CdA_median',NaN,NaN,'(m^2)', ...
        'SKIP: faltan rho/v_rel/Fuerza Drag', 'SKIP','');
else
    [~, rho]   = getSig(logsout,SIG.rho,   t); rho = rho(:);
    [~, vrel]  = getSig(logsout,SIG.v_rel, t);
    [~, Fdrag] = getSig(logsout,SIG.F_drag,t);

    vmag = vecnorm(vrel,2,2);
    Fmag = vecnorm(Fdrag,2,2);

    m = isfinite(rho)&(rho>0)&isfinite(vmag)&(vmag>0)&isfinite(Fmag)&(Fmag>0);

    if nnz(m) < 50 || max(Fmag(m)) < ZERO_TOL_ABS
        rows = add_row(rows,'DRAG','CdA_median',NaN,NaN,'(m^2)', ...
            'SKIP: drag trivial o pocas muestras validas','SKIP','');
        CdA_med = NaN;
    else
        CdA = 2*Fmag(m)./(rho(m).*vmag(m).^2);
        CdA_med  = median(CdA);
        CdA_mean = mean(CdA);
        CdA_CV   = std(CdA)/max(abs(CdA_mean),eps);

        rows = add_row(rows,'DRAG','CdA_median',CdA_med,NaN,'(m^2)', ...
            'INFO: CdA desde logs', 'INFO', sprintf('N=%d',nnz(m)));

        if SAT_SHAPE == 1
            crit = sprintf('Cannonball: CdA ~const. OK si CV<=%.3g',TOL.drag_CdA_CV_max);
            status = pass(CdA_CV <= TOL.drag_CdA_CV_max);
        else
            crit = 'INFO: si area efectiva varia, CV no aplica como FAIL.';
            status = 'INFO';
        end

        rows = add_row(rows,'DRAG','CdA_CV',CdA_CV,NaN,'(-)',crit,status,'');
    end
end
% =====================================================================
% TEST: DRAG_DA
% =====================================================================
print_header('DRAG_DA','da/dt: sim (orbit-mean) vs (i) energia del drag (consistencia) y (ii) teoria burda.');

if ~doDrag || ~isfinite(CdA_med) || CdA_med<=0
    rows = add_row(rows,'DRAG_DA','da_dt',NaN,NaN,'(m/s)', ...
        'SKIP: no hay drag o CdA invalido','SKIP','');
else
    m_sat = try_eval(base,'m',NaN);
    if ~isfinite(m_sat) || m_sat<=0
        rows = add_row(rows,'DRAG_DA','da_dt',NaN,NaN,'(m/s)', ...
            'SKIP: falta m en workspace','SKIP','');
    else
        B = m_sat/CdA_med;

        % --- ventana efectiva (no Inf) ---
        t1 = max(fit_window(1), t(1));
        t2 = min(fit_window(2), t(end));
        tspan = max(t2 - t1, 0);

        % --- periodo orbital aproximado ---
        a0 = median_finite(oe.a);
        Torb = 2*pi*sqrt(a0^3/mu);
        Norb = tspan / Torb;

        % --- sim: orbit-mean de a(t) + slope robusto (median de slopes) ---
        [tc, a_orb] = orbit_mean(t, oe.a, Torb, t1, t2);

        if isempty(tc) || numel(tc) < 10
            rows = add_row(rows,'DRAG_DA','da_dt',NaN,NaN,'(m/s)', ...
                'SKIP: no hay suficientes orbit-means','SKIP', ...
                sprintf('Torb=%.1f s | Norb=%.2f',Torb,Norb));
        else
            slopes = diff(a_orb) ./ diff(tc);
            slopes = slopes(isfinite(slopes));
            da_dt_sim = median(slopes);  % robusto (no depende de linealidad global)

            % --- (i) "expected" por consistencia energética del drag ---
            % E = -mu/(2a)  =>  da/dt = (2a^2/mu) * (v·a_drag)
            [~, Fdrag] = getSig(logsout, SIG.F_drag, t);   % asumido en el MISMO frame que v_eci (ideal: ECI)
            a_drag = Fdrag / m_sat;
            pow = sum(v_eci .* a_drag, 2);                % [m^2/s^3] = dE/dt por unidad de masa
            da_dt_inst = (2*(oe.a.^2)/mu) .* pow;          % [m/s]
            [~, da_orb_energy] = orbit_mean(t, da_dt_inst, Torb, t1, t2);

            da_dt_energy = mean(da_orb_energy(isfinite(da_orb_energy)));

            % --- (ii) teoria burda (solo para referencia) ---
            [~, rho_all] = getSig(logsout, SIG.rho, t); rho_all=rho_all(:);
            mw = (t>=t1)&(t<=t2)&isfinite(rho_all)&rho_all>0&isfinite(oe.a);
            rho_mean = mean(rho_all(mw));
            a_mean   = mean(oe.a(mw));
            da_dt_th = -(rho_mean/B)*sqrt(mu*a_mean);

            % --- métrica de consistencia (sim vs energía) ---
            err_pct_energy = 100*(da_dt_sim - da_dt_energy)/max(abs(da_dt_energy), eps);

            % Criterio: si energía no cierra, hay bug (frames/signo/factor)
            crit = 'OK si |sim - energia| <= 10% (consistencia energetica del drag).';
            status = pass(abs(err_pct_energy) <= 10);

            note = sprintf(['sim=median(diff(a_orb)/diff(t)) | energy=mean(orbit-mean(da_dt_inst)) | ' ...
                            'err_energy=%.1f%% | th=%.3e | rho_mean=%.3e | B=%.2f | a0=%.3f km | Torb=%.1f s | Norb=%.2f'], ...
                            err_pct_energy, da_dt_th, rho_mean, B, a0/1e3, Torb, Norb);

            % Reporto: sim vs ENERGY como "expected" (es el chequeo real)
            rows = add_row(rows,'DRAG_DA','da_dt', da_dt_sim, da_dt_energy,'(m/s)', crit, status, note);

            % Y además dejo la teoría burda como INFO separado (para que la veas, sin mezclar)
            rows = add_row(rows,'DRAG_DA','da_dt_th', da_dt_th, NaN,'(m/s)', ...
                'INFO: teoria burda -(rho/B)*sqrt(mu*a) (solo referencia).', 'INFO','');
        end
    end
end


% =====================================================================
% TEST: RV_SIGN
% =====================================================================
print_header('RV_SIGN','Heurístico: fracción de instantes dr>0 y dv>0.');

rnorm = vecnorm(r_eci,2,2);
vnorm = vecnorm(v_eci,2,2);
dr = diff(rnorm); dv = diff(vnorm);

m = isfinite(dr)&isfinite(dv);
frac = mean((dr(m)>0) & (dv(m)>0));

crit = sprintf('OK si frac <= %.3g (heuristico)',TOL.rv_sign_max_frac);
rows = add_row(rows,'RV_SIGN','frac_drpos_dvpos',frac,NaN,'(-)',crit,pass(frac<=TOL.rv_sign_max_frac), ...
    'Si es alto: sospecha integracion/frames o fuerza mal aplicada.');

% =====================================================================
% TEST: GG_TORQUE
% =====================================================================
print_header('GG_TORQUE','tau = 3*mu/r^3 * (rhat_B x (I*rhat_B))');

doGG = hasSig(logsout,SIG.x_angular) && hasSig(logsout,SIG.M_grav);
if ~doGG
    rows = add_row(rows,'GG_TORQUE','rel_rms',NaN,NaN,'(-)', ...
        'SKIP: faltan x_angular y/o Momento Gravedad','SKIP','');
else
    [~, xang]  = getSig(logsout,SIG.x_angular,t);
    [~, Mgrav] = getSig(logsout,SIG.M_grav,t);

    q = xang(:,1:4);          % [N x 4]
    M_sim = Mgrav;            % [N x 3]

    Ibody = try_eval(base,'I',[]);
    if isempty(Ibody)
        rows = add_row(rows,'GG_TORQUE','rel_rms',NaN,NaN,'(-)', 'SKIP: falta I en workspace','SKIP','');
    else
        N = numel(t);
        stride = max(1, ceil(N/GG_MAXPTS));
        idx = 1:stride:N;

        rhatI = r_eci(idx,:);
        rnorm = vecnorm(rhatI,2,2);
        rhatI = rhatI ./ max(rnorm,eps);

        Mk_th = zeros(numel(idx),3);

        for k=1:numel(idx)
            qk = q(idx(k),:).';      % 4x1
            rhI = rhatI(k,:).';      % 3x1

            % ACA se usa tu script. Sin wrappers, sin fallback.
            rhB = MT_ECI2BODY(qk, rhI);

            Mk_th(k,:) = (3*mu/(rnorm(k)^3) * cross(rhB, Ibody*rhB)).';
        end

        Ms = M_sim(idx,:);
        d = vecnorm(Ms - Mk_th,2,2);
        a = vecnorm(Ms,2,2);
        rel_rms = sqrt(mean(d.^2)) / max(sqrt(mean(a.^2)),eps);

        crit = sprintf('OK si rel_rms <= %.3g',TOL.gg_rel_rms_max);
        rows = add_row(rows,'GG_TORQUE','rel_rms',rel_rms,NaN,'(-)',crit,pass(rel_rms<=TOL.gg_rel_rms_max), ...
            sprintf('stride=%d (N=%d)',stride,N));
    end
end

% =====================================================================
% TEST: MAG_STEP
% =====================================================================
print_header('MAG_STEP','Saltos en torque geomagnético: magnitud absoluta y relativa (entre muestras).');

if ~hasSig(logsout,SIG.M_geo)
    rows = add_row(rows,'MAG_STEP','step_abs_mean',NaN,NaN,'(N*m)', ...
        'SKIP: falta Momento Geomagnético TOTAL','SKIP','');
else
    [~, Mgeo] = getSig(logsout,SIG.M_geo,t);
    Mgeo = Mgeo(:,1:3);

    Mnorm = vecnorm(Mgeo,2,2);

    if max(Mnorm) < ZERO_TOL_ABS
        rows = add_row(rows,'MAG_STEP','step_abs_mean',0,NaN,'(N*m)', ...
            'SKIP: M_geo ~ 0','SKIP','');
    else
        % Saltos absolutos entre muestras
        dM = vecnorm(diff(Mgeo,1,1),2,2);      % [N-1 x 1]
        Mprev = Mnorm(1:end-1);               % |M_k| para normalizar ΔM(k)=M(k+1)-M(k)

        mask_abs = isfinite(dM) & (dM > STEP_TOL_ABS);

        if ~any(mask_abs)
            rows = add_row(rows,'MAG_STEP','step_abs_mean',NaN,NaN,'(N*m)', ...
                'SKIP: sin escalones detectables (umbral)','SKIP', ...
                sprintf('umbral |ΔM|>%.3g', STEP_TOL_ABS));
        else
            dM_use = dM(mask_abs);
            step_abs_mean = mean(dM_use);
            step_abs_max  = max(dM_use);

            % Referencia de escala para interpretar "relativo"
            M_med = median_finite(Mnorm);

            % Evitar explosión del relativo cuando |M_k| es muy chico:
            % se excluyen muestras con |M_k| <= denom_min
            denom_min = max(1e-3*M_med, ZERO_TOL_ABS);

            mask_rel = mask_abs & isfinite(Mprev) & (Mprev > denom_min);
            rel = dM(mask_rel) ./ Mprev(mask_rel);

            % Default si no hay relativos válidos
            step_rel_mean = NaN;
            step_rel_max  = NaN;
            M_at_rel_max  = NaN;
            dM_at_rel_max = NaN;

            if any(mask_rel)
                step_rel_mean = mean(rel);
                [step_rel_max, jmax] = max(rel);

                idx_rel = find(mask_rel);      % índices en 1..N-1
                kmax = idx_rel(jmax);          % salto entre kmax y kmax+1
                M_at_rel_max  = Mprev(kmax);
                dM_at_rel_max = dM(kmax);
            end

            % --- Reporte (simple y con contexto) ---
            rows = add_row(rows,'MAG_STEP','M_med', M_med, NaN,'(N*m)', ...
                'Referencia: median(|M_geo|).', 'INFO','');

            rows = add_row(rows,'MAG_STEP','step_abs_mean', step_abs_mean, NaN,'(N*m)', ...
                'INFO: mean(|ΔM|) entre muestras', 'INFO', ...
                sprintf('N=%d | umbral |ΔM|>%.3g', nnz(mask_abs), STEP_TOL_ABS));

            rows = add_row(rows,'MAG_STEP','step_abs_max', step_abs_max, NaN,'(N*m)', ...
                'INFO: max(|ΔM|) entre muestras', 'INFO','');

            rows = add_row(rows,'MAG_STEP','step_rel_mean', step_rel_mean, NaN,'(-)', ...
                'INFO: mean(|ΔM|/|M_k|) relativo local', 'INFO', ...
                sprintf('excluye |M_k|<=%.3e', denom_min));

            rows = add_row(rows,'MAG_STEP','step_rel_max', step_rel_max, NaN,'(-)', ...
                'INFO: max(|ΔM|/|M_k|) relativo local', 'INFO', ...
                sprintf('|M_k|@max=%.3e N*m | |ΔM|@max=%.3e N*m | excluye |M_k|<=%.3e', ...
                        M_at_rel_max, dM_at_rel_max, denom_min));
        end
    end
end
% =========================
% Reporte final + export
% =========================
T = struct2table(rows);
assignin(base,'pps_validation_table',T);

fprintf('\n=== REPORTE (pps_validation_table) ===\n');
print_report(T);
fprintf('\nTabla guardada en workspace: pps_validation_table\n');
fprintf('Listo.\n');

end

% ==========================
% Helpers 
% ==========================

function print_header(name, desc)
fprintf('\n[%s]\n  %s\n', name, desc);
end

function s = pass(tf)
if tf, s='OK'; else, s='FAIL'; end
end

function R = blank_row()
R = struct('test','','metric','','sim',NaN,'expected',NaN,'err_abs',NaN,'err_pct',NaN, ...
           'units','','status','','criteria','','note','');
end

function rows = add_row(rows,test,metric,sim,exp,units,crit,status,note)
R = blank_row();
R.test = char(test); R.metric = char(metric);
R.sim = sim; R.expected = exp;
R.units = char(units); R.criteria = char(crit);
R.status = char(status); R.note = char(note);

if isfinite(sim) && isfinite(exp)
    R.err_abs = sim - exp;
    if abs(exp)>0, R.err_pct = 100*(sim-exp)/exp; end
end
rows(end+1,1) = R; 
end

function tf = hasSig(logsout, name)
% True si existe un elemento con ese Name en logsout (sin warnings)
tf = ~isempty(getElQuiet(logsout, name));
end

function [t, D] = getSig(logsout, name, t_master)
% Devuelve D como [N x M]. Si t_master no vacío, interpola a t_master.
el = getElQuiet(logsout, name);
if isempty(el)
    error('No encontre "%s" en logsout.', char(name));
end

ts = el.Values;
t  = ts.Time(:);
D  = ts_data_to_mat(ts.Data, numel(t));

if ~isempty(t_master)
    if numel(t)==numel(t_master) && all(abs(t - t_master(:)) < 1e-12)
        t = t_master(:);
        return
    end
    D = interp1(t, D, t_master(:), 'linear', 'extrap');
    t = t_master(:);
end
end

function el = getElQuiet(logsout, name)
% Lookup silencioso (sin warnings) por Name exacto y luego case-insensitive
el = [];
target = strtrim(char(name));

nEl = logsout.numElements;

% match exacto
for i = 1:nEl
    nm = '';
    try
        nm = strtrim(logsout{i}.Name);
    catch
        nm = '';
    end
    if strcmp(nm, target)
        el = logsout{i};
        return
    end
end

% match case-insensitive
for i = 1:nEl
    nm = '';
    try
        nm = strtrim(logsout{i}.Name);
    catch
        nm = '';
    end
    if strcmpi(nm, target)
        el = logsout{i};
        return
    end
end
end

function D = ts_data_to_mat(raw, N)
% Convierte ts.Data a matriz [N x M] independientemente de dónde esté el eje tiempo.
% Ejemplos válidos:
%   raw = [N x M]
%   raw = [M x N]
%   raw = [M x 1 x N]
%   raw = [1 x M x N]
% etc.

raw = squeeze(raw);

% Caso vector
if isvector(raw) && numel(raw) == N
    D = raw(:);
    return
end

sz = size(raw);
idx = find(sz == N, 1, 'last');   % toma la última dimensión que coincide con N
if isempty(idx)
    error('Data incompatible con Time. size(Data)=%s, Ntime=%d', mat2str(sz), N);
end

% Lleva el eje tiempo a la primera dimensión
perm = [idx, setdiff(1:ndims(raw), idx, 'stable')];
tmp = permute(raw, perm);
tmp = squeeze(tmp);

if size(tmp,1) ~= N
    error('No pude alinear Data con Time. size=%s, Ntime=%d', mat2str(size(tmp)), N);
end

D = reshape(tmp, N, []);
end

function oe = oe_from_rv_vec(r, v, mu)
R = vecnorm(r,2,2);
V = vecnorm(v,2,2);

h = cross(r, v, 2);
hn = vecnorm(h,2,2);

n = [-h(:,2), h(:,1), zeros(size(h,1),1)];
rv = sum(r.*v,2);

evec = (1/mu)*((V.^2 - mu./R).*r - rv.*v);
e = vecnorm(evec,2,2);

a = 1 ./ (2./R - (V.^2)/mu);

ci = h(:,3)./max(hn,eps);
ci = max(min(ci,1),-1);
inc = acos(ci);

RAAN = atan2(n(:,2), n(:,1));
RAAN(RAAN<0) = RAAN(RAAN<0) + 2*pi;

oe.a=a; oe.e=e; oe.i=inc; oe.RAAN=RAAN;
end

function [slope, R2] = fit_slope(t, y, win)
t=t(:); y=y(:);
m = (t>=win(1))&(t<=win(2))&isfinite(y);
tt=t(m); yy=y(m);
if numel(tt)<10, slope=NaN; R2=NaN; return; end
p = polyfit(tt,yy,1);
yhat = polyval(p,tt);
SSres = sum((yy-yhat).^2);
SStot = sum((yy-mean(yy)).^2);
R2 = 1 - SSres/max(SStot,eps);
slope = p(1);
end

function w = raan_rate_J2(a,e,i,mu,Re,J2)
p = a*(1-e^2);
n = sqrt(mu/a^3);
w = -1.5*J2*n*(Re/p)^2*cos(i);
end

function med = median_finite(x)
x = x(isfinite(x));
if isempty(x), med=NaN; else, med=median(x); end
end

function print_report(T)
tests = unique(T.test,'stable');
for k=1:numel(tests)
    tk = tests{k};
    Tk = T(strcmp(T.test,tk),:);

    fprintf('\n%s\n',tk);
    for i=1:height(Tk)
        sim = Tk.sim(i); ex = Tk.expected(i);
        if isfinite(ex)
            fprintf('  %-18s : sim=%s  exp=%s  err=%s  (%s)  [%s]\n', ...
                Tk.metric{i}, fmt(sim), fmt(ex), fmt(Tk.err_abs(i)), Tk.units{i}, Tk.status{i});
        else
            fprintf('  %-18s : sim=%s  (%s)  [%s]\n', ...
                Tk.metric{i}, fmt(sim), Tk.units{i}, Tk.status{i});
        end
        if ~isempty(Tk.criteria{i}), fprintf('    criterio: %s\n',Tk.criteria{i}); end
        if ~isempty(Tk.note{i}),     fprintf('    note: %s\n',Tk.note{i}); end
    end
end
end

function s = fmt(x)
if ~isfinite(x)
    s='NaN';
elseif abs(x) >= 1e4 || abs(x) < 1e-3
    s=sprintf('%.4e',x);
else
    s=sprintf('%.6g',x);
end
end

function v = getOpt(opts,name,default)
if isfield(opts,name), v=opts.(name); else, v=default; end
end

function x = try_eval(ws,varname,default)
try, x=evalin(ws,varname); catch, x=default; end
end

function x = try_eval_param(ws,varname,default)
try
    p = evalin(ws,varname);
    if isa(p,'Simulink.Parameter'), x=double(p.Value);
    else, x=double(p);
    end
catch
    x=default;
end
end

function s = grav_label(g)
if ~isfinite(g), s='NaN';
elseif g==1, s='1 (N4)';
elseif g==2, s='2 (J2)';
elseif g==3, s='3 (SIMPLE)';
else, s=sprintf('%g (UNKNOWN)',g);
end
end

function s = sel_label(v,names)
if ~isfinite(v), s='NaN'; return; end
ii = round(v);
if ii>=1 && ii<=numel(names)
    s = sprintf('%d (%s)',ii,names{ii});
else
    s = sprintf('%d (UNKNOWN)',ii);
end
end
function [tc, ymean] = orbit_mean(t, y, Torb, t1, t2)
% Promedia y(t) por órbita (bins de tamaño Torb).
% Devuelve:
%   tc    : tiempo medio de cada bin
%   ymean : promedio de y en ese bin

m = (t>=t1) & (t<=t2) & isfinite(y);
tt = t(m); yy = y(m);

if numel(tt) < 10
    tc = []; ymean = [];
    return
end

k = floor((tt - tt(1))/Torb);
K = unique(k);

tc = zeros(numel(K),1);
ymean = zeros(numel(K),1);

for i=1:numel(K)
    mi = (k==K(i));
    tc(i) = mean(tt(mi));
    ymean(i) = mean(yy(mi));
end
end
