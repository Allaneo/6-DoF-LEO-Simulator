%% TEST_torque_gravitacional.m
clear; clc;

INIT_parametros;

%% Caso de prueba
rng(1);
q = randn(1,4); q = q/norm(q);

% Posición 

posicion_inicial = 6900e3*[0.7;0.7;0.1415];

%% Pre-cálculos independientes de delta
R_ECEF = posicion_inicial; % era=0
[r, lambda, sinphi, cosphi, s, rhoP, phi] = MT_ECEF2GEOCENTRIC(R_ECEF);
a0_ECEF = MD_aceleracion_geopotencial(r, lambda, sinphi, cosphi, s, rhoP, Cbar, Sbar, K, GM, phi);

% Torque extendido para un delta (forward difference Jacobian)
compute_tau = @(delta) MD_momento_gravitatorio( ...
    MT_ECEF2BODY_matrix(q, 0, ...
        MD_Ja_forward(R_ECEF, Sbar, Cbar, K, delta, GM, a0_ECEF) ...
    ), I);

%% Barrido de delta
delta_vec = logspace(-6, 5, 131);  % 1e-5 ... 1e5 m
N = numel(delta_vec);

tauE      = zeros(3,N);
tauE_half = zeros(3,N);
dTau_abs  = nan(1,N);             

eps0 = 1e-15;

for i = 1:N
    d  = delta_vec(i);
    d2 = d/2;

    tau      = compute_tau(d);   tau  = tau(:);
    tau_half = compute_tau(d2);  tau_half = tau_half(:);

    tauE(:,i)      = tau; 
    tauE_half(:,i) = tau_half;

    dTau_abs(i) = norm(tau - tau_half); 
    % la diferencia se calcula con d/2 porque de esta forma el ratio entre
    % los delta es constante independientemente de la cantidad de puntos
    % elegida para delta_vec
end

%% Métrica de "no-colapso" direccional (L1-share)
absTau = abs(tauE); %valor absoluto de cada componente de la matriz que contiene cada torque para cada delta
sumAbs = sum(absTau,1) + eps0;
share  = absTau ./ sumAbs;

minor_share = 1 - max(share,[],1);  % 0 colapso; máximo teórico 2/3

%% Selección automática de delta (mínimo de salto absoluto en la zona interior)
valid = ~isnan(dTau_abs) & isfinite(dTau_abs);
[~, idx_min] = min(dTau_abs(valid));
idxs = find(valid);
idx_min = idxs(idx_min);

% Elegí un rango cercano al mínimo para evitar elegir un punto aislado
good = valid & (dTau_abs <= 2*dTau_abs(idx_min));
idx_good = find(good);

if isempty(idx_good)
    idx_pick = idx_min;
else
    idx_pick = idx_good(round(numel(idx_good)/2));
end

delta_pick = delta_vec(idx_pick);

fprintf('Delta recomendado (según salto absoluto): %.6g m\n', delta_pick);
fprintf('DeltaTau(δ)=||τ(δ)-τ(δ/2)|| = %.3e N·m\n', dTau_abs(idx_pick));
fprintf('minor_share(δ)= %.3f (0 colapso, 2/3 reparto uniforme)\n', minor_share(idx_pick));

absTau = abs(tauE);
sumAbs = sum(absTau,1) + eps0;

tau_max = max(absTau,[],1);          % N·m
s_max   = tau_max ./ sumAbs;         % adimensional, entre 1/3 y 1


%% =========================
% FIG 1: |tau_i| vs delta (log-log)
% =========================
figure('Color','w','Name','Fig1 |tau_i| vs delta');
grid on; hold on;

loglog(delta_vec, abs(tauE(1,:)) + eps0, '-', 'LineWidth', 1.6);
loglog(delta_vec, abs(tauE(2,:)) + eps0, '-', 'LineWidth', 1.6);
loglog(delta_vec, abs(tauE(3,:)) + eps0, '-', 'LineWidth', 1.6);

xline(delta_pick, '--', sprintf('\\delta^* = %.3g m', delta_pick), ...
    'LabelOrientation','horizontal', 'LabelVerticalAlignment','middle');

xlabel('\delta (m)');
ylabel('|\tau_i| (N\cdotm)');
title('Torque extendido: magnitud por componente');
legend('|\tau_x|','|\tau_y|','|\tau_z|','Location','best');

%% =========================
% FIG 2: DeltaTau absoluto vs delta (log-log) — sin normalización
% =========================
figure('Color','w','Name','Fig2 Step refinement jump (points+line)');
grid on; hold on;

loglog(delta_vec, dTau_abs + eps0, '-o', 'LineWidth', 1.2, 'MarkerSize', 4);

set(gca,'XScale','log','YScale','log');
xlim([min(delta_vec) max(delta_vec)]);

xline(delta_pick, '--', sprintf('\\delta^* = %.3g m', delta_pick), ...
    'LabelOrientation','horizontal', 'LabelVerticalAlignment','middle');

xlabel('\delta (m)');
ylabel('||\tau(\delta)-\tau(\delta/2)|| (N\cdotm)');
title('Sensibilidad al paso: salto absoluto al refinar \delta \rightarrow \delta/2');

%% =========================
% FIG 3: minor_share vs delta (arranca en 0)
% =========================
figure('Color','w','Name','Fig3 Dominancia normalizada');
grid on; hold on;

semilogx(delta_vec, s_max, '-', 'LineWidth', 1.6);

yline(1/3, ':', 'reparto uniforme (1/3)');
yline(1.0,  ':', 'colapso total (1)');
ylim([0.3 1.02]);

xline(delta_pick, '--', sprintf('\\delta^* = %.3g m', delta_pick), ...
    'LabelOrientation','horizontal', 'LabelVerticalAlignment','middle');

xlabel('\delta (m)');
ylabel('s_{max} = max(|\tau_i|)/\sum|\tau_i| (-)');
title('Distribución de componentes: dominancia del eje principal');

%% =========================
% FIG 4: Trade-off limpio (Δτ vs minor_share)
% =========================
figure('Color','w','Name','Fig4 Trade-off (Δτ vs s_max)');
grid on; hold on;
set(gca,'XScale','log');

yyaxis left
loglog(delta_vec, dTau_abs + eps0, '-', 'LineWidth', 1.6);
ylabel('||\tau(\delta)-\tau(\delta/2)|| (N\cdotm)');

yyaxis right
semilogx(delta_vec, s_max, '-', 'LineWidth', 1.6);
ylabel('s_{max} (-)');
ylim([0.3 1.02]);

xline(delta_pick, '--', sprintf('\\delta^* = %.3g m', delta_pick), ...
    'LabelOrientation','horizontal', 'LabelVerticalAlignment','middle');

xlabel('\delta (m)');
title('Selección de \delta: estabilidad numérica vs dominancia (distribución)');
legend('Salto absoluto','s_{max}', 'Location','best');
