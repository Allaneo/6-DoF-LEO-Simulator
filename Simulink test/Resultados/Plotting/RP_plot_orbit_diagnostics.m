function RP_plot_orbit_diagnostics(D, cfg, H)
shadowSegments = D.derived.shadowSegments;

tx = D.tx(:);
r  = D.rECI;
v  = D.vECI;

nT = size(r,1);
a = nan(nT,1); e = nan(nT,1); inc = nan(nT,1);
raan = nan(nT,1); argp = nan(nT,1); nu = nan(nT,1);

% self-check (comparación con forma tipo post_orbitas para e)
e_alt = nan(nT,1);
a_alt = nan(nT,1);

for ii = 1:nT
    ri = r(ii,:).'; vi = v(ii,:).';
    [a(ii), e(ii), inc(ii), raan(ii), argp(ii), nu(ii)] = H.coe_from_rv(ri, vi, cfg.muEarth);

    % alt:
    h = cross(ri,vi);
    R = norm(ri);
    evec2 = (cross(vi,h)/cfg.muEarth) - (ri/R);
    e_alt(ii) = norm(evec2);

    p = (norm(h)^2)/cfg.muEarth;
    if (1 - e(ii)^2) > 0
        a_alt(ii) = p/(1 - e(ii)^2);
    end
end

% Reporte diferencias (si están bien, deberían ser ~numéricas)
de = max(abs(e - e_alt));
da = max(abs(a - a_alt));
fprintf('\n[Orbit self-check] max |e - e_alt| = %.3e\n', de);
fprintf('[Orbit self-check] max |a - a_alt| = %.3e m\n', da);

% Mean elements opcionales (movmean aprox 1 órbita)
if cfg.PLOT_MEAN_ELEMENTS
    dt = median(diff(tx));
    a_med = median(a(isfinite(a)));
    Torb = 2*pi*sqrt((a_med^3)/cfg.muEarth);
    W = max(3, round(Torb/dt));

    a_m   = local_movmean(a, W);
    e_m   = local_movmean(e, W);
    inc_m = local_movmean(inc, W);

    raan_u = unwrap(raan);
    raan_m = local_movmean(raan_u, W);
else
    a_m = []; e_m = []; inc_m = []; raan_m = [];
end

figure('Name','Órbita - Semieje mayor','Color','w');
ax = axes; hold(ax,'on'); grid(ax,'on');
plot(ax, tx, a/1e3, 'LineWidth', 1.2);
if ~isempty(a_m)
    plot(ax, tx, a_m/1e3, 'LineWidth', 2.0);
    legend(ax, {'Osculante','Media ~1 órbita'}, 'Location','best');
end
xlabel(ax,'Tiempo [s]'); ylabel(ax,'a [km]');
title(ax,'Semieje mayor (osculante)');
H.apply_shadow(ax, shadowSegments);

figure('Name','Órbita - Excentricidad','Color','w');
ax = axes; hold(ax,'on'); grid(ax,'on');
plot(ax, tx, e, 'LineWidth', 1.2);
if ~isempty(e_m)
    plot(ax, tx, e_m, 'LineWidth', 2.0);
    legend(ax, {'Osculante','Media ~1 órbita'}, 'Location','best');
end
xlabel(ax,'Tiempo [s]'); ylabel(ax,'e [-]');
title(ax,'Excentricidad (osculante)');
H.apply_shadow(ax, shadowSegments);

figure('Name','Órbita - Inclinación','Color','w');
ax = axes; hold(ax,'on'); grid(ax,'on');
plot(ax, tx, inc*180/pi, 'LineWidth', 1.2);
if ~isempty(inc_m)
    plot(ax, tx, inc_m*180/pi, 'LineWidth', 2.0);
    legend(ax, {'Osculante','Media ~1 órbita'}, 'Location','best');
end
xlabel(ax,'Tiempo [s]'); ylabel(ax,'i [deg]');
title(ax,'Inclinación (osculante)');
H.apply_shadow(ax, shadowSegments);

figure('Name','Órbita - RAAN','Color','w');
ax = axes; hold(ax,'on'); grid(ax,'on');
plot(ax, tx, H.unwrap_deg(raan*180/pi), 'LineWidth', 1.2);
if ~isempty(raan_m)
    plot(ax, tx, H.unwrap_deg(raan_m*180/pi), 'LineWidth', 2.0);
    legend(ax, {'Osculante','Media ~1 órbita'}, 'Location','best');
end
xlabel(ax,'Tiempo [s]'); ylabel(ax,'\Omega [deg]');
title(ax,'RAAN (osculante)');
H.apply_shadow(ax, shadowSegments);
end

function y = local_movmean(x, W)
% movmean compatible: si no existe movmean en tu MATLAB, hace conv simple.
if exist('movmean','file') == 2
    y = movmean(x, W, 'omitnan');
else
    w = ones(W,1)/W;
    y = conv(x, w, 'same');
end
end
