function RP_plot_forces(D, cfg, H)
shadowSegments = D.derived.shadowSegments;

forceKeys   = {'F_drag','F_solar','F_grav_nonrad'};
forceLabels = {'Drag','Solar','Gravedad NO radial'};

figure('Name','Perturbaciones - Fuerzas (Magnitud)','Color','w');
ax = axes; hold(ax,'on'); grid(ax,'on');

hF = []; lblF = {};
allVals = [];

for i = 1:numel(forceKeys)
    k = forceKeys{i};
    if ~isfield(D.data,k) || ~D.data.(k).ok, continue; end

    t = D.data.(k).t(:);
    X = D.data.(k).X;

    mag_raw  = H.mag_timeseries_samplewise(X,t);
    mag_plot = mag_raw;

    % Solar: cortar en sombra si shadow existe
    if strcmp(k,'F_solar') && isfield(D.data,'shadow') && D.data.shadow.ok
        sh = interp1(D.data.shadow.t(:), double(D.data.shadow.X(:)>0.5), t, 'nearest', 'extrap') > 0.5;
        mag_plot(sh) = NaN;
    end

    % Log plot: ceros / negativos -> NaN
    mag_raw(mag_raw<=0)   = NaN;
    mag_plot(mag_plot<=0) = NaN;

    h = semilogy(ax, t, mag_plot, 'LineWidth', 1.4);

    % Plot punteado del raw (si difiere por sombra)
    if strcmp(k,'F_solar')
        try
            col = get(h,'Color');
            semilogy(ax, t, mag_raw, ':', 'LineWidth', 1.0, 'Color', col, 'HandleVisibility','off');
        catch
            semilogy(ax, t, mag_raw, ':', 'LineWidth', 1.0, 'HandleVisibility','off');
        end
    end

    hF(end+1) = h; %#ok<AGROW>
    lblF{end+1} = forceLabels{i}; %#ok<AGROW>

    allVals = [allVals; mag_plot(:)]; %#ok<AGROW>
end

if isempty(hF)
    text(0.02,0.5,'No hay señales suficientes para fuerzas.','Units','normalized');
    return;
end

xlabel(ax,'Tiempo [s]');
ylabel(ax,'|F|');
title(ax,'Fuerzas perturbadoras (sin fuerza total)');

H.apply_shadow(ax, shadowSegments);
legend(ax, hF, lblF, 'Location','best');

% Evitar que la escala log "se vaya" a valores absurdos por datos muy chicos
vals = allVals(isfinite(allVals) & allVals>0);
if ~isempty(vals)
    vmax = max(vals);
    vmin = min(vals);
    ylo = max([vmin, vmax*1e-12, eps]);  % piso relativo y eps, pero NO realmin
    set(ax,'YLim',[ylo, vmax*1.05]);
end
end
