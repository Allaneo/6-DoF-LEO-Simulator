function RP_plot_moments(D, cfg, H)
shadowSegments = D.derived.shadowSegments;

momentKeys   = {'M_drag','M_solar','M_grav','M_geo'};
momentLabels = {'Drag','Solar','Gravedad','Geomagnético'};

figure('Name','Perturbaciones - Momentos (Magnitud)','Color','w');
ax = axes; hold(ax,'on'); grid(ax,'on');

hM = []; lblM = {};
allVals = [];

for i = 1:numel(momentKeys)
    k = momentKeys{i};
    if ~isfield(D.data,k) || ~D.data.(k).ok, continue; end

    t = D.data.(k).t(:);
    X = D.data.(k).X;

    mag_raw  = H.mag_timeseries_samplewise(X,t);
    mag_plot = mag_raw;

    % Solar: cortar en sombra si shadow existe
    if strcmp(k,'M_solar') && isfield(D.data,'shadow') && D.data.shadow.ok
        sh = interp1(D.data.shadow.t(:), double(D.data.shadow.X(:)>0.5), t, 'nearest', 'extrap') > 0.5;
        mag_plot(sh) = NaN;
    end

    mag_raw(mag_raw<=0)   = NaN;
    mag_plot(mag_plot<=0) = NaN;

    h = semilogy(ax, t, mag_plot, 'LineWidth', 1.4);

    if strcmp(k,'M_solar')
        try
            col = get(h,'Color');
            semilogy(ax, t, mag_raw, ':', 'LineWidth', 1.0, 'Color', col, 'HandleVisibility','off');
        catch
            semilogy(ax, t, mag_raw, ':', 'LineWidth', 1.0, 'HandleVisibility','off');
        end
    end

    hM(end+1) = h; %#ok<AGROW>
    lblM{end+1} = momentLabels{i}; %#ok<AGROW>

    allVals = [allVals; mag_plot(:)]; %#ok<AGROW>
end

if isempty(hM)
    text(0.02,0.5,'No hay señales suficientes para momentos.','Units','normalized');
    return;
end

xlabel(ax,'Tiempo [s]');
ylabel(ax,'|M|');
title(ax,'Momentos perturbadores (BODY) - sin momento total');

H.apply_shadow(ax, shadowSegments);
legend(ax, hM, lblM, 'Location','best');

% Piso log razonable (sin realmin)
vals = allVals(isfinite(allVals) & allVals>0);
if ~isempty(vals)
    vmax = max(vals);
    vmin = min(vals);
    ylo = max([vmin, vmax*1e-12, eps]);
    set(ax,'YLim',[ylo, vmax*1.05]);
end

%% Barras max por componente SOLO si la señal es vector 3D (si es escalar, se omite)
momentKeys_bar   = {'M_drag','M_solar','M_grav','M_geo'};   % NO total
momentLabels_bar = {'Drag','Solar','Gravedad','Geomagnético'};

MmaxMat = [];
MmaxNames = {};

for i = 1:numel(momentKeys_bar)
    k = momentKeys_bar{i};
    if ~isfield(D.data,k) || ~D.data.(k).ok, continue; end

    t = D.data.(k).t(:);
    X3 = H.vec3_timeseries_samplewise(D.data.(k).X, t);

    if isempty(X3)
        % Es escalar o <3 componentes -> no hay "max por componente".
        continue;
    end

    mx = max(abs(X3), [], 1);
    MmaxMat = [MmaxMat; mx]; %#ok<AGROW>
    MmaxNames{end+1} = momentLabels_bar{i}; %#ok<AGROW>
end

if ~isempty(MmaxMat)
    figure('Name','Perturbaciones - Momentos (Max componentes)','Color','w');
    axB = axes; grid(axB,'on'); hold(axB,'on');
    bar(axB, MmaxMat');
    set(axB,'XTick',1:3,'XTickLabel',{'|M_x|_{max}','|M_y|_{max}','|M_z|_{max}'});
    ylabel(axB,'Max |M| [N·m]');
    title(axB,'Máximos por componente (solo señales vectoriales)');
    legend(axB, MmaxNames, 'Location','best');
end
end
