function RP_plot_groundtrack(D, cfg, H)
% Ground track:
% - si hay phi/lambda logueados => usa eso
% - si no => ECI -> (ERA) -> ECEF -> lat/lon (WGS84)

useLogged = isfield(D.data,'phi') && isfield(D.data,'lambda') && D.data.phi.ok && D.data.lambda.ok;

if useLogged
    tpos = D.data.phi.t(:);
    phi  = D.data.phi.X(:);
    lam  = interp1(D.data.lambda.t(:), D.data.lambda.X(:), tpos, 'linear', 'extrap');

    if strcmpi(cfg.PHI_LAMBDA_UNITS,'rad')
        lat_deg = phi*180/pi;
        lon_deg = lam*180/pi;
    else
        lat_deg = phi;
        lon_deg = lam;
    end
else
    tpos = D.tx(:);

    % ERA
    if isfield(D.data,'era') && D.data.era.ok
        theta = interp1(D.data.era.t(:), D.data.era.X(:), tpos, 'linear', 'extrap');
    else
        if ~isempty(cfg.jd_ut_inicial)
            if cfg.has_MD_era
                theta = zeros(size(tpos));
                for i = 1:numel(tpos)
                    theta(i) = MD_era(tpos(i), cfg.jd_ut_inicial);
                end
            else
                theta = H.era_from_jd_ut1(cfg.jd_ut_inicial, tpos);
            end
        else
            theta = cfg.omegaE * tpos;
        end
    end

    ct = cos(theta); st = sin(theta);
    xeci = D.rECI(:,1); yeci = D.rECI(:,2); zeci = D.rECI(:,3);

    xecef = ct.*xeci + st.*yeci;
    yecef = -st.*xeci + ct.*yeci;
    zecef = zeci;

    % C_FIX si aplica
    if ~isempty(cfg.XFORM_MODE_val) && (cfg.XFORM_MODE_val == 2) && ~isempty(cfg.C_FIX)
        try
            rtmp = [xecef yecef zecef];
            rfix = rtmp * (cfg.C_FIX.');
            xecef = rfix(:,1); yecef = rfix(:,2); zecef = rfix(:,3);
        catch
        end
    end

    [lat_deg, lon_deg] = H.ecef2lla_wgs84_deg(xecef, yecef, zecef);
end

lon_deg = mod(lon_deg + 180, 360) - 180;
jump = [false; abs(diff(lon_deg)) > 180];
lat_deg(jump) = NaN;
lon_deg(jump) = NaN;

% separar iluminado/sombra
shadowOn = false(size(tpos));
if isfield(D.data,'shadow') && D.data.shadow.ok
    shadowOn = interp1(D.data.shadow.t(:), double(D.data.shadow.X(:)>0.5), tpos, 'nearest', 'extrap') > 0.5;
end

lat_sun = lat_deg; lon_sun = lon_deg;
lat_sh  = lat_deg; lon_sh  = lon_deg;
lat_sun(shadowOn) = NaN; lon_sun(shadowOn) = NaN;
lat_sh(~shadowOn) = NaN; lon_sh(~shadowOn) = NaN;

figure('Name','Ground Track - Sombra vs Iluminado','Color','w');
ax = axes; hold(ax,'on'); grid(ax,'on');

didMap = false;
try
    load coastlines
    plot(ax, coastlon, coastlat, 'LineWidth', 0.8, 'HandleVisibility','off');
    didMap = true;
catch
    try
        load coast
        plot(ax, long, lat, 'LineWidth', 0.8, 'HandleVisibility','off');
        didMap = true;
    catch
        didMap = false;
    end
end

hSun = plot(ax, lon_sun, lat_sun, 'LineWidth', 1.6);
hSh  = plot(ax, lon_sh,  lat_sh,  'LineWidth', 1.6, 'LineStyle','--');

xlabel(ax,'Longitud [deg]');
ylabel(ax,'Latitud [deg]');
if didMap
    title(ax, sprintf('Ground track - pos de %s(:,1:3) en ECI', D.xNameUsed));
else
    title(ax, sprintf('Ground track (sin mapa base) - pos de %s(:,1:3) en ECI', D.xNameUsed));
end
xlim(ax,[-180 180]); ylim(ax,[-90 90]);
legend(ax, [hSun hSh], {'Iluminado','Sombra'}, 'Location','best');
end
