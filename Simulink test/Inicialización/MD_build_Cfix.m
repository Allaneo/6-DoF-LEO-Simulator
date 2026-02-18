function [Cfix, info] = MD_build_Cfix(utc0, deltaAT)
%MD_build_Cfix_ECI2ECEF  Construye una corrección fija Cfix tal que:
%   C_used(t) = R3(-ERA(t)) * Cfix
% y en t0 se cumpla: C_used(t0) = C_full(t0) (modelo IAU-2000/2006 AeroTB).
%
% Inputs:
%   utc0    : [Y M D h m s] en UTC
%   deltaAT : TAI-UTC [s] (ej: 37)
%
% Outputs:
%   Cfix : 3x3 matriz fija (aplicada a la derecha de R3(-ERA))
%   info : struct con valores de debug + error en t0

    if nargin < 2
        deltaAT = 37;
    end

    % --- JD UTC en t0
    jd_utc0 = juliandate(utc0(1),utc0(2),utc0(3),utc0(4),utc0(5),utc0(6));

    % --- EOP en t0 (necesitás tener esta función en tu proyecto)
    % devuelve: deltaUT1 [s], xp [rad], yp [rad]
    [deltaUT1, xp_rad, yp_rad] = MD_getEOP_finals2000A(utc0);

    % --- JD UT1 para ERA
    jd_ut10 = jd_utc0 + deltaUT1/86400.0;

    % --- DCM "full" en t0 (Aerospace Toolbox)
    polarmotion = [xp_rad yp_rad];
    C_full0 = dcmeci2ecef('IAU-2000/2006', utc0, deltaAT, deltaUT1, polarmotion); % ECI->ECEF

    % --- Tu ERA en t0
    % asumido: MD_era(t_sim, jd_ut1) retorna ERA en rad
    [ERA0, ~] = MD_era(0, jd_ut10);

    % --- Tu rotación "ERA-only" (ECI->ECEF) = R3(-ERA)
    R0 = R3_minus(ERA0);

    % --- Cfix para usar como: C_used(t) = R(t) * Cfix
    % Impone: R0*Cfix = C_full0  ==>  Cfix = R0' * C_full0
    Cfix = R0.' * C_full0;

    % --- Ortonormalizar por robustez numérica (opcional pero recomendado)
    Cfix = orthoDCM(Cfix);

    % --- Info debug
    info = struct();
    info.deltaUT1  = deltaUT1;
    info.xp_rad    = xp_rad;
    info.yp_rad    = yp_rad;
    info.jd_utc0   = jd_utc0;
    info.jd_ut10   = jd_ut10;
    info.ERA0      = ERA0;

    C_used0 = R0 * Cfix;
    info.fro_err_t0 = norm(C_used0 - C_full0, 'fro'); % debería ser ~1e-12..1e-9 típico

end

% ===== Rotación alrededor de Z con signo "ECI->ECEF" =====
function R = R3_minus(era)
    c = cos(era); s = sin(era);
    R = [ c  s  0;
         -s  c  0;
          0  0  1];
end

% ===== Ortonormalización (SVD) =====
function C = orthoDCM(C)
    [U,~,V] = svd(C);
    C = U*V';
    if det(C) < 0
        U(:,3) = -U(:,3);
        C = U*V';
    end
end

function [dut1_sec, xp_rad, yp_rad] = MD_getEOP_finals2000A(utc_vec)
% Lee UT1-UTC y polar motion (xp, yp) para la fecha utc_vec desde:
% https://maia.usno.navy.mil/ser7/finals2000A.data
% Formato de columnas según MathWorks (aeroReadIERSData). :contentReference[oaicite:3]{index=3}

    url = 'https://maia.usno.navy.mil/ser7/finals2000A.data';
    local = fullfile(tempdir, 'finals2000A.data');

    if ~isfile(local)
        try
            websave(local, url);
        catch
            error('No pude descargar finals2000A.data. Bajalo manualmente de USNO y ponelo en: %s', local);
        end
    end

    y = utc_vec(1); m = utc_vec(2); d = utc_vec(3);

    fid = fopen(local,'r');
    if fid < 0, error('No pude abrir %s', local); end

    dut1_sec = 0; xp_rad = 0; yp_rad = 0;
    found = false;

    % finals2000A.data es fixed-width. Usamos las columnas oficiales:
    % Year(1-2), Month(3-4), Day(5-6)
    % PM-x (19-27) arcsec, PM-y (38-46) arcsec, UT1-UTC (59-68) sec :contentReference[oaicite:4]{index=4}
    while ~feof(fid)
        line = fgetl(fid);
        if ~ischar(line) || numel(line) < 68, continue; end

        yy = str2double(line(1:2));
        mm = str2double(line(3:4));
        dd = str2double(line(5:6));

        if isnan(yy) || isnan(mm) || isnan(dd), continue; end

        % Regla del archivo: año real depende de MJD; para fechas modernas es 2000+yy.
        % Para 2026 esto es correcto.
        yyyy = 2000 + yy;

        if yyyy==y && mm==m && dd==d
            pmx_asec = str2double(line(19:27));
            pmy_asec = str2double(line(38:46));
            dut1_sec = str2double(line(59:68)); % UT1-UTC [s]

            asec2rad = (pi/180.0)/3600.0;
            xp_rad = pmx_asec * asec2rad;
            yp_rad = pmy_asec * asec2rad;

            found = true;
            break;
        end
    end
    fclose(fid);

    if ~found
        error('No encontré EOP para %04d-%02d-%02d en finals2000A.data', y,m,d);
    end
end