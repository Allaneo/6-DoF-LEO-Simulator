function H = RP_helpers()
% Biblioteca de helpers accesible desde módulos (sin duplicar código).
H = struct();

H.get_param_value = @get_param_value;

H.get_dataset_timeseries_exact = @get_dataset_timeseries_exact;
H.get_sdi_timeseries = @get_sdi_timeseries;
H.get_sdi_timeseries_by_name = @get_sdi_timeseries_by_name;
H.read_ts = @read_ts;

H.mag_timeseries_samplewise = @mag_timeseries_samplewise;
H.vec3_timeseries_samplewise = @vec3_timeseries_samplewise;

H.apply_shadow = @apply_shadow;

H.ecef2lla_wgs84_deg = @ecef2lla_wgs84_deg;

H.quat_continuity = @quat_continuity;
H.quat_to_dcm_cib = @quat_to_dcm_cib;
H.dcm_to_euler_zyx = @dcm_to_euler_zyx;
H.unwrap_deg = @unwrap_deg;

H.coe_from_rv = @coe_from_rv;

H.era_from_jd_ut1 = @era_from_jd_ut1;
end

function val = get_param_value(x)
% Si es Simulink.Parameter devuelve .Value, si no devuelve x.
try
    if isa(x,'Simulink.Parameter')
        val = x.Value;
    else
        val = x;
    end
catch
    val = x;
end
end

function ts = get_dataset_timeseries_exact(logs, name)
% Importante: NO usar logs.get(name) porque dispara warnings si no existe.
ts = [];
if isempty(logs), return; end

try
    nEl = logs.numElements;
catch
    return;
end

% 1) match exacto (case-sensitive)
for i = 1:nEl
    try
        nm = logs{i}.Name;
    catch
        nm = '';
    end
    if ~isempty(nm) && strcmp(nm, name)
        try
            ts = logs{i}.Values;
        catch
            ts = [];
        end
        return;
    end
end

% 2) match case-insensitive
for i = 1:nEl
    try
        nm = logs{i}.Name;
    catch
        nm = '';
    end
    if ~isempty(nm) && strcmpi(nm, name)
        try
            ts = logs{i}.Values;
        catch
            ts = [];
        end
        return;
    end
end
end

function ts = get_sdi_timeseries(sigObj)
ts = [];
try, ts = sigObj.Values; return; catch, end
try, ts = sigObj.getValues; return; catch, end
end

function ts = get_sdi_timeseries_by_name(sdiNames, sdiTS, name)
ts = [];
if isempty(sdiNames), return; end

idx = find(strcmp(sdiNames, name), 1, 'first');
if ~isempty(idx) && ~isempty(sdiTS{idx}), ts = sdiTS{idx}; return; end

idx = find(strcmpi(sdiNames, name), 1, 'first');
if ~isempty(idx) && ~isempty(sdiTS{idx}), ts = sdiTS{idx}; return; end
end

function [t, X] = read_ts(ts)
t = ts.Time(:);
Xraw = squeeze(ts.Data);

% Asegurar NxM con N = numel(t)
if size(Xraw,1) == numel(t)
    X = Xraw;
elseif size(Xraw,2) == numel(t)
    X = Xraw.';
else
    % fallback: intentar reshape si es compatible
    try
        X = reshape(Xraw, numel(t), []);
    catch
        X = Xraw;
    end
end

if isrow(X), X = X(:); end
end

function mag = mag_timeseries_samplewise(X, t)
% Devuelve magnitud por muestra.
% Si X ya es escalar => abs(X) (asumo que es magnitud o componente).
% Si X tiene >=3 columnas => norma de las 3 primeras.
% Si tiene 2 columnas => norma 2D.
n = numel(t);
Xs = squeeze(X);

if size(Xs,1) ~= n && size(Xs,2) == n
    Xs = permute(Xs, [2 1 3:ndims(Xs)]);
end

if size(Xs,1) ~= n
    Xs = reshape(Xs, n, []);
else
    Xs = reshape(Xs, n, []);
end

m = size(Xs,2);

if m == 1
    mag = abs(Xs(:,1));
elseif m >= 3
    mag = sqrt(sum(Xs(:,1:3).^2, 2));
else
    mag = sqrt(sum(Xs(:,1:m).^2, 2));
end
end

function X3 = vec3_timeseries_samplewise(X, t)
% Intenta devolver Nx3. Si no puede (porque la señal es escalar o <3),
% devuelve [] (NO error), así el caller decide qué hacer.
n = numel(t);
Xs = squeeze(X);

if size(Xs,1) ~= n && size(Xs,2) == n
    Xs = permute(Xs, [2 1 3:ndims(Xs)]);
end

if size(Xs,1) ~= n
    Xs = reshape(Xs, n, []);
else
    Xs = reshape(Xs, n, []);
end

if size(Xs,2) < 3
    X3 = [];
    return;
end

X3 = Xs(:,1:3);
end

function apply_shadow(ax, segs)
if isempty(segs), return; end
yl = get(ax,'YLim');

for r = 1:size(segs,1)
    patch(ax, ...
        [segs(r,1) segs(r,2) segs(r,2) segs(r,1)], ...
        [yl(1) yl(1) yl(2) yl(2)], ...
        [0.85 0.85 0.85], ...
        'EdgeColor','none', 'FaceAlpha',0.25, 'HandleVisibility','off');
end
uistack(findall(ax,'Type','line'),'top');
end

function [lat_deg, lon_deg] = ecef2lla_wgs84_deg(x, y, z)
a  = 6378137.0;
f  = 1/298.257223563;
e2 = f*(2-f);
b  = a*sqrt(1-e2);
ep2 = (a^2 - b^2)/b^2;

lon = atan2(y, x);
p   = sqrt(x.^2 + y.^2);

theta = atan2(z*a, p*b);
st = sin(theta);
ct = cos(theta);

lat = atan2(z + ep2*b*st.^3, p - e2*a*ct.^3);

lat_deg = lat * 180/pi;
lon_deg = lon * 180/pi;
end

function q = quat_continuity(q)
for i = 2:size(q,1)
    if sum(q(i,:).*q(i-1,:)) < 0
        q(i,:) = -q(i,:);
    end
end
end

function [c11,c12,c13,c21,c22,c23,c31,c32,c33] = quat_to_dcm_cib(q)
% C_IB (BODY->INERTIAL) desde cuaternión scalar-first.
q0 = q(:,1); q1 = q(:,2); q2 = q(:,3); q3 = q(:,4);

c11 = q0.^2 + q1.^2 - q2.^2 - q3.^2;
c12 = 2*(q1.*q2 - q0.*q3);
c13 = 2*(q1.*q3 + q0.*q2);

c21 = 2*(q1.*q2 + q0.*q3);
c22 = q0.^2 - q1.^2 + q2.^2 - q3.^2;
c23 = 2*(q2.*q3 - q0.*q1);

c31 = 2*(q1.*q3 - q0.*q2);
c32 = 2*(q2.*q3 + q0.*q1);
c33 = q0.^2 - q1.^2 - q2.^2 + q3.^2;
end

function [yaw, pitch, roll] = dcm_to_euler_zyx(c11,c21,c31,c32,c33)
s = -c31;
s = max(min(s,1),-1);
pitch = asin(s);
roll  = atan2(c32, c33);
yaw   = atan2(c21, c11);
end

function a_deg = unwrap_deg(a_deg)
for i = 2:numel(a_deg)
    da = a_deg(i) - a_deg(i-1);
    if da > 180
        a_deg(i:end) = a_deg(i:end) - 360;
    elseif da < -180
        a_deg(i:end) = a_deg(i:end) + 360;
    end
end
end

function [a, e, inc, raan, argp, nu] = coe_from_rv(rI, vI, mu)
tol_n = 1e-12;
tol_e = 1e-12;

r = rI(:); v = vI(:);
R = norm(r);
V = norm(v);

h = cross(r, v);
hnorm = norm(h);

k = [0;0;1];
n = cross(k, h);
nnorm = norm(n);

evec = (1/mu)*((V^2 - mu/R)*r - (dot(r,v))*v);
e = norm(evec);

energy = V^2/2 - mu/R;
if abs(energy) > 1e-16
    a = -mu/(2*energy);
else
    a = Inf;
end

if hnorm < tol_n
    inc = 0;
else
    c = h(3)/hnorm;
    c = max(min(c,1),-1);
    inc = acos(c);
end

if nnorm < tol_n
    raan = 0;
else
    raan = atan2(n(2), n(1));
    if raan < 0, raan = raan + 2*pi; end
end

if e < tol_e
    argp = 0;
    if nnorm < tol_n
        nu = atan2(r(2), r(1));
        if nu < 0, nu = nu + 2*pi; end
    else
        sin_u = dot(cross(n, r), h) / (nnorm*R*hnorm);
        cos_u = dot(n, r) / (nnorm*R);
        nu = atan2(sin_u, cos_u);
        if nu < 0, nu = nu + 2*pi; end
    end
else
    sin_nu = dot(cross(evec, r), h) / (e*R*hnorm);
    cos_nu = dot(evec, r) / (e*R);
    nu = atan2(sin_nu, cos_nu);
    if nu < 0, nu = nu + 2*pi; end

    if nnorm < tol_n
        argp = atan2(evec(2), evec(1));
        if argp < 0, argp = argp + 2*pi; end
        raan = 0;
    else
        sin_w = dot(cross(n, evec), h) / (nnorm*e*hnorm);
        cos_w = dot(n, evec) / (nnorm*e);
        argp = atan2(sin_w, cos_w);
        if argp < 0, argp = argp + 2*pi; end
    end
end
end

function theta = era_from_jd_ut1(jd_ut1_0, tsec)
jd_ut1 = jd_ut1_0 + tsec(:)/86400;
fraction = jd_ut1 - 2451545.0;
rev = 0.7790572732640 + 1.00273781191135448 * fraction;
rev_frac = rev - floor(rev);
theta = 2*pi*rev_frac;
end
