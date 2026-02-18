function D = RP_add_derived(D, cfg, H)
% Derivados:
%  - shadowSegments
%  - F_grav_nonrad (ECI)
%  - actitud derivada (Euler, |omega|, nadir_err)

D.derived = struct();
D.derived.shadowSegments = [];

%% Shadow segments
if isfield(D.data,'shadow') && D.data.shadow.ok
    tt = D.data.shadow.t(:);
    s  = D.data.shadow.X(:) > 0.5;

    if numel(tt) == numel(s) && numel(tt) > 1
        ds = diff([false; s; false]);
        i1 = find(ds == 1);
        i2 = find(ds == -1) - 1;
        D.derived.shadowSegments = [tt(i1) tt(i2)];
    end
end

%% Gravedad no radial
D.data.F_grav_nonrad = struct('ok',false,'name','Fuerza Gravedad NO radial (ECI)','t',[],'X',[]);
if isfield(D.data,'F_grav') && D.data.F_grav.ok
    tFg = D.data.F_grav.t(:);
    Fg3 = H.vec3_timeseries_samplewise(D.data.F_grav.X, tFg);

    % r interpolada a tFg
    rx = interp1(D.tx, D.rECI(:,1), tFg, 'linear', 'extrap');
    ry = interp1(D.tx, D.rECI(:,2), tFg, 'linear', 'extrap');
    rz = interp1(D.tx, D.rECI(:,3), tFg, 'linear', 'extrap');
    rI = [rx ry rz];

    rnorm = sqrt(sum(rI.^2,2));
    rhat  = rI ./ [rnorm rnorm rnorm];

    Frad_mag = sum(Fg3 .* rhat, 2);
    F_rad = rhat .* [Frad_mag Frad_mag Frad_mag];
    F_perp = Fg3 - F_rad;

    D.data.F_grav_nonrad.ok = true;
    D.data.F_grav_nonrad.t  = tFg;
    D.data.F_grav_nonrad.X  = F_perp;
end

%% Actitud derivada
D.derived.att = struct('ok',false);

if D.haveAtt
    q_raw = D.q_raw;
    wB    = D.wB;
    tA    = D.tA;

    qn = sqrt(sum(q_raw.^2,2));
    omega_mag = sqrt(sum(wB.^2,2));

    q_plot = H.quat_continuity(q_raw);

    qn_safe = qn; qn_safe(qn_safe < 1e-12) = 1;
    q_used = q_plot ./ [qn_safe qn_safe qn_safe qn_safe];

    [c11,c12,c13,c21,c22,c23,c31,c32,c33] = H.quat_to_dcm_cib(q_used);
    [yaw, pitch, roll] = H.dcm_to_euler_zyx(c11,c21,c31,c32,c33);

    yaw_deg   = H.unwrap_deg(yaw*180/pi);
    pitch_deg = H.unwrap_deg(pitch*180/pi);
    roll_deg  = H.unwrap_deg(roll*180/pi);

    % Eje body a inercial y error a nadir
    axb = cfg.body_axis_for_nadir(:);
    axb = axb / max(norm(axb), eps);

    axI = [c11*axb(1) + c12*axb(2) + c13*axb(3), ...
           c21*axb(1) + c22*axb(2) + c23*axb(3), ...
           c31*axb(1) + c32*axb(2) + c33*axb(3)];

    rx = interp1(D.tx, D.rECI(:,1), tA, 'linear', 'extrap');
    ry = interp1(D.tx, D.rECI(:,2), tA, 'linear', 'extrap');
    rz = interp1(D.tx, D.rECI(:,3), tA, 'linear', 'extrap');
    rI = [rx ry rz];
    rnorm = sqrt(sum(rI.^2,2));
    rhat = rI ./ [rnorm rnorm rnorm];

    nadir = -rhat;
    cosang = sum(axI .* nadir, 2);
    cosang = max(min(cosang,1),-1);
    nadir_err = acos(cosang) * 180/pi;

    D.derived.att.ok = true;
    D.derived.att.tA = tA;
    D.derived.att.yaw_deg = yaw_deg;
    D.derived.att.pitch_deg = pitch_deg;
    D.derived.att.roll_deg = roll_deg;
    D.derived.att.wB = wB;
    D.derived.att.omega_mag = omega_mag;
    D.derived.att.q_norm_err = qn - 1;
    D.derived.att.nadir_err = nadir_err;
end
end
