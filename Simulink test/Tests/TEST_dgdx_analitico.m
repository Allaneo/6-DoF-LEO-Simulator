function TEST_dgdx_analitico()
% Valida dgdx analítico comparándolo con FD usando MD_aceleracion_geopotencial.

INIT_parametros;
Re = a_world;

N = size(Cbar,1)-1;

% puntos de test (ECEF) alejados de polos
alts = [200e3, 500e3, 900e3];                 % altitudes típicas
lats = deg2rad([0, 30, -45]);                 % geocentric lat (aprox)
lons = deg2rad([0, 60, 150]);                 % lon
R0   = Re + alts;

% step de FD
h = 5; % [m]

modes = [1 2 3]; % 1=Full (tu N=4), 2=J2, 3=Central

for mode = modes
    fprintf('\n=== gravity_mode = %d ===\n', mode);

    for k = 1:numel(R0)
        r_ecef = [ R0(k)*cos(lats(k))*cos(lons(k));
                   R0(k)*cos(lats(k))*sin(lons(k));
                   R0(k)*sin(lats(k)) ];

        % evita polos por las protecciones en tu aceleración
        if abs(r_ecef(3)/norm(r_ecef)) > 0.999
            continue;
        end

        % Analítico
        vars = ecef2vars(r_ecef, Re, N);
        dgdx_ana = MD_Ja_analitico( ...
            vars.r, vars.lambda, vars.sinphi, vars.cosphi, vars.s, vars.rhoP, ...
            Cbar, Sbar, K, GM, mode);

        % FD
        dgdx_fd = fd_jac(@(rr) accel_wrapper(rr, Re, N, Cbar, Sbar, K, GM, mode), r_ecef, h);

        % métricas
        abs_error = norm(dgdx_fd - dgdx_ana, 'fro');
        sym = norm(dgdx_ana - dgdx_ana.', 'fro');
        tr  = trace(dgdx_ana);

        fprintf('Caso %d: Error(FD)=%9.2e | sym=%9.2e | trace=%9.2e\n', k, abs_error, sym, tr);
    end
end

rvec = r_ecef;
r    = norm(rvec);
dgdx_ref = GM * ( 3*(rvec*rvec.')/r^5 - eye(3)/r^3 );

err_simple = norm(dgdx_ana - dgdx_ref,'fro') / norm(dgdx_ref,'fro');
fprintf("SIMPLE vs cerrado: %e\n", err_simple);

end

function J = fd_jac(fun, x, h)
J = zeros(3,3);
I = eye(3);
for i = 1:3
    fp = fun(x + h*I(:,i));
    fm = fun(x - h*I(:,i));
    J(:,i) = (fp - fm) / (2*h);
end
end

function a = accel_wrapper(r_ecef, Re, N, Cbar, Sbar, K, GM, mode)
% fuerza coherente con el mode:
% - Full: usa tu MD_aceleracion_geopotencial con Cbar/Sbar completos
% - J2:  arma Cbar/Sbar con solo C20
% - Central: usa fórmula cerrada (tu MD_aceleracion_geopotencial no es ideal para esto)

if mode == 3
    mu = GM;
    r  = norm(r_ecef);
    a = -mu * r_ecef / (r^3);
    return;
end

Cuse = Cbar; Suse = Sbar;
if mode == 2
    Cuse(:)=0; Suse(:)=0;
    Cuse(3,1) = Cbar(3,1); % C20
end

v = ecef2vars(r_ecef, Re, N);
a = MD_aceleracion_geopotencial(v.r, v.lambda, v.sinphi, v.cosphi, v.s, v.rhoP, Cuse, Suse, K, GM, v.phi);
end

function v = ecef2vars(r_ecef, Re, N)
x = r_ecef(1); y = r_ecef(2); z = r_ecef(3);

r = norm(r_ecef);
lambda = atan2(y,x);

rho = hypot(x,y);
phi = atan2(z, rho);      % geocentric latitude

sinphi = sin(phi);
cosphi = cos(phi);
s = cosphi;

% (Re/r)^n, n=0..N
rhoP = (Re/r).^(0:N).';

v = struct('r',r,'lambda',lambda,'phi',phi,'sinphi',sinphi,'cosphi',cosphi,'s',s,'rhoP',rhoP);
end
