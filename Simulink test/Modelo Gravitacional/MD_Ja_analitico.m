function dgdx_ECEF = MD_Ja_analitico(r_ecef, a_world, Cbar, Sbar, K, GM, GRAV_MODEL)
% MD_Ja_analitico
% Jacobiano dg/dx (3x3) de la aceleración gravitatoria en ECEF respecto de r_ECEF.
% Modelos: GRAV_MODEL = "SIMPLE" | "J2" | "N4"

X = r_ecef(1);  Y = r_ecef(2);  Z = r_ecef(3);
r = norm(r_ecef);

% --- Caso SIMPLE: a = -mu r / r^3 => dg/dx = mu*(3 rr^T/r^5 - I/r^3)
if GRAV_MODEL == 3
    mu = GM;
    r2 = r*r;
    r3 = r2*r;
    r5 = r3*r2;
    dgdx_ECEF = mu*(3*(r_ecef*r_ecef.')/r5 - eye(3)/r3);
    return;
end

% --- Geometría geocéntrica
rho_xy = hypot(X, Y);
lambda = atan2(Y, X);
phi    = atan2(Z, rho_xy);     % lat geocéntrica

sinphi = Z / r;
cosphi = rho_xy / r;
s      = cosphi;

% protecciones
eps_s = 1e-12;
s_safe = max(abs(s), eps_s);
rho2 = X*X + Y*Y;

% --- Potencias rhoP = (a/r)^n, n=0..Nmax
Nmax = size(Cbar,1) - 1;     % con 4x4 -> 3
rhoP = (a_world/r).^(0:Nmax).';

% --- Switch explícito J2
onlyJ2 = (GRAV_MODEL == 2);

% --- Legendre fully-normalized y derivadas (x = sinphi)
[Pbar, dPdx, d2Pdx2] = MD_legendreN3_bar(sinphi, cosphi, K);

% --- Sumas (bloque "núcleo")
S_r              = 0;
S_rr             = 0;
S_r_x            = 0;
S_r_lambda       = 0;

S_phi_x          = 0;
S_phi_x_n        = 0;
S_phi_x_x        = 0;
S_phi_x_lambda   = 0;

S_lambda         = 0;
S_lambda_n       = 0;
S_lambda_x       = 0;
S_lambda_lambda  = 0;

for n = 2:Nmax
    rho = rhoP(n+1);

    for m = 0:n
        if onlyJ2 && ~(n==2 && m==0)
            continue;
        end

        C = Cbar(n+1,m+1);
        S = Sbar(n+1,m+1);

        cml = cos(m*lambda);
        sml = sin(m*lambda);

        f    = C*cml + S*sml;                 % combinación trig
        fl   = m*(-C*sml + S*cml);            % d/dlambda de f
        fll  = -m*m*f;                        % d2/dlambda2 de f

        P    = Pbar(n+1,m+1);
        dP   = dPdx(n+1,m+1);
        d2P  = d2Pdx2(n+1,m+1);

        fac   = rho * f;
        fac_l = rho * fl;

        % términos para a_r
        S_r        = S_r        + (n+1)*fac*P;
        S_rr       = S_rr       + n*(n+1)*fac*P;
        S_r_x      = S_r_x      + (n+1)*fac*dP;
        S_r_lambda = S_r_lambda + (n+1)*fac_l*P;

        % términos para a_phi (trabajando con S_phi_x y luego S_phi = s*S_phi_x)
        S_phi_x        = S_phi_x        + fac*dP;
        S_phi_x_n      = S_phi_x_n      + n*fac*dP;
        S_phi_x_x      = S_phi_x_x      + fac*d2P;
        S_phi_x_lambda = S_phi_x_lambda + fac_l*dP;

        % términos para a_lambda
        S_lambda        = S_lambda        + fac_l*P;
        S_lambda_n      = S_lambda_n      + n*fac_l*P;
        S_lambda_x      = S_lambda_x      + fac_l*dP;
        S_lambda_lambda = S_lambda_lambda + rho*fll*P;
    end
end

% --- a_sph (geocéntrico)
Kgrav = GM / (r*r);

S_phi = s * S_phi_x;

a_r   = -Kgrav*(1 + S_r);
a_phi =  Kgrav*S_phi;
a_lam =  Kgrav*(S_lambda / s_safe);

a_sph = [a_r; a_phi; a_lam];

% --- Derivadas parciales de a_sph respecto de (r, x=sinphi, lambda)
dS_r_dr     = -(S_rr)/r;
dS_phi_dr   = s * (-(S_phi_x_n)/r);
dS_lam_dr   = -(S_lambda_n)/r;

dS_r_dx     = S_r_x;
dS_phi_dx   = (-(sinphi/s_safe))*S_phi_x + s*S_phi_x_x;
dS_lam_dx   = S_lambda_x;

dS_r_dlam    = S_r_lambda;
dS_phi_dlam  = s*S_phi_x_lambda;
dS_lam_dlam  = S_lambda_lambda;

dK_dr = -2*Kgrav/r;

da_r_dr   = (2*Kgrav/r)*(1+S_r) - Kgrav*dS_r_dr;
da_r_dx   = -Kgrav*dS_r_dx;
da_r_dlam = -Kgrav*dS_r_dlam;

da_phi_dr   = dK_dr*S_phi + Kgrav*dS_phi_dr;
da_phi_dx   = Kgrav*dS_phi_dx;
da_phi_dlam = Kgrav*dS_phi_dlam;

da_lam_dr   = dK_dr*(S_lambda/s_safe) + Kgrav*(dS_lam_dr/s_safe);
da_lam_dx   = Kgrav*(dS_lam_dx/s_safe + S_lambda*(sinphi)/(s_safe^3));
da_lam_dlam = Kgrav*(dS_lam_dlam/s_safe);

da_dr   = [da_r_dr;   da_phi_dr;   da_lam_dr];
da_dx   = [da_r_dx;   da_phi_dx;   da_lam_dx];
da_dlam = [da_r_dlam; da_phi_dlam; da_lam_dlam];

% --- Jacobiano de (r, x=sinphi, lambda) respecto de (X,Y,Z)
r2 = r*r;
r3 = r2*r;

drdX = X/r;  drdY = Y/r;  drdZ = Z/r;

dxdX = -Z*X/r3;
dxdY = -Z*Y/r3;
dxdZ = (rho2)/r3;

if rho2 < 1e-24
    dlmdX = 0; dlmdY = 0; dlmdZ = 0;   % singular en el eje Z
else
    dlmdX = -Y/rho2;
    dlmdY =  X/rho2;
    dlmdZ =  0;
end

da_sph_dX = da_dr*drdX + da_dx*dxdX + da_dlam*dlmdX;
da_sph_dY = da_dr*drdY + da_dx*dxdY + da_dlam*dlmdY;
da_sph_dZ = da_dr*drdZ + da_dx*dxdZ;

% --- Derivadas de la transformación (dE/dX, etc.) para llevar a ECEF
cl = cos(lambda);
sl = sin(lambda);

% d/dx (x=sinphi) de E (columnas er,ephi,elam)
dsdx = -(sinphi/s_safe);

E_x = [ dsdx*cl,   -cl,   0;
        dsdx*sl,   -sl,   0;
        1,         dsdx,  0 ];

% d/dlambda de E
E_lam = [ -s*sl,    sinphi*sl,  -cl;
           s*cl,   -sinphi*cl,  -sl;
           0,       0,           0 ];

E_dX = E_x*dxdX + E_lam*dlmdX;
E_dY = E_x*dxdY + E_lam*dlmdY;
E_dZ = E_x*dxdZ;

% --- dg/dx en ECEF
dgdx_ECEF = zeros(3,3);

dgdx_ECEF(:,1) = E_dX*a_sph + MT_GEOCENTRIC2ECEF(da_sph_dX, phi, lambda);
dgdx_ECEF(:,2) = E_dY*a_sph + MT_GEOCENTRIC2ECEF(da_sph_dY, phi, lambda);
dgdx_ECEF(:,3) = E_dZ*a_sph + MT_GEOCENTRIC2ECEF(da_sph_dZ, phi, lambda);

end
