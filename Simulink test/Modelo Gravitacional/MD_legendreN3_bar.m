function [Pbar, dPdx, d2Pdx2] = MD_legendreN3_bar(sinphi, cosphi, K)
% MD_legendreN3_bar
% Pbar, dPdx, d2Pdx2 para n<=3, m<=n, con x = sin(phi)=cos(theta).
% CONVENCIÓN FIJA: NO incluye fase Condon–Shortley en el resultado final.
% (equivalente a tomar la definición con CS y multiplicar por (-1)^m).
%
% Entradas:
%   sinphi : x = sin(phi)
%   cosphi : s = cos(phi) (geocéntrica, >=0 fuera de polos)
%   K      : matriz (N+1)x(N+1) con factores de normalización (N<=3)
%
% Salidas:
%   Pbar, dPdx, d2Pdx2: (N+1)x(N+1) con índices (n+1,m+1)

Nmax = size(K,1) - 1;
if Nmax > 3
    error('MD_legendreN3_bar: Nmax=%d > 3. Rutina hardcodeada para N<=3.', Nmax);
end

x = sinphi;
s = cosphi;

% protección numérica cerca de polos
eps_s  = 1e-12;
s_safe = max(s, eps_s);
s3     = s_safe^3;

U      = zeros(Nmax+1, Nmax+1);
dUdx   = zeros(Nmax+1, Nmax+1);
d2Udx2 = zeros(Nmax+1, Nmax+1);

% -------------------------
% Base "con CS" (tipo física)
% -------------------------
% n=0
U(1,1) = 1;

% n=1
if Nmax >= 1
    U(2,1)      = x;
    dUdx(2,1)   = 1;

    U(2,2)      = -s;
    dUdx(2,2)   = x / s_safe;
    d2Udx2(2,2) = 1 / s3;
end

% n=2
if Nmax >= 2
    U(3,1)      = 0.5*(3*x^2 - 1);
    dUdx(3,1)   = 3*x;
    d2Udx2(3,1) = 3;

    U(3,2)      = -3*x*s;
    dUdx(3,2)   = -3*(1 - 2*x^2) / s_safe;
    d2Udx2(3,2) = 3*x*(3 - 2*x^2) / s3;

    U(3,3)      = 3*s^2;
    dUdx(3,3)   = -6*x;
    d2Udx2(3,3) = -6;
end

% n=3
if Nmax >= 3
    U(4,1)      = 0.5*(5*x^3 - 3*x);
    dUdx(4,1)   = 0.5*(15*x^2 - 3);
    d2Udx2(4,1) = 15*x;

    U(4,2)      = 1.5*(1 - 5*x^2)*s;
    dUdx(4,2)   = -(3*x*(11 - 15*x^2)) / (2*s_safe);
    d2Udx2(4,2) = -(3*(30*x^4 - 45*x^2 + 11)) / (2*s3);

    U(4,3)      = 15*x*s^2;
    dUdx(4,3)   = 15*(1 - 3*x^2);
    d2Udx2(4,3) = -90*x;

    U(4,4)      = -15*s^3;
    dUdx(4,4)   = 45*x*s;
    d2Udx2(4,4) = 45*(1 - 2*x^2) / s_safe;
end

% -------------------------
% Normalización
% -------------------------
Pbar   = K .* U;
dPdx   = K .* dUdx;
d2Pdx2 = K .* d2Udx2;

% -------------------------
% Convención FINAL: NO-CS  -> multiplicar por (-1)^m
% -------------------------
phase_m = (-1).^(0:Nmax);     % m = 0..Nmax
Pbar    = bsxfun(@times, Pbar,    phase_m);
dPdx    = bsxfun(@times, dPdx,    phase_m);
d2Pdx2  = bsxfun(@times, d2Pdx2,  phase_m);

end
