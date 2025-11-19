function a_ecef = MD_aceleracion_geopotencial(r, lambda, sinphi, cosphi, s, rhoP, Cbar, Sbar, K, GM, phi)
% Calcula componentes de aceleración en el sistema geocéntrico (a_r,a_phi,a_lambda).
% Entradas:
%  r, phi, lambda, sinphi, cosphi, s, rhoP : provenientes del bloque anterior (MT_ECEF2GEOCENTRIC)
%  Cbar,Sbar,K : matrices (N+1 x N+1) (4x4)

% umbral para evitar divisiones por s cero
eps_s = 1e-8;
s_safe = max(s, eps_s); % toma el valor mas grande para evitar tomar s en caso de que sea cero

% convenciones
N = size(Cbar,1)-1;

% --- construir U_nm y dU/dx
U = zeros(N+1,N+1);
dUdx = zeros(N+1,N+1);

%% Coeficientes de normalización U y derivadas du/dphi
% n=0
U(1,1)=1; dUdx(1,1)=0;
% n=1
  U(2,1)= sinphi;                 dUdx(2,1)= 1;
  U(2,2)= -s;                dUdx(2,2)= sinphi / s_safe;
% n=2
  U(3,1)= (3*sinphi^2 - 1)/2;     dUdx(3,1)= 3*sinphi;
  U(3,2)= -3*sinphi*s;            dUdx(3,2)= -3*(1 - 2*sinphi^2) / s_safe;
  U(3,3)= 3*s^2;             dUdx(3,3)= -6*sinphi;
% n=3
  U(4,1)= (5*sinphi^3 - 3*sinphi)/2;   dUdx(4,1)= (15*sinphi^2 - 3)/2;
  U(4,2)= (3/2)*(1 - 5*sinphi^2)*s; dUdx(4,2)= - (3*sinphi*(11 - 15*sinphi^2)) / (2*s_safe);
  U(4,3)= 15*sinphi*s^2;          dUdx(4,3)= 15*(1 - 3*sinphi^2);
  U(4,4)= -15*s^3;           dUdx(4,4)= 45*sinphi*s;

% aplicar K
Pbar = K .* U;
dPbar_dx = K .* dUdx;

% acumuladores
dU_dr_term = 0;
dU_dphi_term = 0;
dU_dlambda_term = 0;

for n = 1:N
    for m = 0:n
        C = Cbar(n+1,m+1); S = Sbar(n+1,m+1);
        Pnm = Pbar(n+1,m+1);
        dPnm_dx = dPbar_dx(n+1,m+1);

        cos_mlambda = cos(m * lambda);
        sin_mlambda = sin(m * lambda);
        factor = rhoP(n+1) * ( C * cos_mlambda + S * sin_mlambda );

        dU_dr_term = dU_dr_term + (n + 1) * factor * Pnm;
        % dx/dphi = cos(phi) with sign -> usar cosphi (firmado)
        dU_dphi_term = dU_dphi_term + factor * (dPnm_dx * cosphi);
        if m ~= 0
            dU_dlambda_term = dU_dlambda_term + rhoP(n+1) * Pnm * ( - m * C * sin_mlambda + m * S * cos_mlambda );
        end
    end
end

% derivadas parciales del potencial Phi
dPhi_dr      = -GM / (r^2) * ( 1 + dU_dr_term );
dPhi_dphi    =  GM / r   * ( dU_dphi_term );
dPhi_dlambda =  GM / r   * ( dU_dlambda_term );

% componentes esféricas
a_r = dPhi_dr;
a_phi = dPhi_dphi / r;  % 1/r factor
if abs(cosphi) < 1e-12
    a_lambda = 0;       % protección en polos
else
    a_lambda = dPhi_dlambda / (r * cosphi);
end
a_geo = [a_r;a_phi;a_lambda];
a_ecef = MT_GEOCENTRIC2ECEF(a_geo, phi, lambda);
end
