function a_ecef = MD_aceleracion_geopotencial(r, lambda, sinphi, cosphi, s, rhoP, Cbar, Sbar, K, GM, phi)
% MD_aceleracion_geopotencial
% Aceleración gravitatoria geopotencial en ECEF mediante armónicos esféricos.
% Convención: Pbar y derivadas vienen de MD_legendreN3_bar (allí se fija CS sí/no).
%
% Entradas:
%  r, phi, lambda, sinphi, cosphi, s, rhoP : del bloque MT_ECEF2GEOCENTRIC (o equivalente)
%  Cbar,Sbar,K : (N+1)x(N+1)
%  GM : mu [m^3/s^2]

% N máximo (grado/orden)
N = size(Cbar,1)-1;

% Legendre fully-normalized y derivadas w.r.t x = sinphi
% (Convención de fase CS ya resuelta adentro de MD_legendreN3_bar)
[Pbar, dPbar_dx, ~] = MD_legendreN3_bar(sinphi, cosphi, K);

% acumuladores
dU_dr_term      = 0;
dU_dphi_term    = 0;
dU_dlambda_term = 0;

% EGM geocéntrico: no hay términos de grado 1 -> arrancar en n=2
for n = 2:N
    for m = 0:n
        C = Cbar(n+1,m+1);
        S = Sbar(n+1,m+1);

        Pnm    = Pbar(n+1,m+1);
        dPnm_dx = dPbar_dx(n+1,m+1);

        cos_mlambda = cos(m * lambda);
        sin_mlambda = sin(m * lambda);

        F = (C * cos_mlambda + S * sin_mlambda);
        fac = rhoP(n+1) * F;

        % término radial (aparece (n+1) por d/dr de (GM/r)*(a/r)^n )
        dU_dr_term = dU_dr_term + (n + 1) * fac * Pnm;

        % d/dphi: x = sinphi => dx/dphi = cosphi
        dU_dphi_term = dU_dphi_term + fac * (dPnm_dx * cosphi);

        % d/dlambda
        if m ~= 0
            dU_dlambda_term = dU_dlambda_term + rhoP(n+1) * Pnm * ( -m * C * sin_mlambda + m * S * cos_mlambda );
        end
    end
end

% derivadas del potencial Φ
dPhi_dr      = -GM / (r^2) * ( 1 + dU_dr_term );
dPhi_dphi    =  GM / r     * ( dU_dphi_term );
dPhi_dlambda =  GM / r     * ( dU_dlambda_term );

% componentes geocéntricas esféricas (er, ephi, elambda)
a_r   = dPhi_dr;
a_phi = dPhi_dphi / r;

% protección de polos
if abs(cosphi) < 1e-12
    a_lambda = 0;
else
    a_lambda = dPhi_dlambda / (r * cosphi);
end

a_geo  = [a_r; a_phi; a_lambda];

% a ECEF usando tu módulo
a_ecef = MT_GEOCENTRIC2ECEF(a_geo, phi, lambda);

end
