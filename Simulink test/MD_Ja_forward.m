function J_a_ECEF = MD_Ja_forward(R_ECEF, Sbar, Cbar, K,delta,GM,a0)
ex=[1;0;0]; ey=[0;1;0]; ez=[0;0;1];

[r, lambda, sinphi, cosphi, s, rhoP,phi] = MT_ECEF2GEOCENTRIC(R_ECEF+ delta*ex);
ap = MD_aceleracion_geopotencial(r, lambda, sinphi, cosphi, s, rhoP, Cbar, Sbar, K, GM,phi);
J_a_ECEF1 = (ap - a0) / delta;

[r, lambda, sinphi, cosphi, s, rhoP,phi] = MT_ECEF2GEOCENTRIC(R_ECEF+ delta*ey);
ap = MD_aceleracion_geopotencial(r, lambda, sinphi, cosphi, s, rhoP, Cbar, Sbar, K, GM,phi);
J_a_ECEF2 = (ap - a0) / delta;

[r, lambda, sinphi, cosphi, s, rhoP,phi] = MT_ECEF2GEOCENTRIC(R_ECEF+ delta*ez);
ap = MD_aceleracion_geopotencial(r, lambda, sinphi, cosphi, s, rhoP, Cbar, Sbar, K, GM,phi);
J_a_ECEF3 = (ap - a0) / delta;

J_a_ECEF=[J_a_ECEF1,J_a_ECEF2,J_a_ECEF3];
end
