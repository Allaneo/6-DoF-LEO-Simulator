function J_a_ECEF = MD_Ja_forward(R_ECEF, Sbar, Cbar, K,GM,a0_ECEF)
ex=[1;0;0]; ey=[0;1;0]; ez=[0;0;1];
delta = 0.1;
% se transforma la entrada de ECEF a GEOCENTRIC para calcular la
% aceleración en el delta
[r, lambda, sinphi, cosphi, s, rhoP,phi] = MT_ECEF2GEOCENTRIC(R_ECEF+ delta*ex); 

% cálculo de la aceleración gravitatoria en el punto delta
ap = MD_aceleracion_geopotencial(r, lambda, sinphi, cosphi, s, rhoP, Cbar, Sbar, K, GM,phi);

% calculo del gradiente de gravedad en adelanto. Se podría mejorar la
% precisión utilizando diferencias centradas pero de esta forma se ahorra
% un calculo al reutilizar la gravedad del punto actual.
J_a_ECEF1 = (ap - a0_ECEF) / delta;

% se repite el proceso en cada dirección
[r, lambda, sinphi, cosphi, s, rhoP,phi] = MT_ECEF2GEOCENTRIC(R_ECEF+ delta*ey);
ap = MD_aceleracion_geopotencial(r, lambda, sinphi, cosphi, s, rhoP, Cbar, Sbar, K, GM,phi);
J_a_ECEF2 = (ap - a0_ECEF) / delta;

[r, lambda, sinphi, cosphi, s, rhoP,phi] = MT_ECEF2GEOCENTRIC(R_ECEF+ delta*ez);
ap = MD_aceleracion_geopotencial(r, lambda, sinphi, cosphi, s, rhoP, Cbar, Sbar, K, GM,phi);
J_a_ECEF3 = (ap - a0_ECEF) / delta;

J_a_ECEF=[J_a_ECEF1,J_a_ECEF2,J_a_ECEF3];
end
