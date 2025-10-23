function rho = MD_densidad(r_eci)
rho0 = 1e-14;   % kg/m^3, densidad a 700 km
h0   = 700;     % km
H    = 60;      % km, altura de escala
h_km = norm(r_eci)/1000;
%rho = rho0 * exp(-(h_km - h0)/H);
rho=rho0 * exp(-(h_km - h0)/H);
end
