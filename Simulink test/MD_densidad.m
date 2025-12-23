function rho = MD_densidad(r_eci)
rho0 = 1.454e-13;   % kg/m^3, densidad a 600 km
h0   = 600;     % km
H    = 71.835;      % km, altura de escala
h_km = norm(r_eci)/1000-6388;
%rho = rho0 * exp(-(h_km - h0)/H);
rho=rho0 * exp(-(h_km - h0)/H);
end