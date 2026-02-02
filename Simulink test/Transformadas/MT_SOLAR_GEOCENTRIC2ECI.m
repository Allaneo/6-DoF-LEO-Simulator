function r_solar = MT_SOLAR_GEOCENTRIC2ECI(lambda,solar_cosep,solar_sinep)
% posición del sol en ECEF
r_solar = [cosd(lambda);
    sind(lambda)*solar_cosep;
    sind(lambda)*solar_sinep];
r_solar = r_solar/norm(r_solar);
end