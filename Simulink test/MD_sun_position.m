function longitud_solar = MD_sun_position(M,solar_omegas)
% longitud eclíptica del sol
longitud_solar = solar_omegas+M+(6892/360)*sind(M)+(72/360)*sind(2*M);
end