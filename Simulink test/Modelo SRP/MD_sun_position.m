function longitud_solar = MD_sun_position(M,solar_omegas)
% Calcula la longitud eclíptica del Sol a partir de los elementos orbitales
%  medios del modelo simplificado: se combina la longitud media con 
% correcciones periódicas de primer y segundo armónico en función de la 
% anomalía media.Esta magnitud se usa luego para obtener la dirección 
% Sol–Tierra en el marco inercial según el modelo astronómico adoptado en 
% el simulador.

%Entradas:
% M: anomalía media del Sol (en grados, consistente con sind).
% solar_omegas: términos iniciales del modelo (en grados).

% Salida:
% longitud_solar: longitud eclíptica (en grados).
longitud_solar = solar_omegas+M+(6892/360)*sind(M)+(72/360)*sind(2*M);
end