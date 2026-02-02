%% DATA_sun.m
% Parámetros del modelo simplificado del Sol de Montenbruck para obtener la 
% dirección Sol–Tierra (o Sol–satélite) en un marco inercial.
% El modelo es lo suficientemente preciso para las necesidades del
% simulador para el período entre el año 1950 y el 2050.

solar_omegas  = 282.94;              % Ω+ω [deg] (parámetro de orientación usado en el modelo solar simplificado)
solar_epsilon = 23.43929111;         % ε [deg] oblicuidad de la eclíptica (inclinación del plano eclíptico respecto del ecuador)

% Pre-cálculo de cos(ε) y sin(ε) para evitar recomputarlo en cada paso
solar_cosep   = cosd(solar_epsilon); % cos(ε) 
solar_sinep   = sind(solar_epsilon); % sin(ε)