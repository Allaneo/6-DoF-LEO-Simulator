%% MD_shadow_function.m
function shadow = MD_shadow_function(sun_position,r_eci)
% Determina si el satélite se encuentra en eclipse (sombra de la Tierra)
% bajo el modelo geométrico de sombra cilíndrica: se proyecta la posición 
% del satélite sobre la dirección Sol–Tierra y se evalúa si la línea 
% Sol–satélite intersecta el disco terrestre. El criterio separa el vector 
% posición del satélite en 
% (i) componente paralela a la dirección del Sol y 
% (ii) distancia perpendicular al eje Sol–Tierra, declarando eclipse cuando
% el satélite está "detrás" de la Tierra y dicha distancia es menor que el radio terrestre.

% Entradas:
% sun_position: vector posicion unitario del sol.
% r_eci: vector posicion del satelite

% Salida:
% shadow: boolean (1 or 0). 1 Si el satélite está en eclipse, 0 si no está
% en eclipse.

    F1 = r_eci'*sun_position;
    F2 = norm(r_eci-F1*sun_position);

    if  F1<0 && F2<6378163
        shadow = 1;
    else
        shadow = 0;
    end

end