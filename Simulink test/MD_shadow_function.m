%% MD_shadow_function.m
function shadow = MD_shadow_function(sun_position,r_eci)
% entradas
% sun_position: vector posicion unitario del sol.
% r_eci: vector posicion del satelite

    F1 = r_eci'*sun_position;
    F2 = norm(r_eci-F1*sun_position);

    if  F1<0 && F2<6378163
        shadow = 1;
    else
        shadow = 0;
    end

end