%% MD_solar_force.m

function F_solar = MD_solar_force(solar_pressure,A_proj,Cr,r_solar,shadow)
    % entradas
    % solar_pressure: presion solar a partir de la intensidad solar SF
    % A_proj: área proyectada en la dirección del solar
    % Cr: coeficiente de reflectividad
    % r_solar: vector de posición del sol respecto a la tierra UNITARIO

    %salidas
    % F_solar: fuerza solar en terna ECI
   
    if shadow == 0
        F_solar = -solar_pressure*Cr*A_proj*r_solar;
    elseif shadow == 1
        F_solar = [0;0;0];
    else
        error("unexpected value for shadow")
    end

end