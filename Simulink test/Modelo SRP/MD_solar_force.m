function F_solar = MD_solar_force(solar_pressure,A_proj,Cr,r_solar,shadow)
% --- MD_solar_force.m ---

% Calcula la fuerza por presión de radiación solar (SRP) mediante un modelo
% concentrado de placa equivalente: la magnitud se obtiene con la presión
% solar, el área proyectada y el coeficiente de reflectividad, mientras que
% la dirección se toma opuesta al vector Tierra→Sol como aproximación del 
% vector Satélite→Sol. La fuerza se anula automáticamente cuando el satélite
%  está en eclipse, usando "shadow" provisto por MD_shadow_function.

    % Entradas:
    % solar_pressure: presion solar a partir de la intensidad solar SF
    % A_proj: área proyectada en la dirección del solar
    % Cr: coeficiente de reflectividad
    % r_solar: vector de posición del sol respecto a la tierra UNITARIO
    % shadow: 1 si el satélite está en condiciones de eclipse, 0 si no lo está.

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