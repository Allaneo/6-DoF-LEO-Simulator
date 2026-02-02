function v_rel_I = MD_velocidad_relativa(v_I,r_I)
% --- MD_velocidad_relativa ---
% Calcula la velocidad relativa del satélite respecto de una atmósfera
% co-rotante con la Tierra en el marco inercial, restando a la velocidad
% inercial del satélite la velocidad inducida por la rotación terrestre. 
% Este modelo corresponde a la aproximación estándar de atmósfera 
% "solidaria" a la Tierra para el cálculo de drag.

% Entradas:
% v_I: velocidad en terna inercial
% r_I: posición en terna inercial

omega_earth = [0;0;7.2921150e-5]; %velocidad de rotación de la Tierra.

% 1) velocidad relativa en ECI (atmósfera co-rotante)
v_rel_I = v_I - cross(omega_earth, r_I);
end