function tau_gg_simple = MD_torque_grav_simple(r_body, I_body,GM)
% código de verificación del modelo más complejo utilizado en la
% simulación. Este código NO influye en los resultados de la simulación.

% Entradas:
%   r_body: Vector posición del satélite expresado en marco BODY [m] (3x1)
%   I_body: Tensor de inercia [kg*m^2] (3x3)

    R_mag = norm(r_body);         % Distancia escalar
    
    % Vector unitario en Body (Nadir invertido local)
    u_body = r_body / R_mag; 

    % Fórmula del gradiente de gravedad
    tau_gg_simple = (3 * GM / R_mag^3) * cross(u_body, I_body * u_body);    
end