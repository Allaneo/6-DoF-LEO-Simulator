function [Omega_dot, tau_rw] = MD_RW(u_cmd, Omega, Jw, u_peak, Omega_max, P_peak)
% MD_RW
% Rueda de inercia ideal 1-eje (torque comandado) con límites:
%   - límite de torque (u_peak)
%   - límite de velocidad (Omega_max)
%   - derating por potencia (P_peak) opcional
%
% ENTRADAS:
%   u_cmd     [N*m]    torque requerido
%   Omega     [rad/s]  velocidad actual (estado, integrado afuera)
%   Jw        [kg*m^2] inercia de la rueda
%   u_peak    [N*m]    torque máximo (Inf desactiva)
%   Omega_max [rad/s]  velocidad máxima (Inf desactiva)
%   P_peak    [W]      potencia mecánica máxima (Inf desactiva)
%
% SALIDAS:
%   Omega_dot [rad/s^2] aceleración de la rueda
%   tau_rw    [3x1 N*m] torque aplicado por la rueda SOBRE EL SATÉLITE (reacción)
%
% CONVENCIÓN:
%   tau_rw(1) = -u_act  (eje X body). Los otros ejes = 0.

% ---- forzar escalares (por si entra vector)
u_cmd = u_cmd(1);
Omega = Omega(1);

% ---- preasignar salidas (tamaños fijos SIEMPRE)
Omega_dot = 0;
tau_rw    = zeros(3,1);

% ---- sanitizar entradas y parámetros
if ~isfinite(u_cmd);    u_cmd = 0; end
if ~isfinite(Omega);    Omega = 0; end

if ~isfinite(Jw) || Jw <= 0
    return; % deja salidas en cero
end

if ~isfinite(u_peak) || u_peak < 0
    u_peak = 0;
end

% ---- 1) límite efectivo por potencia (opcional)
u_lim = u_peak;

if isfinite(P_peak)
    Omega_eps   = 1e-3; % evita singularidad
    denom       = max(abs(Omega), Omega_eps);
    u_lim_power = P_peak / denom;
    u_lim       = min(u_lim, u_lim_power);
end

% si u_lim quedó raro
if ~isfinite(u_lim) || u_lim < 0
    u_lim = 0;
end

% ---- 2) saturación por torque
u_act = max(min(u_cmd, u_lim), -u_lim);

% ---- 3) saturación por velocidad (no acelerar más allá de Omega_max)
if isfinite(Omega_max)
    if (Omega >=  Omega_max && u_act > 0) || (Omega <= -Omega_max && u_act < 0)
        u_act = 0;
    end
end

% ---- 4) dinámica + torque de reacción
Omega_dot = u_act / Jw;
tau_rw(1) = -u_act;

end
