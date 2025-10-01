function a_drag = AtmosphericDrag(r_ECI, v_ECI, Cd, A_m, rho)
% AtmosphericDrag  Drag acceleration in ECI via LVLH projection
% 
% Inputs:
%   r_ECI  : 3×1 posición satélite en ECI [m]
%   v_ECI  : 3×1 velocidad satélite en ECI [m/s]
%   Cd     : coeficiente de arrastre (adimensional)
%   A_m    : razón área/masa del satélite [m²/kg]
%   rho    : densidad atmosférica local [kg/m³]
%
% Output:
%   a_drag : 3×1 aceleración por arrastre en ECI [m/s²]

%#codegen

% 1. Velocidad del viento atmosférico en ECI (rotación terrestre)
Omega_Earth = 7.2921150e-5;           % velocidad angular de la Tierra [rad/s]
omega      = [0; 0; Omega_Earth];    

v_atm = cross(omega, r_ECI);          
v_rel = v_ECI - v_atm;                
v_rel_norm = norm(v_rel);

if v_rel_norm < 1e-8
    a_drag = zeros(3,1);
    return;
end

% 2. Base LVLH: radial (er), normal al plano (h), along-track (etheta)
r_norm = norm(r_ECI);
er     = r_ECI / r_norm;

h_vec  = cross(r_ECI, v_ECI);
h_norm = norm(h_vec);

if h_norm < 1e-8
    % Órbita degenerada: arrastre opuesto a v_rel
    drag_dir = -v_rel / v_rel_norm;
else
    h      = h_vec / h_norm;
    etheta = cross(h, er);
    drag_dir = -etheta;  
end

% 3. Magnitud de la fuerza de arrastre
D = 0.5 * rho * Cd * A_m * v_rel_norm^2;

% 4. Aceleración de arrastre en ECI
a_drag = D * drag_dir;
end
