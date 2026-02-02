function tau_gg = MD_momento_gravitatorio(Ja_body,I)
% Calcula el torque por gradiente gravitatorio (gravity-gradient) en la 
% terna Body a partir del tensor de gradiente gravitatorio expresado en Body
% y de la matriz de inercia del satélite.

% Inputs:
%  Ja_body : 3x3, matriz de gradiente gravitatorio
%  I_body : 3x3, matriz de inercia en body

% Output:
%  tau_gg : 3x1, torque gravitacional por gradiente en body


% extraer componentes de Ja (coherencia de indices g_ij = G(i,j))
g11 = Ja_body(1,1); g12 = Ja_body(1,2); g13 = Ja_body(1,3);
g21 = Ja_body(2,1); g22 = Ja_body(2,2); g23 = Ja_body(2,3);
g31 = Ja_body(3,1); g32 = Ja_body(3,2); g33 = Ja_body(3,3);

% extraer componentes de I
Ixx = I(1,1);
Iyy = I(2,2);
Izz = I(3,3);
Ixy = I(1,2); % = I_body(2,1)
Ixz = I(1,3); % = I_body(3,1)
Iyz = I(2,3); % = I_body(3,2)

% ecuación 8-30 
tau1 = g23*(Izz - Iyy) + g13*Ixy - g12*Ixz + Iyz*(g33 - g22);
tau2 = g13*(Ixx - Izz) - g23*Ixy + g12*Iyz + Ixz*(g11 - g33);
tau3 = g12*(Iyy - Ixx) + g23*Ixz - g13*Iyz + Ixy*(g22 - g11);

tau_gg = [tau1; tau2; tau3];
end
