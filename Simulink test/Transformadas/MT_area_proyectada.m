function [Aproj, w, vhat] = MT_area_proyectada(v_rel_I, q, A_transversal, SAT_SHAPE)
% MT_area_proyectada
% SAT_SHAPE:
%   1 = Cannonball: Fuerza con A constante = A_transversal(1)
%       pero mantiene w para CP (torque aero activo).
%   2 = Prisma:     Fuerza con A proyectada = sum(A_transversal .* |vhat_B|)
%
% A_transversal SIEMPRE 3x1.

% v_rel en BODY
v_rel_B = MT_ECI2BODY(q, v_rel_I);
v = norm(v_rel_B);
vhat = v_rel_B / v;

% Pesos para CP (SIEMPRE, para mantener torque aero activo)
w = A_transversal(:) .* abs(vhat(:));

% modo
shp = SAT_SHAPE;
if isa(shp,'Simulink.Parameter'), shp = shp.Value; end
shp = double(shp);

if shp == 1
    % Cannonball: A constante (pero w se mantiene para CP)
    Aproj = A_transversal(1);

elseif shp == 2
    % Prisma: A proyectada
    Aproj = w(1) + w(2) + w(3);

else
    error('MT_area_proyectada: SAT_SHAPE=%g inválido. Usar 1=Cannonball, 2=Prisma.', shp);
end
end
