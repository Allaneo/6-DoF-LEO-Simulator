function r_CP = MD_cp_calculator(w, vhat, L,r_CG_in_G)
% --- MD_cp_calculator ---
% Calcula un centro de presión equivalente para un satélite modelado como paralelepípedo, 
% suponiendo que la fuerza aerodinámica resultante puede representarse como 
% una fuerza única aplicada en un punto equivalente. El punto se obtiene como
%  un promedio ponderado por contribuciones de área/proyección sobre los ejes
%  principales, seleccionando la cara "expuesta" segúnd la dirección del flujo
%  incidente (vhat). El resultado se entrega referido al centro de masas,
% incorporando el desplazamiento del centro de masas respecto del centro
% geométrico.

% Entradas:
% r_CG_in_G: desplazamiento del centro de masas respecto del centro geométrico
% L = [Lx; Ly; Lz] en metros
% w: pesos relativos de cada área
% vhat: vector unitario que indica la dirección de la velocidad relativa

% Salida:
% r_CP: posición del centro de presiones en la terna body.

W = w(1) + w(2) + w(3);

sx = -sign(vhat(1));
sy = -sign(vhat(2));
sz = -sign(vhat(3));

r_CP_G = [ (sx*L(1)/2)*w(1);
           (sy*L(2)/2)*w(2);
           (sz*L(3)/2)*w(3) ] / W;
r_CP   = r_CP_G - r_CG_in_G;
end
