function q_norm = MD_normalize_q(q_raw)
% Normaliza el cuaternión de actitud para evitar deriva numérica durante la
% integración, manteninedo la convención scalar first y dejando la
% velocidad angular sin cambios.

% Entradas:
% q_raw: vector de cuaterniones sin normalizar

% Salida:
% q_norm: vector de cuaterniones normalizado.

q = q_raw(1:4);
w = q_raw(5:7);

n = sqrt(q.'*q);
q = q./n;


q_norm = [q; w];
end
