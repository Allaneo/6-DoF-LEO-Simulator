function tau_mag = MD_torque_magnetico(B_body, m_res, m_cmd)
%#codegen
% Torque magnético en BODY: tau = (m_res + m_cmd) x B
% Entradas aceptadas:
%   - 3 elementos (cualquier forma) o escalar (se replica).
% Salida:
%   - tau_mag [3x1] fijo (codegen-friendly)

B  = to3x1(B_body);
mr = to3x1(m_res);
mc = to3x1(m_cmd);

tau_mag = cross(mr + mc, B);
end

function v = to3x1(x)
%#codegen
v = zeros(3,1);

% Nota: usar numel() + indexado lineal mantiene salida de tamaño fijo.
if numel(x) == 3
    v(1) = x(1);
    v(2) = x(2);
    v(3) = x(3);
elseif numel(x) == 1
    v(:) = x(1);
else
    % Fallo duro: si llega algo distinto, el modelo está mal cableado
    coder.internal.error('MD_torque_magnetico: input debe tener 1 o 3 elementos.');
end
end
