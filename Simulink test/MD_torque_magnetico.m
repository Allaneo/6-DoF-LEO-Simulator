function tau_mag = MD_torque_magnetico(B_body, m_res, m_cmd)
% MD_torque_magnetico
% Torque magnético en BODY: tau = (m_res + m_cmd) x B
%
% Entradas:
%   B_body_T    : [3x1] campo geomagnético en BODY [Tesla]
%   m_res_body  : [3x1] dipolo residual en BODY [A*m^2]
%   m_cmd_body  : [3x1] dipolo comandado en BODY [A*m^2]
%
% Salidas:
%   tau_mag : [3x1] torque magnético en BODY [N*m]

B_body = reshape(B_body, 3, 1);
m_res  = reshape(m_res,  3, 1);
m_cmd  = reshape(m_cmd,  3, 1);

m_tot  = m_res + m_cmd;
tau_mag = cross(m_tot, B_body);
end
