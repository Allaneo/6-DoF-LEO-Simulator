function v_rel_I = MD_velocidad_relativa(v_I,r_I)
% constantes
omega_earth = [0;0;7.2921150e-5];

% 1) velocidad relativa en ECI (atmósfera co-rotante)
v_rel_I = v_I - cross(omega_earth, r_I);
end