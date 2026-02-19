function w_ref_eci = MT_LVLH_OMEGA_ECI(r_eci, v_eci)
% w_ref_eci = omega(LVLH/I) expressed in ECI
r = r_eci(:);
v = v_eci(:);
w_ref_eci = cross(r,v) / (dot(r,r));
end
