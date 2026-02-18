function RP_plot_attitude(D, cfg, H)
if ~isfield(D.derived,'att') || ~D.derived.att.ok
    return;
end

shadowSegments = D.derived.shadowSegments;
A = D.derived.att;

figure('Name','Actitud - Euler ZYX (deg)','Color','w');
axE = axes; hold(axE,'on'); grid(axE,'on');
plot(axE, A.tA, A.yaw_deg,   'LineWidth', 1.3);
plot(axE, A.tA, A.pitch_deg, 'LineWidth', 1.3);
plot(axE, A.tA, A.roll_deg,  'LineWidth', 1.3);
xlabel(axE,'Tiempo [s]'); ylabel(axE,'Ángulo [deg]');
title(axE,'Euler ZYX desde q (C_{IB} fijo) - osculante');
H.apply_shadow(axE, shadowSegments);
legend(axE, {'Yaw \psi','Pitch \theta','Roll \phi'}, 'Location','best');

figure('Name','Actitud - Rates (rad/s)','Color','w');
axW = axes; hold(axW,'on'); grid(axW,'on');
plot(axW, A.tA, A.wB(:,1), 'LineWidth', 1.2);
plot(axW, A.tA, A.wB(:,2), 'LineWidth', 1.2);
plot(axW, A.tA, A.wB(:,3), 'LineWidth', 1.2);
plot(axW, A.tA, A.omega_mag, 'LineWidth', 1.6);
xlabel(axW,'Tiempo [s]'); ylabel(axW,'\omega [rad/s]');
title(axW,'Rates (body) + norma');
H.apply_shadow(axW, shadowSegments);
legend(axW, {'\omega_x','\omega_y','\omega_z','|\omega|'}, 'Location','best');

figure('Name','Actitud - Error a nadir (deg)','Color','w');
axP = axes; hold(axP,'on'); grid(axP,'on');
plot(axP, A.tA, A.nadir_err, 'LineWidth', 1.6);
xlabel(axP,'Tiempo [s]'); ylabel(axP,'Error [deg]');
title(axP, sprintf('Error apuntamiento a nadir (axis=[%.0f %.0f %.0f]_B, C_{IB} fijo)', ...
    cfg.body_axis_for_nadir(1), cfg.body_axis_for_nadir(2), cfg.body_axis_for_nadir(3)));
H.apply_shadow(axP, shadowSegments);
end
