function RP_print_summary(D)
fprintf('\n=== RESUMEN DISPONIBILIDAD DE SEÑALES ===\n');

keys = fieldnames(D.data);
for i = 1:numel(keys)
    k = keys{i};
    if D.data.(k).ok
        fprintf('[OK]   %s\n', D.data.(k).name);
    else
        fprintf('[MISS] %s\n', D.data.(k).name);
    end
end

fprintf('[OK]   Estado usado (órbita):  %s\n', D.xNameUsed);

if D.haveAtt
    fprintf('[OK]   Estado usado (actitud): %s\n', 'x_angular');
    if isfield(D.derived,'att') && D.derived.att.ok
        fprintf('[INFO] Max |omega| = %.6g rad/s\n', max(D.derived.att.omega_mag));
        fprintf('[INFO] Max (||q_raw||-1) = %.6g\n', max(abs(D.derived.att.q_norm_err)));
        fprintf('[INFO] Max nadir err = %.6g deg\n', max(D.derived.att.nadir_err));
    end
else
    fprintf('[MISS] Estado de actitud: %s\n', 'x_angular');
end
end
