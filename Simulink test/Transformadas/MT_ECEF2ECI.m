function r_eci = MT_ECEF2ECI(r_ecef, era, C_FIX_val)
%   Convierte posición de ECEF a ECO
%   r_eci: vector de posición en ECI (3x1)
%   r_ecef: vector de posicion en ECEF (3x1)
%   era: earth rotation angle [0,2pi)

    c = cos(era); s = sin(era);
    R = [c, -s, 0;
          s, c, 0; 
          0, 0, 1];

    C = C_FIX_val'*R;
    r_eci = C * r_ecef;
end