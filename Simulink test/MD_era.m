function [era,M] = MD_era(t)
% ERA_FROM_JD_UT1  Calcula Earth Rotation Angle (rad) a partir de JD(UT1)
% Basado en IERS/IAU: ERA = 2*pi*(0.7790572732640 + 1.00273781191135448*(JD_UT1 - 2451545.0))
% Devuelve valor reducido en [0,2*pi).

    jd_ut_inicial = 2451545.0;
    jd_ut1 =  jd_ut_inicial+ t/86400;
    fraction = jd_ut1 - 2451545.0;
    % compute ERA in revolutions
    rev = 0.7790572732640 + 1.00273781191135448 * fraction;
    % keep fractional part
    rev_frac = rev - floor(rev);
    era = 2*pi*rev_frac;

    T = (jd_ut1-jd_ut_inicial)/36525; 
    M = 357.5256+35999.049*T;
end
