function [era,M] = MD_era(t,jd_ut_inicial)
% Calcula el ángulo de rotación de la Tierra (ERA) a partir de la fecha 
% inicial en tiempo juliano y el tiempo transcurrido de simulación, de modo
% que pueda usarse como entrada para transformar entre marcos inercial y 
% terrestre. Además, calcula la anomalía media del Sol asociada al instante
% actual, para utilizarla en el modelo simplificado de posición solar.

    jd_ut1 =  jd_ut_inicial+ t/86400;
    fraction = jd_ut1 - 2451545.0;

    % compute ERA in revolutions
    rev = 0.7790572732640 + 1.00273781191135448 * fraction;
    % keep fractional part
    rev_frac = rev - floor(rev);
    era = 2*pi*rev_frac;

    T = (jd_ut1 - jd_ut_inicial)/36525; %cantidad de siglos
    M = 357.5256+35999.049*T;
end
