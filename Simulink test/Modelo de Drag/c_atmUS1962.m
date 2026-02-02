
% -------------------------------------------------------------------------
%                           c_atmUS1962
% -------------------------------------------------------------------------

% DESCRIPCIÓN:
% Función que calcula los parámetros atmosféricos para las capas de la
% atmósfera superiores a los 86 km hasta los 700 km geométricos 
% según el estándar US-1962.

% ENTRADAS:
% h: altitud geométrica      [m]
% dataAtm: estructura de datos de los modelos atmosféricos

% SALIDAS:
% Atm: estructura que contiene:
%           Atm.T: temperatura             [K]
%           Atm.P: presión                 [Pa]
%           Atm.rho: densidad              [kg/m3]
%           Atm.c: velocidad del sonido    [m/s]
%           Atm.L: mean free path          [m]
%           Atm.flag: indica si el cálculo se encuentra fuera del rango de
%           cálculo de la norma            [boolean]

% LIMITACIONES:
% La función es válida para altitudes geométricas mayores a 86 km. A
% partir de los 700 km considera que no hay más atmósfera. Presenta
% diferencias significativas para altitudes mayores a 86 km respecto la 
% norma US 1976 (posterior), sin embargo, el modelo es ampliamente 
% utilizado en la bibliografía ya que simplifica las leyes de variación.
% Se justifica el uso de este modelo considerando que las cargas aerodinámicas 
% son muy bajas para altitudes mayores a 120 km, y que para un modelado más
% realista se deberían incluir otros efectos atmosféricos que no tiene
% sentido contemplar para este tipo de simulador.

function rho = c_atmUS1962(h,dataAtm) 

% Altitud geométrica
Z = h-6378000;	% [m]

% Si la altitud es menor o igual al límite superior de la última capa 
% atmosférica se calcula
if Z <= dataAtm.US62.Zb(end)
    
% Determinación de factor beta, Tb y Pb para cada capa
% Valor inicial si H no está en una capa
idx = 0;

% Busca la capa para el Z
for i=1:length(dataAtm.US62.Zb)-1  
    % Para cada capa
    if (Z>=dataAtm.US62.Zb(i)) && (Z<dataAtm.US62.Zb(i+1))
        idx = i;
    end
    % Último valor de la última capa es <=
    if (Z>=dataAtm.US62.Zb(end-1)) && (Z<=dataAtm.US62.Zb(end))
        idx = i;
    end
end

% Asigna valores
if idx ~= 0
    Lb = dataAtm.US62.Lb(idx);	%[K/m]
    Tb = dataAtm.US62.Tb(idx);	%[K]
    Pb = dataAtm.US62.Pb(idx);	%[Pa]
    Zb = dataAtm.US62.Zb(idx);	%[m]    
    flag = 0;
else 
    Lb = NaN;	%[K/m]
    Tb = NaN;	%[K]
    Pb = NaN;	%[Pa]
    Zb = NaN;	%[m] 
    flag = 1;
end

% Determinación de M en función de la altura
Mint = interp1(dataAtm.US62.Zb,dataAtm.US62.M,Z);
 
% Determinación de temperatura
TM = Tb+Lb*(Z-Zb);                          % Temperatura molecular [K]
T = (Tb+Lb*(Z-Zb))*(Mint/dataAtm.US62.Mo);  % Temperatura cinética  [K]

% Determinación de presión
C1 = Lb/Tb;
C2 = -((dataAtm.US62.go/(dataAtm.US62.R*Lb))*(1+dataAtm.US62.B*((Tb/Lb)-Zb)));
C3 = dataAtm.US62.go*dataAtm.US62.B/(dataAtm.US62.R*Lb);

P = Pb*((C1*(Z-Zb)+1)^C2)*exp(C3*((Z-Zb)));	%[Pa]

%Determinación de densidad
rho = P*dataAtm.US62.Mo/(dataAtm.US62.Res*TM);	%[kg/m3]

%Determinación de la velocidad del sonido
c = (dataAtm.US62.gamma*dataAtm.US62.Res*TM/dataAtm.US62.Mo)^(1/2); %[m/s]

% Determinación de mean free path
L = (dataAtm.US62.Res*TM*Mint)/((2^0.5)*pi*dataAtm.US62.Na*(dataAtm.US62.sigma^2)*P*dataAtm.US62.Mo);

% Si la altitud es mayor al límite superior de la última capa se considera
% que no hay atmósfera
else
    
T = NaN;    % No definida
P = 0;      % [Pa]
rho = 0;    % [kg/m3]
c = NaN;    % No definida
L = NaN;    % No definida
flag = 0;
    
end

%Salidas
Atm.T = T ;
Atm.P = P ;
Atm.rho = rho ;
Atm.c = c ;
Atm.L = L ;
Atm.flag = flag ;

end