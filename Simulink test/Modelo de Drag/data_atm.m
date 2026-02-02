
% -------------------------------------------------------------------------
%                                data_atm
% -------------------------------------------------------------------------

% DESCRIPCIÓN:
% Función que establece los parámetros de atmosféricos requeridos para 
% ejecutar lo modelos de atmósfera. Tambien establece que modelo sera
% utilizado para bajas altitutes

% ENTRADAS
% modelo: Modelo atmosferico a utilizar para bajas altitudes.
% ("ISO2533"/"US1976")
% 
% SALIDAS:
% dataAtm: estructura de datos para cada norma que contiene 
%           - dataAtm.ISO: datos del modelo ISO-2533
%           - dataAtm.US76: datos del modelo US-1976
%           - dataAtm.US62: datos del modelo US-1962
%           - dataAtm.modelo: booleano (1 = ISO2533) / (0=US1976)
% Dentro contiene:
% - Po=Pn: Presión a nivel del mar [Pa];
% - To=Tn: Temperatura a nivel del mar [K];
% - go=gn: Constante gravitatoria a nivel del mar [m/s2];
% - Mo: Peso molecular del aire a nivel del mar [adim];
% - R: Constante molar específica del aire [J/kg.K];
% - Res: Constante molar específica del aire [J/kg.mol.K];
% - gamma=k: Relación de Cp y Cv del aire [adim];
% - B: Inversa de radio terrestre esférico [1/m];
% - Na: Número de Avogadro [1/kg.mol];
% - sigma: Diámetro de colisión del aire [m];
% - Hb: Vector de altitudes geopotenciales de las capas [m];
% - Zb: Vector de altitudes geométricas de las capas [m];
% - Tb: Vector de temperaturas moleculares de las capas [K];
% - Lb: Vector de pendiente de TM de las capas [K/m];
% - Pb: Vector de presiones de las capas [Pa];
% - M: Vector de peso molecular del aire para Zb [adim];
% - MMo: Vector de relación M/Mo [adim];
% - HMMo: Vector de altitudes  de MMo [m];

function dataAtm = data_atm(modelo)


%% Determinacion del modelo a utilizar

switch modelo
    case 'ISO2533'
        dataAtm.modelo = 1;
    case 'US1976'
        dataAtm.modelo = 0;
end

%% Datos de los modelos

% ---------------
% Modelo ISO-2533
% ---------------

% Constantes modelo ISO-2533
Pn = 101325;        % [Pa]
Tn = 288.15;        % [K]
gn = 9.80665;       % [m/s2]
R = 287.05287;      % [J/kg.K]
Res = 8314.32;      % [J/kmol.K]
Na = 6.02257e26;    % [1/kmol]
sigma = 3.65e-10;   % [m]
k = 1.4;            % [adim]

% Datos de las capas
Ncapas = 8;
Hb = [-2 0 11 20 32 47 51 71 80]*1e3;
Tb = [301.15 Tn 216.65 216.65 228.65 270.65 270.65 214.65];
Lb = [-6.5 -6.5 0 1 2.8 0 -2.8 -2]*1e-3;
Pb = [127773.73015 Pn 0 0 0 0 0 0];

% Determinación de Pb de cada capa
for i=1:Ncapas-1
    
    if Lb(i) ~= 0
        Pb(i+1) = Pb(i)*(1+(Lb(i)/Tb(i))*(Hb(i+1)-Hb(i)))^(-gn/(Lb(i)*R));    %[Pa]
    else
        Pb(i+1) = Pb(i)*exp((-gn*(Hb(i+1)-Hb(i)))/(R*Tb(i)));                 %[Pa]
    end
    
end

% Guardo datos
dataAtm.ISO.Pn = Pn;
dataAtm.ISO.Tn = Tn;
dataAtm.ISO.gn = gn;
dataAtm.ISO.R = R;
dataAtm.ISO.Res = Res;
dataAtm.ISO.k = k;
dataAtm.ISO.Na = Na;
dataAtm.ISO.sigma = sigma;
dataAtm.ISO.Hb = Hb;
dataAtm.ISO.Tb = Tb;
dataAtm.ISO.Lb = Lb;
dataAtm.ISO.Pb = Pb;

clear Pn Tn gn R Res k Na sigma Ncapas Hb Tb Lb Pb i

% --------------
% Modelo US-1976
% --------------

% Constantes modelo US-1976
Po = 101325;        % [Pa]
To = 288.15;        % [K]
go = 9.80665;       % [m/s2]
Res = 8314.32;      % [J/kmol.K]
Mo = 28.9644;       % [kg/kmol]
gamma = 1.4;        % [adim]
Na = 6.022169e26;   % [1/kmol]
sigma = 3.65e-10;   % [m]

% Datos de las capas
Ncapas = 7;
Hb = [0 11 20 32 47 51 71 84.852]*1e3;
Tb = [To 216.65 216.65 228.65 270.65 270.65 214.65];
Lb = [-6.5 0 1 2.8 0 -2.8 -2]*1e-3;
Pb = [Po 0 0 0 0 0 0];
MMo   = [1 0.999996 0.999989 0.999971 0.999941 0.999909 0.999870 0.999829 0.999786 0.999741 0.999694 0.999641 0.999579];
HMMo = [79005.7 79493.3 79980.8 80468.2 80955.7 81443.0 81930.2 82417.3 82904.4 83391.4 83878.4 84365.2 84852.0]; 

% Determinación de Pb de cada capa
for i=1:Ncapas-1
    
    if Lb(i) ~= 0
        Pb(i+1) = Pb(i)*(Tb(i)/(Tb(i)+Lb(i)*(Hb(i+1)-Hb(i))))^(go*Mo/(Res*Lb(i)));     %[Pa]
    else
        Pb(i+1) = Pb(i)*exp(-go*Mo*(Hb(i+1)-Hb(i))/(Res*Tb(i)));                       %[Pa]
    end
    
end

% Guardo datos
dataAtm.US76.Po = Po;
dataAtm.US76.To = To;
dataAtm.US76.go = go;
dataAtm.US76.Res = Res;
dataAtm.US76.Mo = Mo;
dataAtm.US76.gamma = gamma;
dataAtm.US76.Na = Na;
dataAtm.US76.sigma = sigma;
dataAtm.US76.Hb = Hb;
dataAtm.US76.Tb = Tb;
dataAtm.US76.Lb = Lb;
dataAtm.US76.Pb = Pb;
dataAtm.US76.MMo = MMo;
dataAtm.US76.HMMo = HMMo;

clear Po To go Res Mo gamma Na sigma Ncapas Hb Tb Lb Pb MMo HMMo i
   
% --------------
% Modelo US-1962
% --------------

% Constantes modelo US-1962
Po = 101325 ;           % [Pa]
To = 288.15;            % [K]
go = 9.80665;           % [m/s2]
Res = 8314.32;          % [J/kmol.K]
Mo = 28.9644;           % [kg/kmol]
gamma = 1.4;            % [adim]
Na = 6.02257e26;        % [1/kmol]
sigma = 3.65e-10;       % [m]
B = 2/(6.3781*1e6);     % [1/m]
R = Res/Mo;             % [J/kg.K]

% Datos de las capas
Ncapas = 13;
Zb = [85.99995 100 110 120 150 160 170 190 230 300 400 500 600 700]*1e3;
Tb = [186.946 210.65 260.65 360.65 960.65 1110.65 1210.65 1350.65 1550.65 1830.65 2160.65 2420.65 2590.65 2700.65];
Lb = [1.6481 5 10 20 15 10 7 5 4 3.3 2.6 1.7 1.1]*1e-3;
Pb = [0.373383589976216 0 0 0 0 0 0 0 0 0 0 0 0];
M = [Mo 28.88 28.56 28.07 26.92 26.66 26.40 25.85 24.70 22.66 19.94 17.94 16.84 16.17];

% Determinación de Pb de cada capa
for i=1:Ncapas-1
    
    C1 = Lb(i)/Tb(i);
    C2 = -((go/(R*Lb(i)))*(1+(B*((Tb(i)/Lb(i))-Zb(i)))));
    C3 = go*B/(R*Lb(i));

    Pb(i+1) = Pb(i)*((C1*(Zb(i+1)-Zb(i))+1)^C2)*exp(C3*((Zb(i+1)-Zb(i))));	 %[Pa]

end

% Guardo datos
dataAtm.US62.Po = Po;
dataAtm.US62.To = To;
dataAtm.US62.go = go;
dataAtm.US62.Res = Res;
dataAtm.US62.Mo = Mo;
dataAtm.US62.gamma = gamma;
dataAtm.US62.Na = Na;
dataAtm.US62.sigma = sigma;
dataAtm.US62.R = R;
dataAtm.US62.B = B;
dataAtm.US62.Zb = Zb;
dataAtm.US62.Tb = Tb;
dataAtm.US62.Lb = Lb;
dataAtm.US62.Pb = Pb;
dataAtm.US62.M = M;

clear Po To go Res Mo gamma Na sigma B R Ncapas Zb Tb Lb Pb M i

end
