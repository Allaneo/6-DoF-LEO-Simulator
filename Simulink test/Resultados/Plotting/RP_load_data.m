function D = RP_load_data(logsout_in, SIG, sigKeys, sigNames, H)
% Carga señales desde logsout (Dataset) o desde SDI (último run).
% Devuelve estructura D con:
%  - D.data.(key).ok/name/t/X
%  - D.tx, D.rECI, D.vECI
%  - D.haveAtt, D.tA, D.q_raw, D.wB
%  - D.src_mode

D = struct();
D.data = struct();
D.src_mode = '';

logs = [];

sdiNames = {};
sdiTS    = {};
sdiRunID = [];

%% 1) logsout
if ~isempty(logsout_in)
    try
        if isa(logsout_in,'Simulink.SimulationData.Dataset') && logsout_in.numElements > 0
            logs = logsout_in;
            D.src_mode = 'workspace.logsout (Dataset)';
        end
    catch
    end
end

%% 2) SDI
if isempty(D.src_mode)
    runIDs = Simulink.sdi.getAllRunIDs;
    if isempty(runIDs)
        error('No encontré logsout y tampoco hay runs en SDI.');
    end

    sdiRunID = runIDs(end);
    r = Simulink.sdi.getRun(sdiRunID);

    ids = r.getSignalIDs;
    nS = numel(ids);

    sdiNames = cell(nS,1);
    sdiTS    = cell(nS,1);

    for i = 1:nS
        sigObj = Simulink.sdi.getSignal(ids(i));
        sdiNames{i} = sigObj.Name;
        sdiTS{i} = H.get_sdi_timeseries(sigObj);
    end

    D.src_mode = sprintf('SDI (último run ID %d)', sdiRunID);
end

%% 3) señales genéricas
for k = 1:numel(sigKeys)
    key  = sigKeys{k};
    name = sigNames{k};

    D.data.(key).ok   = false;
    D.data.(key).name = name;
    D.data.(key).t    = [];
    D.data.(key).X    = [];

    ts = [];
    if ~isempty(logs)
        ts = H.get_dataset_timeseries_exact(logs, name);
    else
        ts = H.get_sdi_timeseries_by_name(sdiNames, sdiTS, name);
    end

    if ~isempty(ts)
        [t, X] = H.read_ts(ts);
        D.data.(key).ok = true;
        D.data.(key).t  = t;
        D.data.(key).X  = X;
    end
end

%% 4) x_lineal obligatorio
xTS = [];
if ~isempty(logs)
    xTS = H.get_dataset_timeseries_exact(logs, SIG.x_lineal);
else
    xTS = H.get_sdi_timeseries_by_name(sdiNames, sdiTS, SIG.x_lineal);
end
if isempty(xTS)
    error('No encontré la señal "%s".', SIG.x_lineal);
end

[tx, Xstate] = H.read_ts(xTS);
if size(Xstate,2) < 6
    error('%s debe tener al menos 6 componentes (r(1:3), v(4:6)).', SIG.x_lineal);
end

D.tx   = tx(:);
D.rECI = Xstate(:,1:3);
D.vECI = Xstate(:,4:6);
D.xNameUsed = SIG.x_lineal;

%% 5) x_angular opcional
xAngTS = [];
if ~isempty(logs)
    xAngTS = H.get_dataset_timeseries_exact(logs, SIG.x_angular);
else
    xAngTS = H.get_sdi_timeseries_by_name(sdiNames, sdiTS, SIG.x_angular);
end

D.haveAtt = false;
D.tA = [];
D.q_raw = [];
D.wB = [];

if ~isempty(xAngTS)
    [tA, Xang] = H.read_ts(xAngTS);
    if size(Xang,2) >= 7
        D.haveAtt = true;
        D.tA = tA(:);
        D.q_raw = Xang(:,1:4);
        D.wB = Xang(:,5:7);
    end
end
end
