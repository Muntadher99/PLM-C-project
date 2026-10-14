%IDEAL_SWITCH_VALIDATION  Ideal-switch Simulink check of the analytic
%   harmonic model (paper, Section 6.2, first sentences).
%
%   Simulates every angle set of the input workbook in NPC_SHM_ideal.slx
%   (ideal switches), computes THD49 of the line-to-line voltage over exactly
%   the 16 counted orders (odd non-triplen 5..49) and compares it with the
%   analytic THD49 of the same angles. The full 2..49 sum is kept as a
%   separate diagnostic (out-of-band content).
%
%   Input : validation_cases.xlsx, sheet 1 (export_validation_cases)
%   Output: validation_cases_sim.xlsx, one row per case (analytic and
%           simulated THD49, absolute and relative error, WTHD49, current
%           THD, neutral-point voltage statistics, status)
%
%   Requires Simulink, Simscape Electrical and the Parallel Computing Toolbox
%   (parsim); MATLAB R2024b.

clc; clear; close all;

%% ------------------------------- settings --------------------------------
excelFile   = 'validation_cases.xlsx';        % input (angles + analytic THD49)
outputExcel = 'validation_cases_sim.xlsx';    % output

maCol            = 'Ma';
thdCol           = 'THD49_pct';
armCol           = 'arm';             % '' if the sheet holds a single method
angleFormat      = 'theta%02d_deg';
nAngles          = 15;
anglesAreDegrees = true;

modelName          = 'NPC_SHM_ideal';
alphaConstantBlock = 'NPC_SHM_ideal/Constant';   % block holding the angle vector (rad)

voltageSignalName = 'vab';
currentSignalName = 'ia';
dVnpSignalName    = 'dVnp';

solverName = 'ode23tb';
maxStep    = '1e-6';

f1              = 50;
maxH            = 49;
simTime         = 1.0;
analysisCycles  = 10;      % an integer number of cycles -> exact FFT bins
resampleStep    = 1e-6;

testOnly    = false;
nTestCases  = 1;

% counted band: odd, non-triplen, 5 <= n <= 49 (P.ORD in plm_core)
ORDERS = 5:2:maxH;
ORDERS = ORDERS(mod(ORDERS,3) ~= 0);

%% ------------------------------- input -----------------------------------
analysisEnd   = simTime;
analysisStart = simTime - analysisCycles/f1;

fprintf('Reading input file: %s\n', excelFile);
data = readtable(excelFile);
if testOnly, data = data(1:min(nTestCases, height(data)), :); end

angleCols    = arrayfun(@(k) sprintf(angleFormat,k), 1:nAngles, 'UniformOutput', false);
requiredCols = [{maCol, thdCol}, angleCols];
missing      = setdiff(requiredCols, data.Properties.VariableNames);
if ~isempty(missing)
    error('Missing columns in %s: %s', excelFile, strjoin(missing, ', '));
end

Ma_list           = data.(maCol);
THD49_from_solver = data.(thdCol);
alpha_all         = data{:, angleCols};
nCases            = height(data);

hasArm = ~isempty(armCol) && ismember(armCol, data.Properties.VariableNames);
if hasArm, armName = string(data.(armCol)); else, armName = repmat("(single)", nCases, 1); end

fprintf('Loaded %d cases.  Angles: %s\n', nCases, ...
        ternary(anglesAreDegrees,'degrees','radians'));
fprintf('Analysis window: %.6f to %.6f s (%d whole cycles at %g Hz)\n', ...
        analysisStart, analysisEnd, analysisCycles, f1);
fprintf('Counted band: %d orders -> %s\n', numel(ORDERS), mat2str(ORDERS));
if hasArm
    u = unique(armName);
    fprintf('Methods present: %s\n', strjoin(cellstr(u).', ', '));
end

%% ------------------------------- simulate --------------------------------
load_system(modelName);
simIn(nCases,1) = Simulink.SimulationInput(modelName);
for k = 1:nCases
    alpha = alpha_all(k,:);
    if anglesAreDegrees, alpha = deg2rad(alpha); end
    simIn(k) = Simulink.SimulationInput(modelName);
    simIn(k) = simIn(k).setVariable('Ma', Ma_list(k));
    simIn(k) = simIn(k).setBlockParameter(alphaConstantBlock,'Value',mat2str(alpha,16));
    simIn(k) = simIn(k).setModelParameter('StopTime',num2str(simTime), ...
                 'SolverType','Variable-step','Solver',solverName,'MaxStep',maxStep);
end

if isempty(gcp('nocreate')), parpool; end
fprintf('\nStarting parallel simulations (%d cases)...\n', nCases);
tic;
simOut = parsim(simIn,'ShowProgress','on','TransferBaseWorkspaceVariables','on');
fprintf('Finished in %.2f minutes.\n', toc/60);

%% ------------------------------- post-process ----------------------------
Vab_THD_band = nan(nCases,1);   Vab_WTHD_band = nan(nCases,1);
Vab_THD_full = nan(nCases,1);   OOB_pct       = nan(nCases,1);
Ia_THD_band  = nan(nCases,1);   Ia_WTHD_band  = nan(nCases,1);
V1_rms = nan(nCases,1);         I1_rms = nan(nCases,1);
dVnp_max_abs = nan(nCases,1);   dVnp_rms_val = nan(nCases,1);
dVnp_mean_val= nan(nCases,1);
simStatus = strings(nCases,1);  errorMessage = strings(nCases,1);

for k = 1:nCases
    fprintf('Processing case %d/%d  (%s, Ma = %.4f)\n', k, nCases, armName(k), Ma_list(k));
    try
        if ~isempty(simOut(k).ErrorMessage)
            simStatus(k) = "SimFailed"; errorMessage(k) = string(simOut(k).ErrorMessage);
            continue;
        end
        simStatus(k) = "Success";

        [t_v , v_ab ] = getTimeAndData(simOut(k).get(voltageSignalName));
        [t_i , i_a  ] = getTimeAndData(simOut(k).get(currentSignalName));
        [t_np, dV_np] = getTimeAndData(simOut(k).get(dVnpSignalName));

        [Vab_THD_band(k), Vab_WTHD_band(k), V1_rms(k), Vab_THD_full(k)] = ...
            harmonics(v_ab, t_v, f1, maxH, ORDERS, analysisStart, analysisEnd, resampleStep);
        [Ia_THD_band(k),  Ia_WTHD_band(k),  I1_rms(k), ~] = ...
            harmonics(i_a , t_i, f1, maxH, ORDERS, analysisStart, analysisEnd, resampleStep);

        % out-of-band content: everything in 2..49 outside the counted band
        OOB_pct(k) = sqrt(max(Vab_THD_full(k)^2 - Vab_THD_band(k)^2, 0));

        idx = (t_np >= analysisStart) & (t_np < analysisEnd);
        d   = dV_np(idx);
        if numel(d) < 10, error('Not enough dVnp samples in the window.'); end
        dVnp_max_abs(k) = max(abs(d));
        dVnp_rms_val(k) = rms(d);
        dVnp_mean_val(k)= mean(d);
    catch ME
        simStatus(k) = "PostProcessFailed"; errorMessage(k) = string(ME.message);
    end
end

%% ------------------------------- results ---------------------------------
absErr = abs(Vab_THD_band - THD49_from_solver);
relErr = 100 * absErr ./ max(abs(THD49_from_solver), eps);

Results = table();
Results.arm            = armName;
Results.Ma             = Ma_list;
Results.THD49_analytic = THD49_from_solver;
Results.THD49_simulink = Vab_THD_band;
Results.abs_error_pp   = absErr;
Results.rel_error_pct  = relErr;
Results.THD_full_2to49 = Vab_THD_full;
Results.out_of_band_pct= OOB_pct;
Results.WTHD49_simulink= Vab_WTHD_band;
Results.Ia_THD49       = Ia_THD_band;
Results.Ia_WTHD49      = Ia_WTHD_band;
Results.V1_rms         = V1_rms;
Results.I1_rms         = I1_rms;
Results.dVnp_max_abs_V = dVnp_max_abs;
Results.dVnp_rms_V     = dVnp_rms_val;
Results.dVnp_mean_V    = dVnp_mean_val;

unitSuffix = ternary(anglesAreDegrees,'_deg','_rad');
for n = 1:nAngles
    Results.(sprintf('theta_%02d%s', n, unitSuffix)) = alpha_all(:,n);
end
Results.Status       = simStatus;
Results.ErrorMessage = errorMessage;

writetable(Results, outputExcel);
fprintf('\nResults saved to: %s\n', outputExcel);

ok = simStatus == "Success";
fprintf('\n%s\n', repmat('=',1,78));
fprintf('ANALYTIC MODEL vs IDEAL-SWITCH SIMULATION (THD49 over the 16 counted orders)\n');
fprintf('%s\n', repmat('=',1,78));
if ~any(ok)
    fprintf('  no successful cases\n');
else
    u = unique(armName(ok));
    fprintf('%-16s %6s %14s %14s %14s\n','method','n','max |err| pp','mean |err| pp','max OOB %');
    for i = 1:numel(u)
        m = ok & (armName == u(i));
        fprintf('%-16s %6d %14.5f %14.5f %14.5f\n', u(i), sum(m), ...
                max(absErr(m)), mean(absErr(m)), max(OOB_pct(m)));
    end
    fprintf('%-16s %6d %14.5f %14.5f %14.5f\n','ALL', sum(ok), ...
            max(absErr(ok)), mean(absErr(ok)), max(OOB_pct(ok)));
end
fprintf('%s\n\n', repmat('=',1,78));

%% ------------------------------- helpers ---------------------------------
function out = ternary(c,a,b), if c, out=a; else, out=b; end, end

function [t,x] = getTimeAndData(sig)
    if isa(sig,'Simulink.SimulationData.Dataset'), sig = sig{1}.Values; end
    if isa(sig,'timeseries')
        t = sig.Time; x = sig.Data;
    elseif isstruct(sig) && isfield(sig,'time')
        t = sig.time; x = sig.signals.values;
    else
        error('Unsupported signal format. Use "Timeseries" or signal logging.');
    end
    t = t(:); if size(x,2) > 1, x = x(:,1); end
    x = x(:);
end

function [THD_band, WTHD_band, fund_rms, THD_full] = ...
         harmonics(x, t, f1, maxH, ORDERS, tStart, tEnd, Ts)
%HARMONICS  THD and WTHD over an explicit order list, plus the 2..maxH sum.
    x = x(:); t = t(:);
    v = isfinite(x) & isfinite(t);  x = x(v); t = t(v);
    [t,i] = sort(t); x = x(i);
    [t,i] = unique(t); x = x(i);
    if numel(x) < 100,                error('Too few samples.'); end
    if t(1) > tStart || t(end) < tEnd, error('Window not covered by the log.'); end

    Nq = round((tEnd - tStart)/Ts);
    tq = tStart + (0:Nq-1).'*Ts;          % integer cycles -> exact FFT bins
    if numel(tq) < 1000, error('Uniform grid too short.'); end
    xq = interp1(t, x, tq, 'linear');
    if any(~isfinite(xq)), error('Interpolation produced NaN/Inf.'); end
    xq = xq - mean(xq);

    N = numel(xq); Fs = 1/Ts; X = fft(xq);
    hr = nan(maxH,1);
    for h = 1:maxH
        fh = h*f1;
        if fh >= Fs/2, continue; end
        bin = round(fh*N/Fs) + 1;
        if bin > numel(X), continue; end
        hr(h) = (2*abs(X(bin))/N)/sqrt(2);
    end

    fund_rms = hr(1);
    if ~isfinite(fund_rms) || fund_rms <= 0
        THD_band = NaN; WTHD_band = NaN; THD_full = NaN; return;
    end

    ob = ORDERS(isfinite(hr(ORDERS)));                       % the counted band
    THD_band  = sqrt(sum(hr(ob).^2))            / fund_rms * 100;
    WTHD_band = sqrt(sum((hr(ob)./ob(:)).^2))   / fund_rms * 100;

    af = 2:maxH;  af = af(isfinite(hr(af)));                 % diagnostic only
    THD_full  = sqrt(sum(hr(af).^2))            / fund_rms * 100;
end
