function T = nonideal_simulink_runs(casesfile, CFG)
%NONIDEAL_SIMULINK_RUNS  Device-level Simulink validation (paper, Section
%   6.2 and Table 8; Reviewer 1 comment 4, Reviewer 2 comment 3).
%
%   T = nonideal_simulink_runs('nonideal_cases.csv')
%
%   casesfile : CSV with columns design, Ma, theta01_deg ... theta15_deg
%               (ten rows: the distortion-minimising design and the compliant
%               witness at Ma = 0.40 0.60 0.80 1.00 1.20; written by
%               build_witness_tables)
%
%   Every angle set is simulated in NPC_SHM.slx under each of six cases:
%     ideal switches | device drops + switching | + 2 us dead time |
%     + DC-link imbalance (2% offset, +-10% capacitance) | + 1 us timer |
%     all combined
%   The line-to-line voltage is analysed over the last CFG.cycles periods;
%   each case is also evaluated with nonideal_switching_model.
%
%   The model reads these base-workspace variables (set here): Ron, Rond,
%   Vf, Vfd, Tf, Tt, Rs, Cs, Td, C1, C2, VC1_0, VC2_0, theta_deg, Ma, and logs
%   the line voltage to the variable named in CFG.log_var.
%   nonideal_model_init sets the same variables for interactive use.
%
%   Output CSV (CFG.out, default nonideal_results.csv), one row per (case,
%   angle set): case, design, Ma, THD49_analytic, THD49_sim, dTHD49_pp,
%   maxq_sim, b1_err_pct_sim, even_max_pct, oob_pct, THD49_swfn, maxq_swfn,
%   b1_err_pct_swfn, run_s.
%
%   Requires Simulink and Simscape Electrical (MATLAB R2024b).

% ======================= model interface ==================================
D.model      = 'NPC_SHM';      % Simulink model (without extension)
D.angle_var  = 'theta_deg';    % 1 x 15 angle vector (degrees) read by the gate logic
D.ma_var     = 'Ma';
D.log_var    = 'vab';          % To Workspace block logging v_ab (Timeseries)
D.Vdc        = 800;            % V
D.f1         = 50;             % Hz
D.stop       = 0.30;           % s simulated
D.cycles     = 10;             % fundamental periods analysed (at the end)
D.Nfft       = 2^16;           % uniform samples per period for the FFT
D.out        = 'nonideal_results.csv';   % output CSV
% ==========================================================================
if nargin < 2, CFG = struct(); end
fn = fieldnames(CFG); for i = 1:numel(fn), D.(fn{i}) = CFG.(fn{i}); end
C = D;

% ---------------- the cases ----------------------------------------------
d0 = device_parameters();
dev.Ron = d0.rce;  dev.Vf = d0.Vce0;  dev.Vfd = d0.Vf0;   % IGBT slope, IGBT and diode knee
dev.Rond = d0.rf;                                         % diode slope
dev.Tf  = d0.Tf;   dev.Tt = d0.Tt;                        % current fall / tail time
ideal   = struct('Ron',1e-3,'Rond',1e-3,'Vf',0,'Vfd',0,'Tf',1e-9,'Tt',1e-9,'Td',0, ...
                 'epsdc',0,'Cmis',0,'tq',0);
K = {};
K{end+1} = setf(ideal, 'name','ideal switches');
K{end+1} = setf(ideal, 'name','device drops + switching', 'Ron',dev.Ron,'Rond',dev.Rond,'Vf',dev.Vf, ...
                'Vfd',dev.Vfd,'Tf',dev.Tf,'Tt',dev.Tt);
K{end+1} = setf(K{2}, 'name','+ dead time 2 us', 'Td',2e-6);
K{end+1} = setf(K{2}, 'name','+ DC-link imbalance', 'epsdc',0.02, 'Cmis',0.10);
K{end+1} = setf(K{2}, 'name','+ 1 us timer', 'tq',1e-6);
K{end+1} = setf(K{2}, 'name','all combined', 'Td',2e-6,'epsdc',0.02,'Cmis',0.10,'tq',1e-6);
% epsdc: static neutral-point offset (fraction of Vdc/2), set through the
%        initial capacitor voltages; Cmis: C1 = C(1-Cmis), C2 = C(1+Cmis)

X = read_cases(casesfile);
Cnom = 470e-6;  Rs = 1e5;  Cs = Inf;                % resistive snubber
rows = {};
for kc = 1:numel(K)
    k = K{kc};
    assignin('base','Ron',k.Ron); assignin('base','Rond',k.Rond);
    assignin('base','Vf',k.Vf);   assignin('base','Vfd',k.Vfd);
    assignin('base','Tf',k.Tf);   assignin('base','Tt',k.Tt); assignin('base','Td',k.Td);
    assignin('base','Rs',Rs);     assignin('base','Cs',Cs);
    assignin('base','C1',Cnom*(1-k.Cmis)); assignin('base','C2',Cnom*(1+k.Cmis));
    assignin('base','VC1_0',C.Vdc/2*(1+k.epsdc)); assignin('base','VC2_0',C.Vdc/2*(1-k.epsdc));
    for r = 1:numel(X.Ma)
        th = X.theta(r,:);
        if k.tq > 0, d = 360*C.f1*k.tq; th = round(th/d)*d; end
        assignin('base', C.angle_var, th);
        assignin('base', C.ma_var, X.Ma(r));
        t0 = tic;
        out = sim(C.model, 'StopTime', num2str(C.stop), 'ReturnWorkspaceOutputs', 'on');
        [t, v] = get_signal(out, C.log_var);
        Msim = spectrum_metrics(t, v, C);
        % the same case in the switching-function model
        NI = struct('drops', k.Vf > 0, 'td', k.Td, 'eps_dc', k.epsdc, ...
                    'Ctot', 2*Cnom, 'tq', k.tq);
        W = nonideal_switching_model(X.theta(r,:), X.Ma(r), NI, struct('Vdc',C.Vdc,'f1',C.f1));
        rows(end+1,:) = {k.name, X.design{r}, X.Ma(r), X.THD49(r), Msim.THD49, ...
                         Msim.THD49 - X.THD49(r), Msim.maxq, Msim.b1_err_pct, ...
                         Msim.even_max, Msim.oob, W.THD49, W.maxq, W.b1_err_pct, toc(t0)}; %#ok<AGROW>
        fprintf('%-26s %-10s Ma=%.2f  THD49 %.4f%% (analytic %.4f, d %+.4f pp)  maxq %.3f  b1err %+.2f%%\n', ...
            k.name, X.design{r}, X.Ma(r), Msim.THD49, X.THD49(r), Msim.THD49 - X.THD49(r), ...
            Msim.maxq, Msim.b1_err_pct);
    end
end
hdr = {'case','design','Ma','THD49_analytic','THD49_sim','dTHD49_pp','maxq_sim', ...
       'b1_err_pct_sim','even_max_pct','oob_pct','THD49_swfn','maxq_swfn','b1_err_pct_swfn','run_s'};
fid = fopen(C.out, 'w');
fprintf(fid, '%s', hdr{1}); fprintf(fid, ',%s', hdr{2:end}); fprintf(fid, '\n');
for r = 1:size(rows,1)
    fprintf(fid, '%s,%s', rows{r,1}, rows{r,2}); fprintf(fid, ',%.6g', rows{r,3:end}); fprintf(fid, '\n');
end
fclose(fid);
fprintf('\nwrote %s  (%d runs)\n', C.out, size(rows,1));
T = struct('hdr', {hdr}, 'rows', {rows});
end

% =========================================================================
function s = setf(s, varargin)
for i = 1:2:numel(varargin), s.(varargin{i}) = varargin{i+1}; end
end

function X = read_cases(f)
%READ_CASES  Angle sets of the cases file and their analytic THD49.
txt = fileread(f);  L = regexp(txt, '\r?\n', 'split');  L = L(~cellfun(@isempty, strtrim(L)));
H = regexp(strtrim(L{1}), ',', 'split');
id = find(strcmp(H,'design'));  im = find(strcmp(H,'Ma'));  it = find(strcmp(H,'THD49_pct')); %#ok<NASGU>
tc = find(strncmp(H,'theta',5));
n = numel(L) - 1;
X.design = cell(n,1); X.Ma = zeros(n,1); X.THD49 = nan(n,1); X.theta = zeros(n, numel(tc));
F = plm_core();  P = F.defaults(numel(tc), 49);
for r = 1:n
    c = regexp(L{r+1}, ',', 'split');
    X.design{r} = c{id};  X.Ma(r) = str2double(c{im});
    X.theta(r,:) = str2double(c(tc));
    R = F.report(deg2rad(X.theta(r,:)), X.Ma(r), P);
    X.THD49(r) = 100*R.THD49;                      % analytic value, from the angles
end
end

function [t, v] = get_signal(out, name)
try
    s = out.get(name);
catch
    s = evalin('base', name);
end
if isa(s, 'timeseries')
    t = s.Time;  v = squeeze(s.Data);
elseif isstruct(s) && isfield(s, 'time')
    t = s.time;  v = s.signals.values;
else
    error('nonideal_simulink_runs:log', 'log %s as Timeseries or Structure With Time', name);
end
t = t(:);  v = v(:,1);
end

function M = spectrum_metrics(t, v, C)
%SPECTRUM_METRICS  THD49, max q_n, fundamental error, largest even order and
%   out-of-band content of v_ab over the last C.cycles periods.
T1 = 1/C.f1;  t1 = t(end);  t0 = t1 - C.cycles*T1;
tu = t0 + (0:C.cycles*C.Nfft-1).'*(T1/C.Nfft);
[tt, iu] = unique(t);  vu = interp1(tt, v(iu), tu, 'previous');   % piecewise-constant PWM
Vh = abs(fft(vu)) / (numel(vu)/2);
k1 = C.cycles + 1;                              % bin of the fundamental
V1 = Vh(k1);
ORD = 5:2:49; ORD = ORD(mod(ORD,3)~=0);
Lim = zeros(size(ORD));
for i = 1:numel(ORD)
    n = ORD(i);
    switch n
        case 5, Lim(i)=0.06; case 7, Lim(i)=0.05; case 11, Lim(i)=0.035; case 13, Lim(i)=0.03;
        case 17, Lim(i)=0.02; case {19,23,25}, Lim(i)=0.015;
        otherwise, Lim(i) = (0.2 + 32.5/n)/100;
    end
end
Vn = Vh(C.cycles*ORD + 1);
M.THD49 = 100*sqrt(sum(Vn.^2))/V1;
M.maxq  = max((Vn(:).'/V1) ./ Lim);
M.b1_err_pct = 100*(V1/(sqrt(3)*C.Vdc/2) - evalin('base', C.ma_var)) / evalin('base', C.ma_var);
ev = 2:2:50;  M.even_max = 100*max(Vh(C.cycles*ev + 1))/V1;
allb = 2:50;  oob = setdiff(allb, [ORD, 3:3:50]);   % non-triplen orders outside the counted set
M.oob = 100*sqrt(sum(Vh(C.cycles*oob + 1).^2))/V1;
end
