function T = run_benchmark(CFG)
%RUN_BENCHMARK  Matched-evaluation benchmark harness (paper, Section 4).
%
%   csvfile = run_benchmark(CFG)
%
%   Runs every requested arm at every budget and seed over the operating-point
%   grid, passes the returned angles through evaluate_candidates (one
%   restoration and acceptance rule for every arm) and appends one row per
%   (arm, budget, seed, Ma) to a long-format CSV.
%
%   CFG fields (defaults in brackets)
%     arms      cell array of arm keys, see ARM REGISTRY below
%               [{'gn_hybrid','gn_multistart','sqp_multistart','multistart_tb','fssa'}]
%     budgets   evaluation allowances per operating point
%               [500 1000 2000 4000 8000 20000]
%     seeds     seed list                                  [1:10]
%     Ma_list   operating points, solved in this order      [1.20:-0.05:0.20]
%     N, nmax   number of angles, highest counted order     [15, 49]
%     minGap_deg  minimum gap between angles, degrees       [0.25]
%                 (0.50 for the minimum-gap study of Section 6.1)
%     out       output CSV name ('' = run_benchmark_<date>.csv) ['']
%     verbose   one line per run                            [true]
%     xlsx      also write <out>.xlsx (csv_to_xlsx)         [true]
%     summary   also write <out>_summary.xlsx (summarize_benchmark) [true]
%
%   Arm keys (the keys are stored in the CSV and are fixed; the names used in
%   the paper are given on the right):
%     gn_hybrid             PLM-C    proposed: PLM corrector + bidirectional continuation
%     gn_multistart         PLM-MS   PLM corrector, multi-start
%     sqp_multistart        SQP-MS   fmincon SQP corrector, multi-start
%     sqp_hybrid            SQP-C    fmincon SQP corrector + bidirectional continuation
%     gn_hybrid_zeroth_fwd  forward-only continuation of [2] with the PLM corrector
%                           (Table 5, row 2)
%     fssa                  FSSA     simulated annealing with continuation [2]
%     fssa_grid             FSSA with the weighted grid-code objective of [2]
%                           (Section 5.6)
%     gwo_sqp_box           GWO-SQP  grey wolf optimiser + SQP [7, 27]
%     multistart_tb         MS-TB    MATLAB MultiStart [28]
%
%   CSV columns: arm, family, engine, strategy, use_spwm, nfe_target,
%   nfe_achieved, nfe_ratio, seed, Ma, accepted, THD49_pct, THD99_pct,
%   WTHD49_pct, WTHD999_pct, worst_harm_above49_pct, worst_harm_order, b3_pct,
%   MAXQ, compliant, b1, b1_err, viol, n_over_limit, n_distinct_minima,
%   time_s, theta01_deg ... thetaNN_deg.
%   nfe_achieved is the evaluations per operating point actually charged
%   (Section 4); time_s is the wall-clock time of the whole sweep.
%
%   Example
%     run_benchmark(struct('arms',{{'gn_hybrid','gn_multistart'}}, ...
%                          'budgets',2000,'seeds',1:3,'out','demo.csv'));

% Shared handles for other drivers (validation export).
if nargin == 1 && (ischar(CFG) || isstring(CFG)) && strcmpi(char(CFG), '-handles')
    T = struct('arm_spec', @arm_spec, 'run_arm', @run_arm, 'append_rows', @append_rows);
    return;
end

D.arms    = {'gn_hybrid','gn_multistart','sqp_multistart','multistart_tb','fssa'};
D.budgets = [500 1000 2000 4000 8000 20000];
D.seeds   = 1:10;
D.Ma_list = 1.20:-0.05:0.20;
D.N = 15;  D.nmax = 49;  D.minGap_deg = 0.25;
D.out = '';
D.verbose = true;
D.summary = true;
D.xlsx    = true;
if nargin < 1, CFG = struct(); end
fn = fieldnames(CFG); C = D;
for i = 1:numel(fn), C.(fn{i}) = CFG.(fn{i}); end
if ischar(C.arms), C.arms = {C.arms}; end

F = plm_core();  P = F.defaults(C.N, C.nmax, C.minGap_deg);
MAS = C.Ma_list(:);  nM = numel(MAS);

if isempty(C.out)
    c = clock(); %#ok<CLOCK>
    C.out = sprintf('run_benchmark_%04d%02d%02d_%02d%02d.csv', c(1),c(2),c(3),c(4),c(5));
end

hdr = [{'arm','family','engine','strategy','use_spwm','nfe_target','nfe_achieved', ...
        'nfe_ratio','seed','Ma','accepted','THD49_pct','THD99_pct', ...
        'WTHD49_pct','WTHD999_pct','worst_harm_above49_pct','worst_harm_order', ...
        'b3_pct','MAXQ','compliant','b1','b1_err','viol','n_over_limit', ...
        'n_distinct_minima','time_s'}, ...
       arrayfun(@(k) sprintf('theta%02d_deg',k), 1:C.N, 'UniformOutput', false)];
fid = fopen(C.out, 'w');
fprintf(fid, '%s', hdr{1});
for i = 2:numel(hdr), fprintf(fid, ',%s', hdr{i}); end
fprintf(fid, '\n');
fclose(fid);

if C.verbose
    fprintf('\n================= run_benchmark =================\n');
    fprintf('arms    : %s\n', strjoin(C.arms, ', '));
    fprintf('budgets : %s\n', mat2str(C.budgets));
    fprintf('seeds   : %s\n', mat2str(C.seeds));
    fprintf('Ma      : %g : %g : %g (%d points)\n', MAS(1), MAS(2)-MAS(1), MAS(end), nM);
    fprintf('out     : %s\n', C.out);
    fprintf('acceptance (all arms): |b1 - Ma| < %.0e and violation < 1e-9\n\n', P.FEAS_B1);
end

nRow = 0;
for ai = 1:numel(C.arms)
    A = arm_spec(C.arms{ai});
    for B = C.budgets(:).'
        for sd = C.seeds(:).'
            [TH, nfe, wall, nmin, failed] = run_arm(A, MAS, B, sd, C, P);
            if failed
                fprintf('  !! %-16s B=%-7g seed=%-3d FAILED TO RUN (%s)\n', ...
                        A.name, B, sd, A.note);
                continue;
            end
            R = evaluate_candidates(TH, MAS, P, F);
            M = [repmat([A.spwm, B, nfe, nfe/max(B,eps), sd], nM, 1), MAS, ...
                 double(R.accepted), 100*R.THD49, 100*R.THD99, ...
                 100*R.WTHD49, 100*R.WTHD999, 100*R.OOB_max, R.OOB_order, ...
                 100*R.b3_ratio, R.MAXQ, double(R.compliant), R.b1, R.b1err, ...
                 R.viol, R.nviol, nmin(:), repmat(wall, nM, 1), R.theta_deg];
            append_rows(C.out, {A.name, A.family, A.engine, A.strategy}, M);
            nRow = nRow + nM;
            if C.verbose
                good = R.accepted & isfinite(R.THD49);
                fprintf(['  %-20s B=%-7g seed=%-3d  accepted %2d/%2d  ' ...
                         'median THD49 %7.4f%%  NFE %8.0f (%.2fx)  %6.1fs\n'], ...
                    A.name, B, sd, sum(good), nM, ...
                    100*nmed(R.THD49(good)), nfe, nfe/max(B,eps), wall);
            end
        end
    end
end

if C.verbose
    fprintf('\nwrote %s  (%d rows)\n', C.out, nRow);
end
if C.xlsx && nRow > 0
    try, csv_to_xlsx(C.out, 'runs'); catch e
        fprintf('  (xlsx skipped: %s -- %s is complete)\n', e.message, C.out);
    end
end
if C.summary && nRow > 0
    try
        summarize_benchmark(C.out);
    catch e
        fprintf('\n  (summary skipped: %s)\n  rerun with summarize_benchmark(''%s'')\n', ...
                e.message, C.out);
    end
end
T = C.out;
end

% =========================================================================
% ARM REGISTRY
% =========================================================================
function A = arm_spec(name)
%ARM_SPEC  Arm key -> driver, settings and the labels written to the CSV.
A.name = name;  A.ovr = struct();  A.note = '';
A.spwm = 0;                       % no SPWM-seeded starts in any arm
switch lower(name)
  % ---- plm_sweep: PLM / SQP corrector x multi-start / continuation ------
  case 'gn_multistart'                                          % PLM-MS
    A.kind='plm_sweep'; A.family='proposed'; A.engine='projected-LM'; A.strategy='multistart';
    A.ovr = struct('method','multistart','engine','plm');
  case 'gn_hybrid'                                              % PLM-C
    A.kind='plm_sweep'; A.family='proposed'; A.engine='projected-LM'; A.strategy='continuation';
    A.ovr = struct('method','continuation','engine','plm');
  case 'sqp_multistart'                                         % SQP-MS
    A.kind='plm_sweep'; A.family='proposed'; A.engine='fmincon-SQP'; A.strategy='multistart';
    A.ovr = struct('method','multistart','engine','sqp'); A.note='needs Optimization Toolbox';
  case 'sqp_hybrid'                                             % SQP-C
    A.kind='plm_sweep'; A.family='proposed'; A.engine='fmincon-SQP'; A.strategy='continuation';
    A.ovr = struct('method','continuation','engine','sqp'); A.note='needs Optimization Toolbox';
  case 'gn_hybrid_zeroth_fwd'                                   % Table 5, row 2
    A.kind='plm_sweep'; A.family='ablation'; A.engine='projected-LM';
    A.strategy='continuation-zeroth-1pass';
    A.ovr = struct('method','continuation','engine','plm','PASSES',1);
  % ---- published SHM method with continuation [2] -----------------------
  case 'fssa'                                                   % FSSA
    A.kind='fssa'; A.family='published'; A.engine='simulated-annealing';
    A.strategy='continuation'; A.ovr = struct('OBJ','ls','RESTART_ON_FAIL',false, ...
                                              'SWEEP_DIR','ascending');
  case 'fssa_grid'                                              % FSSA, objective of [2] (Sec. 5.6)
    A.kind='fssa'; A.family='published'; A.engine='simulated-annealing';
    A.strategy='continuation-gridobj'; A.ovr = struct('OBJ','grid','RESTART_ON_FAIL',false, ...
                                                      'SWEEP_DIR','ascending');
  % ---- toolbox-based baselines -------------------------------------------
  case 'gwo_sqp_box'                                            % GWO-SQP
    A.kind='gwo'; A.family='modern';
    A.engine='fmincon-SQP'; A.strategy='gwo-hybrid-boxspace';
    A.ovr = struct('SPACE','box'); A.note='needs Optimization Toolbox';
  case 'multistart_tb'                                          % MS-TB
    A.kind='mstb'; A.family='toolbox';
    A.engine='fmincon-SQP'; A.strategy='multistart';
    A.note='needs Optimization + Global Optimization Toolboxes';
  otherwise
    error('run_benchmark:arm', 'unknown arm "%s"', name);
end
end

% =========================================================================
% RUN ONE ARM AND COLLECT ITS ANGLES
% =========================================================================
function [TH, nfe, wall, nmin, failed] = run_arm(A, MAS, B, seed, C, P)
%RUN_ARM  Run one arm over the grid at one budget and seed. Returns the
%   arm's angles (radians, NaN rows where it returned nothing), the achieved
%   NFE per operating point, the wall-clock time and the number of distinct
%   minima per point.
global NFE_COUNT %#ok<GVMIS>
nM = numel(MAS);
TH = nan(nM, P.N);  nmin = nan(nM,1);  nfe = NaN;  wall = NaN;  failed = false;
NFE_COUNT = 0;

switch A.kind
  case 'plm_sweep'
    R = A.ovr;
    R.Ma_list = MAS(:).';  R.NFE_TARGET = B;  R.seed = seed;
    R.N = C.N;  R.nmax = C.nmax;  R.minGap_deg = C.minGap_deg;  R.verbose = false;
    t0 = tic;
    try, O = plm_sweep(R); catch, failed = true; return; end
    wall = toc(t0);
    TH   = deg2rad(O.theta_deg);
    nfe  = mean(O.NFE);
    nmin = O.n_minima;

  case 'fssa'
    R = A.ovr;
    R.Ma_list = MAS(:).';  R.NFE_TARGET = B;  R.seed = seed;
    R.N = C.N;  R.nmax = C.nmax;  R.minGap_deg = C.minGap_deg;  R.verbose = false;
    t0 = tic;
    try, O = fssa_sweep(R); catch, failed = true; return; end
    wall = toc(t0);
    TH   = deg2rad(O.theta_deg);
    nfe  = mean(O.NFE);
    nmin = O.n_minima;

  case {'gwo', 'mstb'}
    R = A.ovr;
    R.Ma_list = MAS(:).';  R.NFE_TARGET = B;  R.seed = seed;  R.minGap_deg = C.minGap_deg;
    t0 = tic;
    try
        if strcmp(A.kind, 'gwo'), O = gwo_sqp_sweep(R); else, O = mstb_sweep(R); end
    catch
        failed = true; return;
    end
    wall = toc(t0);
    % The global counter is reset at every operating point by these arms, so
    % after the sweep it holds the evaluations charged at the final point.
    nfe  = NFE_COUNT;
    for k = 1:numel(O.Ma)
        j = find(abs(MAS - O.Ma(k)) < 1e-9, 1);
        if isempty(j), continue; end
        TH(j,:) = deg2rad(O.theta_deg(k,:));
        nmin(j) = O.n_minima(k);
    end
end
end

% =========================================================================
function append_rows(fn, labels, M)
%APPEND_ROWS  Append numeric rows with four text labels; NaN -> empty field.
fid = fopen(fn, 'a');
for r = 1:size(M,1)
    fprintf(fid, '%s,%s,%s,%s', labels{1}, labels{2}, labels{3}, labels{4});
    for c = 1:size(M,2)
        if isnan(M(r,c)), fprintf(fid, ',');
        else,             fprintf(fid, ',%.10g', M(r,c));
        end
    end
    fprintf(fid, '\n');
end
fclose(fid);
end

function y = nmed(x)
x = sort(x(isfinite(x)));  n = numel(x);
if n == 0, y = NaN; elseif mod(n,2), y = x((n+1)/2); else, y = 0.5*(x(n/2)+x(n/2+1)); end
end
