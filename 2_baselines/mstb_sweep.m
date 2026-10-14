function OUT = mstb_sweep(CFG)
%MSTB_SWEEP  MS-TB baseline: MultiStart of the MATLAB Global Optimization
%   Toolbox [28] with fmincon SQP as the local solver (paper, Section 5.1 and
%   Table 2).
%
%   OUT = mstb_sweep(CFG)
%
%   Per operating point (allowance B = NFE_TARGET):
%     1. max(20, floor(0.88*B/NFE_PER_SOLVE)) start points drawn uniformly
%        from the feasible set (order-statistic map) and passed to MultiStart
%        as a CustomStartPointSet; every local solve is capped at
%        MAXFE_PER_SOLVE evaluations, and an OutputFcn aborts the remaining
%        solves once RUN_SHARE*B evaluations have been charged.
%     2. Returned minima are re-checked with the same acceptance rule and
%        archive separation as GWO-SQP.
%     3. The best minimum is polished by SQP within B, followed by up to
%        N_COMPLIANT compliance-constrained refinements within B.
%   The reported angle set is the polished design of step 3.
%   Objective and constraints are those of GWO-SQP (baseline_common).
%
%   CFG fields (defaults in brackets)
%     Ma_list           operating points                       [1.20:-0.05:0.20]
%     NFE_TARGET        allowance per operating point B         [20000]
%     seed              random seed s; point i uses rng(s + i)   [1]
%     NFE_PER_SOLVE     assumed cost of one local solve, used only to size
%                       the number of starts                    [505]
%     MAXFE_PER_SOLVE   evaluation cap of one local solve       [3000]
%     RUN_SHARE         share of B after which MultiStart is stopped [0.88]
%     RHO               penalty weight of the objective         [30]
%     minGap_deg        minimum gap between angles, degrees     [0.25]
%     ARCH_SEP_DEG      minimum mean |d theta| between archived minima [1.5]
%     REPORT_COMPLIANT  run the compliance-constrained refinements [true]
%     N_COMPLIANT       number of those refinements             [6]
%     verbose           print one line per operating point      [false]
%
%   OUT fields (one row per operating point that returned a solution)
%     Ma, theta_deg (reported design, degrees), theta_compliant_deg,
%     n_minima, NFE (evaluations charged at that point), rows / row_names
%     (per-point diagnostics), time_s.
%   The evaluation counter is the global NFE_COUNT; it is reset at every
%   operating point.
%
%   Requires baseline_common.m, the Optimization Toolbox and the Global
%   Optimization Toolbox.

if nargin < 1, CFG = struct(); end
seed = 1;  if isfield(CFG, 'seed') && ~isempty(CFG.seed), seed = CFG.seed; end
rng(seed, 'twister');

global NFE_COUNT %#ok<GVMIS>
NFE_COUNT = 0;

F = baseline_common();
P = F.defaults();

% ---------------- configuration (Table 2) ----------------
C.Ma_list          = 1.20:-0.05:0.20;
C.NFE_TARGET       = 20000;
C.NFE_PER_SOLVE    = 505;
C.MAXFE_PER_SOLVE  = 3e3;
C.RUN_SHARE        = 0.88;
C.RHO              = 30;
C.minGap_deg       = 0.25;
C.ARCH_SEP_DEG     = 1.5;
C.REPORT_COMPLIANT = true;
C.N_COMPLIANT      = 6;
C.verbose          = false;
fo = fieldnames(CFG);
for i = 1:numel(fo)
    if isfield(C, fo{i}), C.(fo{i}) = CFG.(fo{i}); end
end

P.RHO = C.RHO;
P.minGap = deg2rad(C.minGap_deg);
N    = P.N;
SEPR = deg2rad(C.ARCH_SEP_DEG);
TGT  = C.NFE_TARGET;

if isempty(ver('globaloptim'))
    error('mstb_sweep:toolbox', 'MultiStart needs the Global Optimization Toolbox.');
end

% Local solver options; OutputFcn stops every local solve once RUN_SHARE*B
% evaluations have been charged at this operating point.
lopts = optimoptions('fmincon','Algorithm','sqp','Display','off', ...
    'MaxIterations',200,'MaxFunctionEvaluations',C.MAXFE_PER_SOLVE, ...
    'OptimalityTolerance',1e-12,'StepTolerance',1e-12, ...
    'ConstraintTolerance',1e-10, ...
    'SpecifyObjectiveGradient',true, ...
    'SpecifyConstraintGradient',true);
CAP_RUN = C.RUN_SHARE * TGT;
lopts = optimoptions(lopts, 'OutputFcn', @(zz,ov,st) nfe_guard(CAP_RUN));
% Polish options: no OutputFcn; each call is bounded by cap_options instead.
pol = optimoptions(lopts,'OutputFcn',[], ...
    'MaxIterations',600,'MaxFunctionEvaluations',3e4, ...
    'OptimalityTolerance',1e-14,'StepTolerance',1e-14);

[Agap, bgap] = F.gap_A(N, P);
LB = P.lb*ones(1,N);  UB = P.ub*ones(1,N);

ROWS = []; XALL = []; XCALL = []; t_all = tic;

for mi = 1:numel(C.Ma_list)
    Ma = C.Ma_list(mi);
    rng(seed + mi);
    NFE_COUNT = 0;  t_ma = tic;

    fun = @(z) F.obj(z, Ma, P);
    con = @(z) F.nonlcon(z, Ma, P, false);

    nStart = max(20, floor(0.88*TGT / C.NFE_PER_SOLVE));

    prob = createOptimProblem('fmincon','objective',fun, ...
        'x0', F.z2theta(rand(1,N), P), 'Aineq',Agap,'bineq',bgap, ...
        'lb',LB,'ub',UB,'nonlcon',con,'options',lopts);

    ARCH = zeros(0,N); ATHD = [];
    try
        ms = MultiStart('Display','off','UseParallel',false, ...
                        'StartPointsToRun','bounds');
        SP = zeros(nStart, N);
        for q = 1:nStart, SP(q,:) = F.z2theta(rand(1,N), P); end
        [~,~,~, sout, sols] = run(ms, prob, CustomStartPointSet(SP));
    catch ME
        fprintf('  Ma=%.2f: MultiStart FAILED (%s)\n', Ma, ME.message);
        sols = []; sout = struct('localSolverTotal',0,'localSolverSuccess',0);
    end

    % ---- re-check every returned minimum with the shared archive rule ----
    for q = 1:numel(sols)
        a = min(max(sort(sols(q).X(:).'), P.lb), P.ub);
        R = F.report(a, Ma, P);
        ok = (R.viol < P.FEAS_VIOL) && (R.b1err < P.FEAS_B1);
        [ARCH, ATHD] = archive_insert(ARCH, ATHD, a, R.THD49, ok, SEPR);
    end

    if isempty(ARCH)
        if C.verbose
            fprintf('  Ma=%.2f: no feasible solution (%d raw solutions).\n', Ma, numel(sols));
        end
        continue;
    end
    [ATHD, o] = sort(ATHD);  ARCH = ARCH(o,:); %#ok<ASGLU>

    % ---------------- polish and compliance-constrained refinements --------
    [xb, okb, Rb] = F.refine(ARCH(1,:), Ma, P, F.cap_options(pol, TGT));
    if ~okb, xb = ARCH(1,:); Rb = F.report(xb, Ma, P); end

    Rc = Rb;  xc = xb;  okc = Rb.compliant;  bestc = inf;
    if C.REPORT_COMPLIANT
        for j = 1:min(C.N_COMPLIANT, size(ARCH,1))
            if NFE_COUNT >= TGT, break; end
            [a2, ok2, R2] = F.refine_compliant(ARCH(j,:), Ma, P, F.cap_options(pol, TGT));
            if ok2 && R2.THD49 < bestc
                bestc = R2.THD49; Rc = R2; okc = true; xc = a2;
            end
        end
    end

    SPX = F.spectrum(xb, 999);
    nLoc = 0; nLocOK = 0;
    if isstruct(sout)
        if isfield(sout,'localSolverTotal'),   nLoc   = sout.localSolverTotal;   end
        if isfield(sout,'localSolverSuccess'), nLocOK = sout.localSolverSuccess; end
    end
    ROWS = [ROWS; Ma, 100*Rb.THD49, Rb.b1, Rb.b1err, Rb.viol, Rb.MAXQ, ...
            double(Rb.compliant), 100*Rc.THD49, Rc.MAXQ, double(okc), ...
            size(ARCH,1), numel(sols), nLoc, nLocOK, NFE_COUNT, NFE_COUNT/TGT, ...
            toc(t_ma), 100*SPX.THD99, 100*SPX.WTHD49, 100*SPX.WTHD999, ...
            100*SPX.OOB_max, SPX.OOB_order, 100*SPX.b3_ratio]; %#ok<AGROW>
    XALL  = [XALL;  rad2deg(xb)]; %#ok<AGROW>
    XCALL = [XCALL; rad2deg(xc)]; %#ok<AGROW>

    if C.verbose
        fprintf(['  Ma=%.2f  THD49=%.4f%%  MAXQ=%.3f | %d kept of %d returned ' ...
                 '| %d local solves (%d ok) | NFE=%d (%.2fx)  %.1fs\n'], ...
            Ma, 100*Rb.THD49, Rb.MAXQ, size(ARCH,1), numel(sols), nLoc, nLocOK, ...
            NFE_COUNT, NFE_COUNT/TGT, toc(t_ma));
    end
end

OUT.row_names = {'Ma','THD49_pct','b1','b1_err','viol','MAXQ','compliant', ...
      'THD49_compliant_pct','MAXQ_compliant','compliant_found', ...
      'n_distinct_minima','n_solver_solutions','n_local_solves', ...
      'n_local_solves_ok','NFE','NFE_over_target','time_s', ...
      'THD99_pct','WTHD49_pct','WTHD999_pct','worst_harm_above49_pct', ...
      'worst_harm_order','b3_pct'};
OUT.rows = ROWS;
if isempty(ROWS)
    OUT.Ma = zeros(0,1); OUT.n_minima = zeros(0,1); OUT.NFE = zeros(0,1);
else
    OUT.Ma = ROWS(:,1); OUT.n_minima = ROWS(:,11); OUT.NFE = ROWS(:,15);
end
OUT.theta_deg = XALL;
OUT.theta_compliant_deg = XCALL;
OUT.time_s = toc(t_all);
end

% =========================================================================
function [ARCH, ATHD] = archive_insert(ARCH, ATHD, a, thd, ok, sep)
%ARCHIVE_INSERT  Keep accepted minima at least sep (mean |d theta|) apart;
%   a closer one replaces its neighbour only if it has lower THD49.
if ok
    if isempty(ARCH) || min(mean(abs(ARCH - a),2)) > sep
        ARCH(end+1,:) = a;  ATHD(end+1,1) = thd;
    else
        [~, jj] = min(mean(abs(ARCH - a),2));
        if thd < ATHD(jj), ARCH(jj,:) = a; ATHD(jj) = thd; end
    end
end
end

function stop = nfe_guard(ceilNFE)
%NFE_GUARD  fmincon OutputFcn: stop the current local solve once the
%   evaluations charged at this operating point reach ceilNFE.
global NFE_COUNT %#ok<GVMIS>
stop = NFE_COUNT >= ceilNFE;
end
