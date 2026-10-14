function OUT = gwo_sqp_sweep(CFG)
%GWO_SQP_SWEEP  GWO-SQP baseline: grey wolf optimiser [27] for exploration
%   followed by fmincon SQP refinement, the explore-then-refine class of [7]
%   with SQP in place of the Newton-Raphson step of the square SHE system
%   (paper, Sections 5.1 and Table 2).
%
%   OUT = gwo_sqp_sweep(CFG)
%
%   Per operating point (allowance B = NFE_TARGET):
%     1. N_RESTART independent GWO runs of POP wolves, initialised uniformly in
%        the box [theta_min, theta_max]^N as in [7]; the coefficient a decreases
%        linearly from 2 to 0. The swarm phase stops at SWARM_SHARE*B.
%        Fitness: baseline_common obj_batch (penalised THD49, 1 NFE per wolf).
%     2. After each run the leader is refined by SQP (baseline_common refine).
%     3. DRAIN: the remaining wolves of all runs are refined best-first until
%        (1 - POLISH_SHARE)*B evaluations have been spent.
%     4. The best archived solution is polished by SQP within B, followed by up
%        to N_COMPLIANT compliance-constrained refinements within B.
%   The reported angle set is the polished design of step 4.
%
%   CFG fields (defaults in brackets)
%     Ma_list           operating points                     [1.20:-0.05:0.20]
%     NFE_TARGET        allowance per operating point B       [20000]
%     seed              random seed s; point i uses rng(s + i) [1]
%     SPACE             'box' (as in [7], used in the paper) | 'z' (order-
%                       statistic coordinates)                ['box']
%     POP               wolves per run                        [30]
%     N_RESTART         GWO runs per operating point          [8]
%     SWARM_SHARE       share of B for the swarm phase        [0.50]
%     DRAIN             refine the remaining wolves           [true]
%     POLISH_SHARE      share of B held back for the polish   [0.10]
%     RHO               penalty weight of the objective       [30]
%     minGap_deg        minimum gap between angles, degrees     [0.25]
%     ARCH_SEP_DEG      minimum mean |d theta| between archived minima [1.5]
%     REPORT_COMPLIANT  run the compliance-constrained refinements [true]
%     N_COMPLIANT       number of those refinements           [6]
%     verbose           print one line per operating point    [false]
%
%   OUT fields (one row per operating point that returned a solution)
%     Ma, theta_deg (reported design, degrees), theta_compliant_deg,
%     n_minima (distinct accepted minima), NFE (evaluations charged at that
%     point), rows / row_names (per-point diagnostics), time_s.
%   The evaluation counter is the global NFE_COUNT; it is reset at every
%   operating point.
%
%   Requires baseline_common.m and the Optimization Toolbox.

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
C.SPACE            = 'box';
C.POP              = 30;
C.SWARM_SHARE      = 0.50;
C.N_RESTART        = 8;
C.DRAIN            = true;
C.RHO              = 30;
C.minGap_deg       = 0.25;
C.ARCH_SEP_DEG     = 1.5;
C.REPORT_COMPLIANT = true;
C.N_COMPLIANT      = 6;
C.POLISH_SHARE     = 0.10;
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
CAP_SWARM = C.SWARM_SHARE * TGT;
CAP_MAIN  = (1 - C.POLISH_SHARE) * TGT;

sqp_opts = optimoptions('fmincon','Algorithm','sqp','Display','none', ...
    'MaxIterations',200,'MaxFunctionEvaluations',3e4, ...
    'OptimalityTolerance',1e-12,'StepTolerance',1e-12, ...
    'ConstraintTolerance',1e-10, ...
    'SpecifyObjectiveGradient',true, ...
    'SpecifyConstraintGradient',true);
pol = optimoptions(sqp_opts,'MaxIterations',600, ...
    'OptimalityTolerance',1e-14,'StepTolerance',1e-14);

ROWS = []; XALL = []; XCALL = []; t_all = tic;

for mi = 1:numel(C.Ma_list)
    Ma = C.Ma_list(mi);
    rng(seed + mi);
    NFE_COUNT = 0;  t_ma = tic;

    ARCH = zeros(0,N); ATHD = [];
    nSwarm = 0; nRef = 0; nRefOK = 0;
    PACK = zeros(0,N); PACKF = [];    % final population of every run

    perRun = CAP_SWARM / max(1, C.N_RESTART);
    iters  = max(10, floor(perRun / C.POP));

    for rs = 1:C.N_RESTART
        if NFE_COUNT >= CAP_SWARM, break; end
        capRun = min(CAP_SWARM, NFE_COUNT + perRun);

        % ---------------- grey wolf optimiser ----------------
        if strcmpi(C.SPACE,'z')
            Wz = rand(C.POP, N);                         % positions in (0,1)^N
            W  = zeros(C.POP, N);
            for q = 1:C.POP, W(q,:) = F.z2theta(Wz(q,:), P); end
        else
            W  = P.lb + (P.ub - P.lb) * rand(C.POP, N);  % box, as in [7]
            Wz = [];
        end
        fit = swarm_fitness(W, Ma, P, F);
        nSwarm = nSwarm + C.POP;
        [fit, ordf] = sort(fit); W = W(ordf,:);
        if ~isempty(Wz), Wz = Wz(ordf,:); end
        Xa = W(1,:); Xb = W(2,:); Xd = W(3,:);
        Fa = fit(1);

        for it = 1:iters
            if NFE_COUNT >= capRun, break; end
            aC = 2 - 2*(it/iters);
            for q = 1:C.POP
                if strcmpi(C.SPACE,'z')
                    pz = Wz(q,:);
                    za = theta2z_local(Xa, P); zb = theta2z_local(Xb, P); zd = theta2z_local(Xd, P);
                    z1 = gwo_step(pz, za, aC, N);
                    z2 = gwo_step(pz, zb, aC, N);
                    z3 = gwo_step(pz, zd, aC, N);
                    pz = min(max((z1+z2+z3)/3, 0), 1);
                    Wz(q,:) = pz;  W(q,:) = F.z2theta(pz, P);
                else
                    p  = W(q,:);
                    p1 = gwo_step(p, Xa, aC, N);
                    p2 = gwo_step(p, Xb, aC, N);
                    p3 = gwo_step(p, Xd, aC, N);
                    W(q,:) = min(max((p1+p2+p3)/3, P.lb), P.ub);
                end
            end
            fit = swarm_fitness(W, Ma, P, F);
            nSwarm = nSwarm + C.POP;
            [fit, ordf] = sort(fit); W = W(ordf,:);
            if ~isempty(Wz), Wz = Wz(ordf,:); end
            if fit(1) < Fa, Fa = fit(1); Xa = W(1,:); end
            Xb = W(2,:); Xd = W(3,:);
        end

        % ---------------- SQP refinement of the leader ----------------
        PACK  = [PACK; W];  PACKF = [PACKF; fit(:)]; %#ok<AGROW>
        x0 = sort(min(max(Xa, P.lb), P.ub));
        [a, ok, R] = F.refine(x0, Ma, P, F.cap_options(sqp_opts, CAP_MAIN));
        nRef = nRef + 1;  nRefOK = nRefOK + double(ok);
        [ARCH, ATHD] = archive_insert(ARCH, ATHD, a, R.THD49, ok, SEPR);
    end

    % ---------------- DRAIN: refine the remaining wolves best-first --------
    if C.DRAIN && ~isempty(PACK)
        [~, op] = sort(PACKF);  PACK = PACK(op,:);
        for q = 1:size(PACK,1)
            if NFE_COUNT >= CAP_MAIN, break; end
            x0 = sort(min(max(PACK(q,:), P.lb), P.ub));
            [a, ok, R] = F.refine(x0, Ma, P, F.cap_options(sqp_opts, CAP_MAIN));
            nRef = nRef + 1;  nRefOK = nRefOK + double(ok);
            [ARCH, ATHD] = archive_insert(ARCH, ATHD, a, R.THD49, ok, SEPR);
        end
    end

    if isempty(ARCH)
        if C.verbose
            fprintf('  Ma=%.2f: no feasible solution (%d refinements, %d ok).\n', ...
                Ma, nRef, nRefOK);
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
    nfeRef = max(0, NFE_COUNT - nSwarm);   % evaluations spent on refinement
    ROWS = [ROWS; Ma, 100*Rb.THD49, Rb.b1, Rb.b1err, Rb.viol, Rb.MAXQ, ...
            double(Rb.compliant), 100*Rc.THD49, Rc.MAXQ, double(okc), ...
            size(ARCH,1), nRef, nRefOK, Rb.gapviol, nSwarm, NFE_COUNT, ...
            NFE_COUNT/TGT, toc(t_ma), 100*SPX.THD99, 100*SPX.WTHD49, ...
            100*SPX.WTHD999, 100*SPX.OOB_max, SPX.OOB_order, 100*SPX.b3_ratio, ...
            nfeRef]; %#ok<AGROW>
    XALL  = [XALL;  rad2deg(xb)]; %#ok<AGROW>
    XCALL = [XCALL; rad2deg(xc)]; %#ok<AGROW>

    if C.verbose
        fprintf(['  Ma=%.2f  THD49=%.4f%%  MAXQ=%.3f | %d minima from %d refinements ' ...
                 '| NFE=%d (%.2fx) = %d swarm + %d refinement  %.1fs\n'], ...
            Ma, 100*Rb.THD49, Rb.MAXQ, size(ARCH,1), nRef, NFE_COUNT, ...
            NFE_COUNT/TGT, nSwarm, nfeRef, toc(t_ma));
    end
end

OUT.row_names = {'Ma','THD49_pct','b1','b1_err','viol','MAXQ','compliant', ...
      'THD49_compliant_pct','MAXQ_compliant','compliant_found', ...
      'n_distinct_minima','n_refines','n_refines_ok','gap_violation', ...
      'n_swarm_evals','NFE','NFE_over_target','time_s', ...
      'THD99_pct','WTHD49_pct','WTHD999_pct','worst_harm_above49_pct', ...
      'worst_harm_order','b3_pct','NFE_refine'};
OUT.rows = ROWS;
if isempty(ROWS)
    OUT.Ma = zeros(0,1); OUT.n_minima = zeros(0,1); OUT.NFE = zeros(0,1);
else
    OUT.Ma = ROWS(:,1); OUT.n_minima = ROWS(:,11); OUT.NFE = ROWS(:,16);
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

function f = swarm_fitness(W, Ma, P, F)
%SWARM_FITNESS  Penalised objective of every wolf, one NFE per wolf.
f = F.obj_batch(W, Ma, P);
f = f(:);
end

function y = gwo_step(p, leader, aC, N)
%GWO_STEP  One encircling update toward a leader [27].
A = 2*aC*rand(1,N) - aC;
C = 2*rand(1,N);
D = abs(C .* leader - p);
y = leader - A .* D;
end

function z = theta2z_local(a, P)
%THETA2Z_LOCAL  Ordered angles -> (0,1)^N for SPACE = 'z' (linear map).
a = sort(a(:).');
Sspan = P.ub - P.lb - (P.N-1)*P.minGap;
z = (a - P.lb - P.minGap*(0:P.N-1)) / max(Sspan, eps);
z = min(max(z, 0), 1);
end
