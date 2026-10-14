function OUT = plm_sweep(CFG)
%PLM_SWEEP  Solve the 15-angle SHM design over a modulation-index grid with
%   the PLM corrector (or fmincon SQP) under a fixed per-point NFE budget.
%
%   OUT = plm_sweep(CFG)
%
%   One driver produces four arms of the paper, selected by two fields:
%
%                         CFG.engine = 'plm'        CFG.engine = 'sqp'
%     method = 'continuation'  PLM-C (proposed)      SQP-C  (ablation, Table 5)
%     method = 'multistart'    PLM-MS                SQP-MS
%
%   MULTISTART    independent multi-start at every operating point; starts
%                 are drawn uniformly from the feasible set (Section 3.1).
%   CONTINUATION  bidirectional continuation in Ma (Section 3.3):
%                 pass 1 (forward, F_FWD of the budget) warm-starts each point
%                 from the K_CARRY best solutions of the previous point;
%                 pass 2 (reverse, the remainder) warm-starts from its own
%                 carried solutions and K_POOL draws from an archive pooled over
%                 the whole grid. Warm starting stops after CARRY_STALL batches
%                 without improvement or at CARRY_SHARE of the pass, and the
%                 rest of the allowance goes to fresh multi-start.
%                 PASSES = 1 gives the published forward-only scheme (Table 5,
%                 row 2), which then receives the whole allowance.
%
%   CFG fields (defaults in brackets; Table 2 of the paper)
%     method        'continuation' | 'multistart'              ['continuation']
%     engine        'plm' | 'sqp'                              ['plm']
%     Ma_list       modulation indices, solved in the given order [1.20:-0.05:0.20]
%     NFE_TARGET    evaluation allowance per operating point    [20000]
%     seed          random seed s; point i uses 1000*s + i (forward/multistart)
%                   and 1000*s + 500 + j (reverse pass, j-th point)    [1]
%     N, nmax       number of angles and highest counted order  [15, 49]
%     minGap_deg    minimum gap delta between angles, degrees   [0.25]
%     K_CARRY       solutions carried between points           [24]
%     K_POOL        draws from the pooled archive in pass 2     [24]
%     F_FWD         share of the allowance for pass 1          [0.60]
%     CARRY_SHARE   maximum share of a pass spent warm-starting [0.55]
%     CARRY_STALL   batches without gain before warm starting stops [2]
%     CARRY_BATCH   carried candidates per batch                [4]
%     PASSES        2 = forward + reverse, 1 = forward only     [2]
%     ARCH_SEP_DEG  minimum mean |d theta| between archived minima [0.75]
%     CHUNK         starts advanced together (batch width)      [16]
%     verbose       print a per-point table                     [true]
%
%   OUT fields
%     Ma, theta_deg (best accepted angles per point, degrees), THD49, MAXQ,
%     b1, b1err, viol, nviol, compliant, n_minima (distinct accepted minima),
%     NFE (achieved evaluations per point), time_s, arch (archives), CFG, P.
%
%   Example
%     OUT = plm_sweep(struct('NFE_TARGET', 4000, 'seed', 3));
%
%   Requires plm_core.m; engine 'sqp' also requires the Optimization Toolbox.
%   Results reported in the paper are produced through run_benchmark.m, which
%   passes every field explicitly and judges the returned angles with
%   evaluate_candidates.m.

cost_est(-1);                     % reset the per-engine cost estimate (see cost_est)
F = plm_core();

% ----------------------------- defaults ----------------------------------
D.method       = 'continuation';
D.engine       = 'plm';
D.Ma_list      = 1.20:-0.05:0.20;
D.NFE_TARGET   = 20000;
D.seed         = 1;
D.N            = 15;
D.nmax         = 49;
D.minGap_deg   = 0.25;
D.K_CARRY      = 24;
D.K_POOL       = 24;
D.F_FWD        = 0.60;
D.CARRY_SHARE  = 0.55;
D.CARRY_STALL  = 2;
D.CARRY_BATCH  = 4;
D.PASSES       = 2;
D.ARCH_SEP_DEG = 0.75;
D.CHUNK        = 16;
D.verbose      = true;
CFG = merge(D, CFG);

P = F.defaults(CFG.N, CFG.nmax, CFG.minGap_deg);
MAS = CFG.Ma_list(:).';
nM  = numel(MAS);
rng_set(CFG.seed);

if CFG.verbose
    fprintf('\n== plm_sweep | method=%s | engine=%s | N=%d | %d Ma points | %g NFE/Ma | seed=%d ==\n', ...
        CFG.method, CFG.engine, P.N, nM, CFG.NFE_TARGET, CFG.seed);
end

t_all = tic;
switch lower(CFG.method)
    case {'multistart', 'baseline'},   [AR, NFE] = run_multistart(F, P, MAS, CFG);
    case {'continuation', 'hybrid'},   [AR, NFE] = run_continuation(F, P, MAS, CFG);
    otherwise, error('plm_sweep:method', 'CFG.method must be continuation or multistart');
end

% ------------------------------ collect ----------------------------------
OUT.CFG = CFG; OUT.P = P; OUT.Ma = MAS(:);
nA = numel(MAS);
OUT.THD49 = nan(nA,1); OUT.MAXQ = nan(nA,1); OUT.b1err = nan(nA,1); OUT.b1 = nan(nA,1);
OUT.viol = nan(nA,1); OUT.nviol = nan(nA,1); OUT.compliant = false(nA,1);
OUT.n_minima = zeros(nA,1); OUT.theta_deg = nan(nA, P.N);
OUT.NFE = NFE(:);
for mi = 1:nA
    A = AR(mi);
    if isempty(A.f)
        warning('plm_sweep:noSolution', ...
            'Ma = %.3f: no accepted solution (%d local solves, %d converged).', ...
            MAS(mi), AR(mi).n_try, AR(mi).n_ok);
        continue;
    end
    x = A.X(1,:);
    R = F.report(x, MAS(mi), P);
    OUT.THD49(mi)=R.THD49; OUT.MAXQ(mi)=R.MAXQ; OUT.b1err(mi)=R.b1err; OUT.b1(mi)=R.b1;
    OUT.viol(mi)=R.viol; OUT.nviol(mi)=R.nviol; OUT.compliant(mi)=R.compliant;
    OUT.n_minima(mi)=numel(A.f); OUT.theta_deg(mi,:)=rad2deg(x);
end
OUT.time_s = toc(t_all);
OUT.arch = AR;

if CFG.verbose
    fprintf('  %-5s %10s %8s %6s %8s %9s\n', 'Ma', 'THD49 %', 'MAXQ', 'comp', 'minima', 'NFE');
    for mi = 1:nA
        fprintf('  %5.2f %10.4f %8.4f %6s %8d %9d\n', OUT.Ma(mi), 100*OUT.THD49(mi), ...
            OUT.MAXQ(mi), tern(OUT.compliant(mi),'yes','no'), OUT.n_minima(mi), OUT.NFE(mi));
    end
    fprintf('  --- median THD49 %.4f%% | NFE %d (%.2fx target) | %.1f s\n\n', ...
        100*median(OUT.THD49(~isnan(OUT.THD49))), sum(OUT.NFE), ...
        mean(OUT.NFE)/CFG.NFE_TARGET, OUT.time_s);
end
end

% =========================================================================
% MULTI-START  (PLM-MS, SQP-MS)
% =========================================================================
function [AR, NFE] = run_multistart(F, P, MAS, CFG)
nM = numel(MAS);
AR = archive_array(CFG.ARCH_SEP_DEG, nM);
NFE = zeros(1,nM);
for mi = 1:nM
    Ma = MAS(mi);
    rng_set(CFG.seed*1000 + mi);
    F.nfe_reset(); B = CFG.NFE_TARGET;
    guard = 0;
    while B - F.nfe() > 1.2*cost_per(CFG)
        guard = guard + 1;
        if guard > max_rounds(CFG, B), warn_stall('multistart'); break; end
        n = min(CFG.CHUNK, max(1, floor((B - F.nfe())/cost_per(CFG))));
        X0 = F.zsample(rand(n, P.N), P);
        [X, I] = local_solve(F, X0, Ma, P, CFG, 'full');
        AR(mi) = archive_add(AR(mi), X, I);
    end
    NFE(mi) = F.nfe();
end
end

% =========================================================================
% BIDIRECTIONAL CONTINUATION  (PLM-C, SQP-C)
% =========================================================================
function [AR, NFE] = run_continuation(F, P, MAS, CFG)
nM = numel(MAS);
AR = archive_array(CFG.ARCH_SEP_DEG, nM);
NFE = zeros(1,nM);
% With a single pass there is no reverse pass to reserve budget for.
if CFG.PASSES >= 2, B1 = CFG.F_FWD * CFG.NFE_TARGET; else, B1 = CFG.NFE_TARGET; end

% ---------------- pass 1: forward continuation ---------------------------
carry = [];
for mi = 1:nM
    Ma = MAS(mi);
    rng_set(CFG.seed*1000 + mi);
    F.nfe_reset();
    if ~isempty(carry)
        AR(mi) = warm_start_phase(F, AR(mi), carry, Ma, P, CFG, B1);
    end
    AR(mi) = fill_with_random_starts(F, AR(mi), Ma, P, CFG, B1);
    NFE(mi) = NFE(mi) + F.nfe();
    nb = min(CFG.K_CARRY, numel(AR(mi).f));
    carry = AR(mi).X(1:nb, :);      % zeroth-order predictor: carry unchanged
end

% ---------------- pass 2: reverse, seeded from the pooled archive --------
if CFG.PASSES < 2, return; end
pool = [];
for mi = 1:nM
    if ~isempty(AR(mi).f)
        pool = [pool; AR(mi).X(1:min(6,end), :)]; %#ok<AGROW>
    end
end
carry = [];
for j = 1:nM
    mi = nM - j + 1;   Ma = MAS(mi);
    rng_set(CFG.seed*1000 + 500 + j);
    F.nfe_reset();
    allowance = CFG.NFE_TARGET - NFE(mi);
    C = [];
    if ~isempty(carry)
        C = carry;
    end
    if ~isempty(pool)
        idx = randperm(size(pool,1), min(CFG.K_POOL, size(pool,1)));
        C = [C; pool(idx,:)];
    end
    if ~isempty(C)
        AR(mi) = warm_start_phase(F, AR(mi), C, Ma, P, CFG, allowance);
    end
    AR(mi) = fill_with_random_starts(F, AR(mi), Ma, P, CFG, allowance);
    NFE(mi) = NFE(mi) + F.nfe();
    nb = min(CFG.K_CARRY, numel(AR(mi).f));
    carry = AR(mi).X(1:nb, :);
end
end

% =========================================================================
% STAGES
% =========================================================================
function A = warm_start_phase(F, A, C, Ma, P, CFG, budget)
%WARM_START_PHASE  Refine carried candidates C in best-first batches of
%   CARRY_BATCH and stop after CARRY_STALL batches that fail to improve the
%   incumbent, or once CARRY_SHARE of the budget is spent.
n0 = F.nfe();
cap = CFG.CARRY_SHARE * budget;
if isempty(A.f)
    % cold probe first, so "improved the incumbent" is a meaningful test
    [X, I] = local_solve(F, F.zsample(rand(4, P.N), P), Ma, P, CFG, 'full');
    A = archive_add(A, X, I);
end
stall = 0; k = 1;
while k <= size(C,1) && (F.nfe() - n0) < cap && stall < CFG.CARRY_STALL
    before = archive_best(A);
    kk = min(size(C,1), k + CFG.CARRY_BATCH - 1);
    [X, I] = local_solve(F, C(k:kk,:), Ma, P, CFG, 'warm');
    A = archive_add(A, X, I);
    if archive_best(A) < before - 1e-12, stall = 0; else, stall = stall + 1; end
    k = kk + 1;
end
end

function A = fill_with_random_starts(F, A, Ma, P, CFG, budget)
%FILL_WITH_RANDOM_STARTS  Spend the remaining budget on fresh random starts.
guard = 0;
while budget - F.nfe() > 1.2*cost_per(CFG)
    guard = guard + 1;
    if guard > max_rounds(CFG, budget), warn_stall('fill_with_random_starts'); break; end
    n = min(CFG.CHUNK, max(1, floor((budget - F.nfe())/cost_per(CFG))));
    [X, I] = local_solve(F, F.zsample(rand(n,P.N),P), Ma, P, CFG, 'full');
    A = archive_add(A, X, I);
end
end

% =========================================================================
% BUDGET CONTROL
% =========================================================================
function c = cost_per(CFG)
%COST_PER  Current estimate of NFE per local solve, used to size the next
%   batch so that the achieved NFE matches the allowance.
c = cost_est(0, 0, strcmpi(CFG.engine,'sqp'));
end

function c = cost_est(update, observed, issqp)
%COST_EST  Running estimate of NFE per local solve, per engine.
%   cost_est(0,0,issqp)   read
%   cost_est(1,observed)  update the SQP estimate with a measured cost
%   cost_est(-1)          reset to the initial values (once per plm_sweep call,
%                         so that no estimate carries over between runs)
persistent c_gn c_sqp
if isempty(c_gn),  c_gn  = 132; end
if isempty(c_sqp), c_sqp = 505; end
if nargin >= 1 && isnumeric(update) && isscalar(update) && update < 0
    c_gn = 132;  c_sqp = 505;  c = NaN;  return;
end
if nargin < 3, issqp = true; end
if update && observed > 0
    if issqp, c_sqp = 0.7*c_sqp + 0.3*observed;
    else,     c_gn  = 0.7*c_gn  + 0.3*observed;
    end
end
if issqp, c = c_sqp; else, c = c_gn; end
end

function m = max_rounds(CFG, budget)
%MAX_ROUNDS  Backstop on the number of batches in one budget loop.
m = ceil(budget / max(1, CFG.CHUNK)) + 100;
end

function warn_stall(where)
warning('plm_sweep:stalled', ...
    ['%s: the local engine reported almost no evaluations; the budget loop ' ...
     'was stopped by its backstop.'], where);
end

% =========================================================================
% LOCAL ENGINES
% =========================================================================
function [X, I] = local_solve(F, X0, Ma, P, CFG, kind)
%LOCAL_SOLVE  Batched PLM corrector, or fmincon SQP one start at a time.
%   kind = 'warm' uses the warm-start weight schedule (Table 2).
if isempty(X0), X = X0; I = struct('ok',false(0,1),'thd',[],'maxq',[]); return; end
if strcmpi(CFG.engine, 'sqp')
    [X, I] = sqp_solve(F, X0, Ma, P);
else
    o = struct();
    if strcmp(kind,'warm'), o.schedule = [1e2 20; 1e5 8]; end
    [X, I] = F.solve(X0, Ma, P, o);
end
end

function [X, I] = sqp_solve(F, X0, Ma, P)
%SQP_SOLVE  fmincon SQP on min sum b_n^2 s.t. b_1 = Ma and the linear gap
%   constraints, with analytic gradients (Table 2). Every solve is finished
%   by the same Newton restoration onto {b_1 = Ma} as the PLM corrector, and
%   that work is charged.
N = P.N;
if exist('fmincon','file') ~= 2
    error('plm_sweep:nofmincon', ...
        'CFG.engine=''sqp'' needs the Optimization Toolbox. Use engine=''plm''.');
end
A = zeros(N-1,N);
for i = 1:N-1, A(i,i) = 1; A(i,i+1) = -1; end
b = -P.minGap*ones(N-1,1);
opt = optimoptions('fmincon','Algorithm','sqp','Display','none', ...
    'MaxIterations',300,'MaxFunctionEvaluations',3e4, ...
    'OptimalityTolerance',1e-12,'StepTolerance',1e-12,'ConstraintTolerance',1e-10, ...
    'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);
X = zeros(size(X0));
I.ok = false(size(X0,1),1); I.thd = inf(size(X0,1),1); I.maxq = inf(size(X0,1),1);
X0 = F.project(X0, P);
for k = 1:size(X0,1)
    x = X0(k,:);  n_k = F.nfe();
    try
        x = fmincon(@(z) sqp_obj(z,P,F), x, A, b, [], [], ...
                    P.lb*ones(1,N), P.ub*ones(1,N), @(z) sqp_con(z,Ma,P), opt);
    catch
    end
    % Newton restoration onto {b_1 = Ma} with projection, charged
    x = F.project(x(:).', P);
    for it = 1:5
        gb = (4/pi)*(cos(x)*P.s(:)) - Ma;
        F.nfe_add(1);
        if abs(gb) < 1e-15, break; end
        G1 = -(4/pi) * sin(x) .* P.s;
        x  = F.project(x - (gb./max(sum(G1.^2,2),1e-300)).*G1, P);
    end
    R = F.report(x, Ma, P);
    X(k,:) = x;
    I.ok(k)  = (R.b1err < P.FEAS_B1) && (R.viol < 1e-9);
    if I.ok(k), I.thd(k) = R.THD49; I.maxq(k) = R.MAXQ; end
    cost_est(1, F.nfe() - n_k);                % update the batch-size estimate
end
end

function [J, g] = sqp_obj(x, P, F)
%SQP_OBJ  sum b_n^2 and its analytic gradient; one NFE per call.
F.nfe_add(1);
[B, JB] = F.bnjac(x(:).', P);
J = sum(B.^2);
g = 2 * reshape(sum(B(:) .* reshape(JB, size(JB,2), size(JB,3)), 1), [], 1);
end

function [c, ceq, gc, gceq] = sqp_con(x, Ma, P)
%SQP_CON  Fundamental equality b_1 = Ma and its gradient.
x = x(:).';
c = []; gc = [];
ceq = (4/pi)*sum(P.s .* cos(x)) - Ma;
gceq = (-(4/pi) * (P.s .* sin(x))).';
end

% =========================================================================
% ARCHIVE OF DISTINCT LOCAL MINIMA (one per operating point)
% =========================================================================
function AR = archive_array(sep_deg, n)
AR = archive_new(sep_deg);
for i = 2:n, AR(i) = archive_new(sep_deg); end
end

function A = archive_new(sep_deg)
A.sep = deg2rad(sep_deg); A.cap = 300;
A.X = zeros(0,0); A.f = zeros(0,1); A.q = zeros(0,1);
A.n_try = 0; A.n_ok = 0;
end

function A = archive_add(A, X, I)
%ARCHIVE_ADD  Insert accepted solutions best-first; a solution within mean
%   |d theta| < sep of an archived one replaces it only if it is better.
if isempty(X), return; end
f = I.thd(:); q = I.maxq(:);
A.n_try = A.n_try + numel(f);
A.n_ok  = A.n_ok  + sum(isfinite(f));
[~, ord] = sort(f);
for t = ord(:).'
    if ~isfinite(f(t)), continue; end
    x = X(t,:);
    placed = false;
    if ~isempty(A.f)
        d = mean(abs(A.X - x), 2);
        [dm, j] = min(d);
        if dm < A.sep
            if f(t) < A.f(j), A.X(j,:) = x; A.f(j) = f(t); A.q(j) = q(t); end
            placed = true;
        end
    end
    if ~placed && numel(A.f) < A.cap
        A.X = [A.X; x]; A.f = [A.f; f(t)]; A.q = [A.q; q(t)];
    end
end
[A.f, o] = sort(A.f); A.X = A.X(o,:); A.q = A.q(o);
end

function b = archive_best(A)
if isempty(A.f), b = Inf; else, b = A.f(1); end
end

% =========================================================================
% UTILITIES
% =========================================================================
function C = merge(D, U)
C = D;
if nargin < 2 || isempty(U), return; end
fn = fieldnames(U);
for i = 1:numel(fn), C.(fn{i}) = U.(fn{i}); end
end

function rng_set(s)
try
    rng(s, 'twister');
catch
    rand('state', s); randn('state', s); %#ok<RAND>
end
end

function s = tern(c,a,b), if c, s = a; else, s = b; end, end
