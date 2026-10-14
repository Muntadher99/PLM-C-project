function OUT = fssa_sweep(CFG)
%FSSA_SWEEP  Followed-solution simulated annealing (FSSA): reimplementation
%   of the published SHM method with continuation [2], run under the matched
%   evaluation protocol of Section 4.
%
%   OUT = fssa_sweep(CFG)
%
%   [2] L. G. Franquelo, J. Napoles, R. C. Portillo Guisado, J. I. Leon and
%       M. A. Aguirre, IEEE Trans. Ind. Electron., 54(6):3022-3029, 2007.
%
%   The method walks the modulation index in steps of MA_STEP = 0.01 and seeds
%   each step with the solution of the previous one (zeroth-order predictor);
%   each step is corrected by simulated annealing. Every intermediate step is
%   charged to the reported operating point it walks toward, so OUT.NFE is per
%   reported point.
%
%   Settings taken from [2]: 15 angles, MA_STEP = 0.01, T0 = 0.1, ordering and
%   minimum pulse width as hard constraints, ascending sweep in Ma.
%   Settings not given in [2] (Table 2 of the paper): geometric cooling from T0
%   to TMIN_RATIO*T0 sized to the available evaluations, MOVES_PER_T proposals
%   per temperature, Gaussian moves of all angles with standard deviation
%   STEP_DEG*(T/T0) degrees, BOOT_STARTS sorted box draws at the first point.
%
%   CFG fields (defaults in brackets)
%     Ma_list      reported operating points            [1.20:-0.05:0.20]
%     NFE_TARGET   evaluations per reported point        [2000]
%     seed         random seed (one stream per sweep)    [1]
%     N, nmax      angles and highest counted order      [15, 49]
%     minGap_deg   minimum gap between angles, degrees   [0.25]
%     T0           initial temperature, relative to the objective at entry [0.10]
%     MA_STEP      internal continuation step            [0.01]
%     MARGIN       target fraction of the limit, OBJ='grid' only [0.80]
%     COOL         fixed cooling factor; [] = sized to the budget []
%     TMIN_RATIO   final / initial temperature           [1e-6]
%     STEP_DEG     proposal standard deviation at T0     [5.0]
%     NPERT        angles moved per proposal, [] = all   []
%     MOVES_PER_T  evaluated proposals per temperature   [12]
%     BOOT_STARTS  sorted box draws at the first point   [12]
%     OBJ          'ls'   sum b_n^2 + FUND_W*|b1 - Ma|   (main benchmark, Tables 3, 4)
%                  'grid' weighted grid-code objective of [2] (Section 5.6) ['ls']
%     FUND_P       'lin' |b1-Ma| (as in [2]) | 'quad'   ['lin']
%     FUND_W       fundamental penalty weight, OBJ='ls'  [1.0]
%     T_ABSOLUTE   true = T0 used literally               [false]
%     RESTART_ON_FAIL  cold restart when the branch is lost (not in [2]) [false]
%     FAIL_TOL, FAIL_RELDEG, FAIL_ABS  branch-loss tests  [1e-2, 1.5, 0.50]
%     SWEEP_DIR    'ascending' (as in [2]) | 'descending' | 'as-given' ['ascending']
%     X0           optional start at the first point     []
%     verbose      print a one-line summary              [false]
%
%   NFE: one per evaluation of Eq. (2). Proposals rejected as infeasible are
%   not evaluated and cost nothing.
%
%   OUT fields: Ma, theta_deg, NFE (per reported point), n_minima (= 1),
%   branch_lost, restarts.

% ---- configuration ------------------------------------------------------
D.Ma_list        = 1.20:-0.05:0.20;
D.NFE_TARGET     = 2000;
D.seed           = 1;
D.N              = 15;
D.nmax           = 49;
D.minGap_deg     = 0.25;
D.T0             = 0.10;
D.MA_STEP        = 0.01;
D.MARGIN         = 0.80;
D.COOL           = [];
D.TMIN_RATIO     = 1e-6;
D.STEP_DEG       = 5.0;
D.NPERT          = [];
D.MOVES_PER_T    = 12;
D.BOOT_STARTS    = 12;
D.OBJ            = 'ls';
D.FUND_P         = 'lin';
D.FUND_W         = 1.0;
D.T_ABSOLUTE     = false;
D.RESTART_ON_FAIL= false;
D.FAIL_TOL       = 1e-2;
D.FAIL_RELDEG    = 1.5;
D.FAIL_ABS       = 0.50;
D.SWEEP_DIR      = 'ascending';
D.X0             = [];
D.verbose        = false;

if nargin < 1, CFG = struct(); end
fn = fieldnames(CFG); C = D;
for i = 1:numel(fn), C.(fn{i}) = CFG.(fn{i}); end

F = plm_core();  P = F.defaults(C.N, C.nmax, C.minGap_deg);
MAS = C.Ma_list(:);  nM = numel(MAS);
rng_set(C.seed);
F.nfe_reset();

% ---- internal continuation ladder at MA_STEP ----------------------------
switch lower(C.SWEEP_DIR)
  case 'ascending',  [~, ord] = sort(MAS, 'ascend');
  case 'descending', [~, ord] = sort(MAS, 'descend');
  otherwise,         ord = (1:nM).';
end
[ladder, owner] = build_ladder(MAS(ord), C.MA_STEP);
owner = ord(owner);          % index of the reported point each step serves
nL = numel(ladder);
share = zeros(nL,1);
for k = 1:nM
    m = (owner == k);
    share(m) = C.NFE_TARGET / max(sum(m), 1);
end

OUT.Ma        = MAS;
OUT.theta_deg = nan(nM, P.N);
OUT.NFE       = zeros(nM, 1);
OUT.n_minima  = ones(nM, 1);
OUT.branch_lost = false(nM,1);
OUT.restarts  = 0;

x = C.X0(:).';  if ~isempty(x), x = sort(x); end
n_lost = 0;  thd_prev = NaN;
for j = 1:nL
    Ma = ladder(j);  k = owner(j);
    n0 = F.nfe();

    if isempty(x)
        x = bootstrap(Ma, P, F, C, min(share(j), C.NFE_TARGET));
    end
    [x, ok] = anneal(x, Ma, P, F, C, share(j) - (F.nfe() - n0));
    [ok, thd_now] = branch_ok(x, Ma, P, F, C, ok, thd_prev);

    if ~ok
        n_lost = n_lost + 1;
        OUT.branch_lost(k) = true;
        if C.RESTART_ON_FAIL
            % optional cold restart (not part of [2])
            x2 = bootstrap(Ma, P, F, C, 0.4*share(j));
            [x2, ~] = anneal(x2, Ma, P, F, C, share(j) - (F.nfe() - n0));
            R2 = F.report(sort(x2), Ma, P);
            if ~isfinite(thd_now) || R2.THD49 < thd_now
                x = x2;  thd_now = R2.THD49;
            end
            OUT.restarts = OUT.restarts + 1;
        end
    end
    thd_prev = thd_now;

    OUT.NFE(k) = OUT.NFE(k) + (F.nfe() - n0);
    if abs(ladder(j) - MAS(k)) < 1e-12
        OUT.theta_deg(k,:) = rad2deg(sort(x));    % a reported point
    end
end

if C.verbose
    fprintf(['fssa_sweep: %d fine steps over %d reported Ma, %d lost branch, ' ...
             '%d restarts, mean NFE %.0f\n'], nL, nM, n_lost, OUT.restarts, ...
             mean(OUT.NFE));
end
end

% =========================================================================
function [ok, thd] = branch_ok(x, Ma, P, F, C, ok_b1, thd_prev)
%BRANCH_OK  The tracked branch is considered lost if the fundamental was not
%   met (ok_b1 false), THD49 grew by more than FAIL_RELDEG times the previous
%   step, or THD49 exceeds FAIL_ABS. Charges no NFE.
R   = F.report(sort(x), Ma, P);
thd = R.THD49;
ok  = ok_b1;
if isfinite(thd_prev) && thd > C.FAIL_RELDEG * max(thd_prev, 1e-6), ok = false; end
if thd > C.FAIL_ABS, ok = false; end
end

% =========================================================================
function [ladder, owner] = build_ladder(MAS, step)
%BUILD_LADDER  Fine Ma sequence and the reported point each step serves.
ladder = MAS(1);  owner = 1;
for k = 2:numel(MAS)
    ns = max(1, round(abs(MAS(k) - MAS(k-1)) / step));
    for j = 1:ns
        ladder(end+1,1) = MAS(k-1) + j*(MAS(k) - MAS(k-1))/ns; %#ok<AGROW>
        owner(end+1,1)  = k;                                    %#ok<AGROW>
    end
end
end

% =========================================================================
function x = bootstrap(Ma, P, F, C, budget)
%BOOTSTRAP  First point of the sweep: best of BOOT_STARTS uniform box draws,
%   sorted, with draws violating the minimum gap rejected.
best = []; bJ = inf;
n = max(1, min(C.BOOT_STARTS, floor(max(budget,1))));
for i = 1:n
    x = [];
    for t = 1:200
        c = sort(P.lb + (P.ub - P.lb)*rand(1, P.N));
        if all(diff(c) >= P.minGap), x = c; break; end
    end
    if isempty(x), continue; end
    J = objective(x, Ma, P, F, C);
    if J < bJ, bJ = J; best = x; end
end
if isempty(best)
    best = P.lb + P.minGap*(0:P.N-1) + (P.ub-P.lb-P.minGap*(P.N-1))*(0.5/P.N)*(0:P.N-1);
end
x = best;
end

% =========================================================================
function [x, ok] = anneal(x, Ma, P, F, C, budget)
%ANNEAL  One simulated-annealing correction at one Ma, started from x.
%   Geometric cooling from T0 to TMIN_RATIO*T0 over the evaluations
%   available; MOVES_PER_T evaluated (feasible) proposals per temperature.
if budget < 1, ok = true; return; end
n0 = F.nfe();
J  = objective(x, Ma, P, F, C);
Jref = max(abs(J), eps);
T0e = C.T0;  if ~C.T_ABSOLUTE, T0e = C.T0 * Jref; end
T   = T0e;
sd0 = deg2rad(C.STEP_DEG);
bestx = x;  bestJ = J;

cool = C.COOL;
if isempty(cool)
    nlev = max(2, floor(budget / max(C.MOVES_PER_T,1)));
    cool = C.TMIN_RATIO ^ (1/nlev);
end
Tmin = C.TMIN_RATIO * T0e;

while (F.nfe() - n0) < budget && T > Tmin
    sd = sd0 * max(T / max(T0e, eps), 1e-3);
    nev = 0;  tries = 0;  cap = 100 * C.MOVES_PER_T;
    while nev < C.MOVES_PER_T && tries < cap
        tries = tries + 1;
        if (F.nfe() - n0) >= budget, break; end
        y = x;
        np = C.NPERT;  if isempty(np), np = P.N; end
        idx = randperm(P.N, min(np, P.N));
        y(idx) = y(idx) + sd*randn(1, numel(idx));
        if ~feasible(y, P), continue; end          % infeasible: not evaluated
        nev = nev + 1;
        Jy = objective(y, Ma, P, F, C);
        d  = Jy - J;
        if d < 0 || rand < exp(-d / max(T, eps))
            x = y;  J = Jy;
            if J < bestJ, bestJ = J; bestx = x; end
        end
    end
    T = cool * T;
end

x = bestx;
b1 = (4/pi) * (cos(x) * P.s(:));
ok = abs(b1 - Ma) <= C.FAIL_TOL;
end

% =========================================================================
function J = objective(x, Ma, P, F, C)
%OBJECTIVE  FSSA objective; charges one NFE.
B  = F.bn(x, P);            % 1 x M harmonic vector
F.nfe_add(1);
b1 = (4/pi) * (cos(x) * P.s(:));
switch lower(C.OBJ)
  case 'ls'
    % sum of squares with the linear fundamental penalty of [2] (Table 2)
    if strcmpi(C.FUND_P, 'quad')
        J = sum(B.^2) + C.FUND_W * (b1 - Ma)^2;
    else
        J = sum(B.^2) + C.FUND_W * abs(b1 - Ma);
    end
  case 'grid'
    % objective of [2]: weighted limit violations plus a THD term, weights
    % raised to 1000 on violation and 15000 on the fundamental
    r   = abs(B) / max(abs(b1), eps);
    E   = r - C.MARGIN * P.LIM;
    c   = ones(size(E));  c(E > 0) = 1000;
    thd = sqrt(sum(B.^2)) / max(abs(b1), eps);
    cT  = 1;  if thd > 0.08, cT = 1000; end
    J   = sum(c .* abs(E)) + cT*thd + 15000*abs(b1 - Ma);
  otherwise
    error('fssa_sweep:obj', 'CFG.OBJ must be ''ls'' or ''grid''');
end
end

% =========================================================================
function ok = feasible(x, P)
%FEASIBLE  Ordering, minimum gap and bounds tested on the raw proposal.
ok = x(1) >= P.lb && x(end) <= P.ub && all(diff(x) >= P.minGap);
end

function rng_set(s)
try, rng(s, 'twister'); catch, rand('state', s); randn('state', s); end %#ok<RAND>
end
