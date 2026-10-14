function S = corrector_ablation(CFG)
%CORRECTOR_ABLATION  Matched-start ablation of the local corrector (paper,
%   Section 5.4 and Table 6).
%
%   S = corrector_ablation()                               % production caps
%   S = corrector_ablation(struct('iter_scale', 10, ...
%           'variants', {{'proposed','gn','lm1','theta'}}))  % run to convergence
%
%   The same K starts per operating point (drawn once from the feasible set
%   with a fixed seed) are solved by the proposed corrector and by variants
%   that change one element each:
%     'proposed'  two-metric projected LM in gap coordinates (= plm_core solve,
%                 'thd' mode; checked start by start against F.solve)
%     'gn'        damping removed: undamped Gauss-Newton with a 12-step
%                 projected backtracking arc search
%     'lm1'       one metric on all variables and the gain ratio of the raw
%                 step (conventional LM)
%     'theta'     the proposed LM in angle coordinates with the Euclidean
%                 projection onto (4) (gap parameterisation removed)
%     'batch1'    'proposed' solved one start at a time (batching removed)
%   Every variant uses the same weight schedule, the same final Newton
%   restoration onto {b1 = Ma}, the same acceptance rule and the same
%   first-order KKT test.
%
%   CFG fields (defaults in brackets)
%     K           starts per operating point        [200]
%     Ma          operating points                  [0.40 0.60 0.80 1.00 1.20]
%     seed        starts at point i use seed + i    [20260930]
%     N           number of angles                  [15]
%     variants    variants to run                   [{'proposed','gn','lm1','theta','batch1'}]
%     kkt_tol     relative KKT residual tolerance   [1e-3]
%     iter_scale  multiplier of every phase's iteration cap [1]
%     out         output CSV                        ['corrector_ablation.csv']
%
%   Output CSV columns: variant, Ma, accepted_pct, kkt_pct, best_THD49_pct,
%   median_THD49_pct, nfe_per_start, time_s.
%   Table 6 of the paper: accepted, median, best and NFE from the run with
%   iter_scale = 1 (corrector_ablation_production.csv); KKT share from the
%   run with iter_scale = 10 (corrector_ablation_converged.csv); averages
%   over the five points. Both runs: run_ablations.

D.K       = 200;
D.Ma      = [0.40 0.60 0.80 1.00 1.20];
D.seed    = 20260930;
D.N       = 15;
D.variants= {'proposed','gn','lm1','theta','batch1'};
D.kkt_tol = 1e-3;
D.out     = 'corrector_ablation.csv';
D.iter_scale = 1;
if nargin < 1, CFG = struct(); end
f = fieldnames(CFG); for i = 1:numel(f), D.(f{i}) = CFG.(f{i}); end
C = D;

F = plm_core();  P = F.defaults(C.N, 49);
nV = numel(C.variants);  nM = numel(C.Ma);
R = struct();
fprintf('\n== corrector_ablation: K = %d matched starts per operating point, iteration caps x%g ==\n', C.K, C.iter_scale);

rows = {};
for mi = 1:nM
    Ma = C.Ma(mi);
    rng_set(C.seed + mi);
    X0 = F.zsample(rand(C.K, P.N), P);          % common starts
    ref = [];
    for vi = 1:nV
        vn = C.variants{vi};
        F.nfe_reset();  t0 = tic;
        switch vn
          case 'batch1'
            X = zeros(size(X0));
            for k = 1:C.K, X(k,:) = corr_solve(X0(k,:), Ma, P, F, 'proposed', C.iter_scale); end
          otherwise
            X = corr_solve(X0, Ma, P, F, vn, C.iter_scale);
        end
        tt = toc(t0);  nfe = F.nfe();
        M = judge(X, Ma, P, F, C.kkt_tol);
        if strcmp(vn, 'proposed') && C.iter_scale == 1
            F.nfe_reset();
            Xs = F.solve(X0, Ma, P);                % plm_core corrector, same starts
            dmax = max(max(abs(sort(Xs,2) - sort(X,2))));
            fprintf('  Ma=%.2f  check vs plm_core solve: max |d theta| = %.2e rad\n', Ma, dmax);
            ref = X;
        end
        if strcmp(vn, 'batch1') && ~isempty(ref)
            fprintf('  Ma=%.2f  batch1 vs batched: max |d theta| = %.2e rad\n', Ma, ...
                    max(max(abs(ref - X))));
        end
        rows(end+1,:) = {vn, Ma, 100*mean(M.acc), 100*mean(M.kkt), ...
                         100*M.best, 100*M.med, nfe/C.K, tt}; %#ok<AGROW>
        R.(vn).acc(mi) = mean(M.acc); R.(vn).kkt(mi) = mean(M.kkt);
        R.(vn).best(mi) = M.best; R.(vn).med(mi) = M.med;
        R.(vn).nfe(mi) = nfe/C.K;  R.(vn).t(mi) = tt;
        R.(vn).thd{mi} = M.thd;
    end
end

% ---- print ---------------------------------------------------------------
fprintf('\n%-9s %5s %8s %7s %10s %10s %9s %8s\n', 'variant','Ma','accept%','KKT%', ...
        'best THD%','med THD%','NFE/start','time s');
for r = 1:size(rows,1)
    fprintf('%-9s %5.2f %8.1f %7.1f %10.4f %10.4f %9.1f %8.2f\n', rows{r,:});
end
fprintf('\npooled over the %d operating points:\n', nM);
fprintf('%-9s %8s %7s %12s %12s %9s %8s\n', 'variant','accept%','KKT%', ...
        'mean best%','median all%','NFE/start','time s');
for vi = 1:nV
    vn = C.variants{vi};  q = R.(vn);
    allthd = [q.thd{:}];  allthd = allthd(isfinite(allthd));
    fprintf('%-9s %8.1f %7.1f %12.4f %12.4f %9.1f %8.2f\n', vn, 100*mean(q.acc), ...
        100*mean(q.kkt), 100*mean(q.best), 100*median(allthd), mean(q.nfe), sum(q.t));
end

% ---- write ---------------------------------------------------------------
fid = fopen(C.out, 'w');
fprintf(fid, 'variant,Ma,accepted_pct,kkt_pct,best_THD49_pct,median_THD49_pct,nfe_per_start,time_s\n');
for r = 1:size(rows,1)
    fprintf(fid, '%s,%.2f,%.4f,%.4f,%.6f,%.6f,%.2f,%.4f\n', rows{r,:});
end
fclose(fid);
fprintf('\nwrote %s\n', C.out);
S = R;
end

% =========================================================================
% THE CORRECTOR, WITH ONE SWITCH PER VARIANT
% =========================================================================
function TH = corr_solve(TH0, Ma, P, F, variant, iscale)
sched = [3e1 30; 1e3 12; 1e5 8];            % plm_core 'thd' schedule
if nargin >= 6, sched(:,2) = round(sched(:,2) * iscale); end
o.narc = 2;  o.gtol = 1e-13;
TH0 = F.project(TH0, P);
if strcmp(variant, 'theta')
    X = TH0;
    for r = 1:size(sched,1)
        X = lm_theta(X, Ma, P, F, sched(r,1), sched(r,2), o);
    end
    TH = newton_b1(X, Ma, P, F);
    return;
end
V = F.projv(F.th2v(TH0, P), P.span);
if strcmp(variant, 'gn'), o.narc = 12; end
for r = 1:size(sched,1)
    V = lm_gap(V, Ma, P, F, sched(r,1), sched(r,2), o, variant);
end
TH = newton_b1(F.v2th(V, P), Ma, P, F);
end

% -------------------------------------------------------------------------
function V = lm_gap(V, Ma, P, F, wg, iters, o, variant)
%LM_GAP  plm_core lm_phase ('thd' mode, unit residual weights) with switches.
[K, N] = size(V);  I = eye(N);  S = P.span;
twometric = ~strcmp(variant, 'lm1');
damped    = ~strcmp(variant, 'gn');
rawgain   =  strcmp(variant, 'lm1');
lam = 1e-3 * ones(K,1);
if ~damped, lam(:) = 0; end
act = true(K,1);
phi = sum(resid(F.v2th(V,P), Ma, P, F, wg, false).^2, 2);
for it = 1:iters
    idx = find(act);
    if isempty(idx), break; end
    v  = V(idx,:);
    [Rr, Jt] = resid(F.v2th(v,P), Ma, P, F, wg, true);
    Jv = flip(cumsum(flip(Jt, 3), 3), 3);
    [H, b] = gram(Jv, Rr);
    if twometric
        A = (v <= 1e-11) & (b < 0);
    else
        A = false(size(v));
    end
    dgv = max(diag3(H), 1e-14);
    Hm = zero_rc(H, A);
    Hm = Hm + pagediag(A + (~A) .* (lam(idx) .* dgv), I);
    bm = b;  bm(A) = 0;
    e1 = ones(size(v)); e1(A) = 0;
    sol = solve2(Hm, bm, e1);
    d = sol(:,:,1); w = sol(:,:,2);
    d(A) = 0; w(A) = 0;
    hitf = (sum(v,2) >= S - 1e-11) & (sum(d,2) > 0);
    den  = sum(w,2);  den(abs(den) < 1e-300) = 1e-300;
    d = d - (hitf .* (sum(d,2)./den)) .* w;
    alpha = ones(numel(idx),1);
    nv = v; nphi = phi(idx); gain = zeros(numel(idx),1);
    pend = true(numel(idx),1);
    for a = 1:o.narc
        sub = find(pend);
        if isempty(sub), break; end
        tv = F.projv(v(sub,:) + alpha(sub).*d(sub,:), S);
        ph = sum(resid(F.v2th(tv,P), Ma, P, F, wg, false).^2, 2);
        if rawgain, dp = alpha(sub).*d(sub,:); else, dp = tv - v(sub,:); end
        pr = 2*sum(b(sub,:).*dp, 2) - quad3(H(:,:,sub), dp);
        pr(abs(pr) < 1e-300) = 1e-300;
        gg = (phi(idx(sub)) - ph) ./ pr;
        hh = ph < phi(idx(sub)) - 1e-20;
        s = sub(hh);
        nv(s,:) = tv(hh,:); nphi(s) = ph(hh); gain(s) = gg(hh);
        pend(s) = false;
        alpha(sub(~hh)) = alpha(sub(~hh)) * 0.2;
    end
    good = ~pend;  gi = idx(good);
    V(gi,:) = nv(good,:);  phi(gi) = nphi(good);
    step = max(abs(nv - v), [], 2);
    if damped
        up = good & (gain > 0.75) & (alpha > 0.9);
        hold_ = good & (gain > 0.10) & ~up;
        lam(idx(up))  = max(lam(idx(up))/6, 1e-13);
        bad = ~(up | hold_);
        lam(idx(bad)) = min(lam(idx(bad))*7, 1e12);
        act(idx((good & step < o.gtol) | (~good & lam(idx) > 1e9))) = false;
    else
        % undamped: stop when the step vanishes or the arc search finds no decrease
        act(idx((good & step < o.gtol) | ~good)) = false;
    end
end
end

% -------------------------------------------------------------------------
function X = lm_theta(X, Ma, P, F, wg, iters, o)
%LM_THETA  The same LM in angle coordinates; feasibility by the Euclidean
%   projection onto (4) (F.project), one metric.
[K, N] = size(X);  I = eye(N);
lam = 1e-3 * ones(K,1);  act = true(K,1);
phi = sum(resid(X, Ma, P, F, wg, false).^2, 2);
for it = 1:iters
    idx = find(act);
    if isempty(idx), break; end
    x = X(idx,:);
    [Rr, Jt] = resid(x, Ma, P, F, wg, true);
    [H, b] = gram(Jt, Rr);
    dgv = max(diag3(H), 1e-14);
    Hm = H + pagediag(lam(idx) .* dgv, I);
    sol = solve2(Hm, b, ones(size(x)));
    d = sol(:,:,1);
    alpha = ones(numel(idx),1);
    nx = x; nphi = phi(idx); gain = zeros(numel(idx),1);
    pend = true(numel(idx),1);
    for a = 1:o.narc
        sub = find(pend);
        if isempty(sub), break; end
        tx = F.project(x(sub,:) + alpha(sub).*d(sub,:), P);
        ph = sum(resid(tx, Ma, P, F, wg, false).^2, 2);
        dp = tx - x(sub,:);
        pr = 2*sum(b(sub,:).*dp, 2) - quad3(H(:,:,sub), dp);
        pr(abs(pr) < 1e-300) = 1e-300;
        gg = (phi(idx(sub)) - ph) ./ pr;
        hh = ph < phi(idx(sub)) - 1e-20;
        s = sub(hh);
        nx(s,:) = tx(hh,:); nphi(s) = ph(hh); gain(s) = gg(hh);
        pend(s) = false;
        alpha(sub(~hh)) = alpha(sub(~hh)) * 0.2;
    end
    good = ~pend;  gi = idx(good);
    X(gi,:) = nx(good,:);  phi(gi) = nphi(good);
    up = good & (gain > 0.75) & (alpha > 0.9);
    hold_ = good & (gain > 0.10) & ~up;
    lam(idx(up))  = max(lam(idx(up))/6, 1e-13);
    bad = ~(up | hold_);
    lam(idx(bad)) = min(lam(idx(bad))*7, 1e12);
    step = max(abs(nx - x), [], 2);
    act(idx((good & step < o.gtol) | (~good & lam(idx) > 1e9))) = false;
end
end

% -------------------------------------------------------------------------
function [R, J] = resid(TH, Ma, P, F, wg, needJ)
%RESID  Residual of objective (8), 'thd' mode; charges K NFE.
K = size(TH,1);
F.nfe_add(K);
if needJ, [B, JB, G1] = F.bnjac(TH, P); else, B = F.bn(TH, P); end
g = (4/pi)*(cos(TH)*P.s(:)) - Ma;
R = [B, wg * g];
if ~needJ, J = []; return; end
J = cat(2, JB, reshape(wg*G1, K, 1, []));
end

function TH = newton_b1(TH, Ma, P, F)
%NEWTON_B1  Newton restoration onto {b1 = Ma}, as in plm_core.
for it = 1:5
    g = (4/pi)*(cos(TH)*P.s(:)) - Ma;
    F.nfe_add(size(TH,1));
    if max(abs(g)) < 1e-15, break; end
    G1 = -(4/pi) * sin(TH) .* P.s;
    nr = max(sum(G1.^2, 2), 1e-300);
    TH = F.project(TH - (g./nr).*G1, P);
end
end

% =========================================================================
% ACCEPTANCE, STATIONARITY AND DISTORTION (identical for every variant)
% =========================================================================
function M = judge(X, Ma, P, F, ktol)
Rp = F.report(X, Ma, P);
M.acc = (Rp.b1err < P.FEAS_B1) & (Rp.viol < 1e-9) & isfinite(Rp.THD49);
M.thd = Rp.THD49.';  M.thd(~M.acc) = Inf;
M.kkt = false(size(M.acc));
for k = find(M.acc(:)).'
    M.kkt(k) = kkt_rel(sort(X(k,:)), P, F) < ktol;
end
f = M.thd(isfinite(M.thd));
if isempty(f), M.best = NaN; M.med = NaN; else, M.best = min(f); M.med = median(f); end
end

function r = kkt_rel(th, P, F)
%KKT_REL  Relative first-order KKT residual of min 1/2||b||^2 s.t. b1 = Ma
%   and (4); the active set is refined by multiplier sign.
N = P.N;
[B, JB, G1] = F.bnjac(th, P);
J = reshape(JB, P.M, N);  g = J.' * B(:);
tolA = 1e-9;
Aall = zeros(N+1, N); slk = zeros(N+1,1);
Aall(1,1) = 1; slk(1) = th(1) - P.lb;
for i = 1:N-1, Aall(i+1,i) = -1; Aall(i+1,i+1) = 1; slk(i+1) = th(i+1)-th(i)-P.minGap; end
Aall(N+1,N) = -1; slk(N+1) = P.ub - th(N);
act = find(slk < tolA);
for pass = 1:N+1
    Cm = [G1(:).'; Aall(act,:)];
    ml = -(Cm.') \ g;
    lam = -ml(2:end);
    if isempty(lam) || all(lam >= -1e-10*max(1,norm(g))), break; end
    [~, j] = min(lam); act(j) = [];
end
r = norm(g + Cm.'*ml) / max(norm(g), 1e-300);
end

% =========================================================================
% BATCHED LINEAR ALGEBRA (as in plm_core)
% =========================================================================
function [H, b] = gram(J, R)
[K, M, N] = size(J);
Jp = permute(J, [2 3 1]);
if exist('pagemtimes', 'builtin') == 5 || exist('pagemtimes','file') == 2
    H = pagemtimes(Jp, 'transpose', Jp, 'none');
    b = -reshape(pagemtimes(Jp, 'transpose', reshape(R.', M, 1, K), 'none'), N, K).';
else
    H = zeros(N, N, K);  b = zeros(K, N);
    for k = 1:K
        Jk = Jp(:,:,k);
        H(:,:,k) = Jk.' * Jk;
        b(k,:)   = -(Jk.' * R(k,:).').';
    end
end
end
function d = diag3(H)
N = size(H,1); K = size(H,3); d = zeros(K, N);
for i = 1:N, d(:,i) = reshape(H(i,i,:), K, 1); end
end
function H = zero_rc(H, A)
Am = permute(A, [2 3 1]);
keep = ~Am & ~permute(Am, [2 1 3]);
H = H .* keep;
end
function D = pagediag(d, I)
D = permute(d, [2 3 1]) .* I;
end
function q = quad3(H, d)
K = size(H,3); q = zeros(K,1); dp = permute(d, [2 3 1]);
if exist('pagemtimes','builtin') == 5 || exist('pagemtimes','file') == 2
    q = reshape(pagemtimes(dp, 'transpose', pagemtimes(H, dp), 'none'), K, 1);
else
    for k = 1:K, q(k) = d(k,:) * H(:,:,k) * d(k,:).'; end
end
end
function X = solve2(H, b1v, b2v)
[K, N] = size(b1v);
RHS = cat(2, reshape(b1v.', N, 1, K), reshape(b2v.', N, 1, K));
if exist('pagemldivide','builtin') == 5 || exist('pagemldivide','file') == 2
    S = pagemldivide(H, RHS);
else
    S = zeros(N, 2, K);
    for k = 1:K
        Hk = H(:,:,k);
        [Lc, p] = chol(Hk, 'lower');
        if p == 0
            S(:,:,k) = Lc.' \ (Lc \ RHS(:,:,k));
        else
            S(:,:,k) = (Hk + 1e-10*eye(N)) \ RHS(:,:,k);
        end
    end
end
X = permute(S, [3 1 2]);
end
function rng_set(s)
try, rng(s, 'twister'); catch, rand('state', s); randn('state', s); end %#ok<RAND>
end
