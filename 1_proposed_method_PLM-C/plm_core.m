function F = plm_core()
%PLM_CORE  Harmonic model, feasible-set maps and the batched projected
%   Levenberg-Marquardt corrector used by PLM-C and its ablation arms.
%
%   F = plm_core() returns a struct of function handles. Unless stated
%   otherwise every routine is batched: TH is a K-by-N matrix whose rows are
%   quarter-wave switching-angle sets in radians.
%
%   Problem (paper, Section 2)
%       b_n = 4/(n*pi) * sum_i (-1)^(i+1) * cos(n*theta_i),  n odd      Eq. (2)
%       THD49 = ||(b_5, b_7, ..., b_49)||_2 / |b_1|,  with b_1 = Ma      Eq. (3)
%       theta_1 >= theta_min,  theta_{i+1} - theta_i >= delta,
%       theta_N <= theta_max                                             Eq. (4)
%
%   Problem data
%       P = F.defaults(N, nmax, minGap_deg)
%                                      constants of Section 2.1 (defaults
%                                      N = 15, nmax = 49, minimum gap 0.25 deg)
%       F.orders(nmax)                 odd non-triplen orders 5..nmax
%       F.limits(orders)               grid-code limits L_n of Table 1
%   Harmonic model
%       F.b1(TH, P), F.bn(TH, P)       fundamental and counted harmonics, Eq. (2)
%       [B, JB, G1] = F.bnjac(TH, P)   harmonics with analytic Jacobian, Eq. (9)
%       R = F.report(TH, Ma, P)        THD49, b1 error, violation, q_max, compliance
%       S = F.spectrum(TH, P, nmax)    extended-band diagnostics (charges no NFE)
%   Feasible set (Section 3.1)
%       F.project(TH, P)               Euclidean projection onto (4)
%       F.th2v(TH, P), F.v2th(V, P)    angle <-> gap coordinates, Eq. (6)
%       F.projv(V, S)                  Euclidean projection onto the simplex (7)
%       F.zsample(Z, P)                uniform map (0,1)^N -> feasible set
%       F.theta2z(TH, P)               inverse of F.zsample
%   Corrector (Sections 3.2 and 3.4)
%       [TH, I] = F.solve(TH0, Ma, P, opts)
%   Evaluation counter: one NFE = one evaluation of Eq. (2) for one candidate,
%   Jacobian included (Section 4)
%       F.nfe(), F.nfe_reset(), F.nfe_add(k)
%   Self test
%       F.selftest()
%
%   Requires no toolbox; runs in MATLAB R2024b and GNU Octave.

    F.defaults   = @defaults;
    F.b1         = @b1_of;
    F.bn         = @bn_of;
    F.bnjac      = @bnjac;
    F.report     = @report;
    F.spectrum   = @spectrum;
    F.project    = @project;
    F.th2v       = @th2v;
    F.v2th       = @v2th;
    F.projv      = @proj_simplex_le;
    F.solve      = @solve;
    F.zsample    = @z2theta;
    F.theta2z    = @theta2z;
    F.orders     = @odd_nontriplen;
    F.limits     = @grid_limits;
    F.nfe        = @nfe_get;
    F.nfe_reset  = @nfe_reset;
    F.nfe_add    = @nfe_bump;
    F.selftest   = @selftest;
end

% =========================================================================
% PROBLEM DATA
% =========================================================================
function P = defaults(N, nmax, minGap_deg)
%DEFAULTS  Problem constants of Section 2.1 and derived quantities.
%   minGap_deg is the minimum gap delta in degrees (0.25 in Section 5;
%   0.50 in the minimum-gap study of Section 6.1).
    if nargin < 1 || isempty(N),    N = 15;    end
    if nargin < 2 || isempty(nmax), nmax = 49; end
    if nargin < 3 || isempty(minGap_deg), minGap_deg = 0.25; end
    P.N       = N;
    P.minGap  = deg2rad(minGap_deg);   % delta
    P.lb      = deg2rad(0.5);          % theta_min
    P.ub      = deg2rad(89.0);         % theta_max
    P.ORD     = odd_nontriplen(nmax);
    P.LIM     = grid_limits(P.ORD);
    P.W       = harm_weights(P.ORD);
    P.WN      = P.W ./ max(P.W);
    P.s       = (-1).^(0:N-1);         % s_i = (-1)^(i+1)
    P.shift   = P.minGap * (0:N-1);
    P.hi      = P.ub - P.minGap*(N-1);
    P.span    = P.ub - P.lb - (N-1)*P.minGap;   % S in Eq. (7)
    P.COMP_TOL= 1e-4;                  % compliance tolerance on q_max (Section 3.4)
    P.FEAS_B1 = 1e-8;                  % acceptance threshold on |b1 - Ma| (Section 4)
    P.M       = numel(P.ORD);
end

function o = odd_nontriplen(nMax)
%ODD_NONTRIPLEN  Counted harmonic orders 5, 7, 11, 13, ..., nMax.
    a = 5:2:nMax;  o = a(mod(a,3) ~= 0);
end

function L = grid_limits(orders)
%GRID_LIMITS  Grid-code limit L_n (fraction of the fundamental), Table 1:
%   EN 50160 for n <= 25, CIGRE WG 36-05 rule (0.2 + 32.5/n)% above.
    L = zeros(size(orders));
    for i = 1:numel(orders)
        n = orders(i);
        switch n
            case 5,  L(i)=0.06;  case 7,  L(i)=0.05;
            case 11, L(i)=0.035; case 13, L(i)=0.03;
            case 17, L(i)=0.02;  case 19, L(i)=0.015;
            case 23, L(i)=0.015; case 25, L(i)=0.015;
            otherwise, L(i) = (0.2 + 32.5/n)/100;
        end
    end
end

function w = harm_weights(orders)
%HARM_WEIGHTS  Relative order weights (stored in P.W / P.WN; not used by solve).
    w = zeros(size(orders));
    for i = 1:numel(orders)
        n = orders(i);
        if     n <= 7,  w(i) = 800;
        elseif n <= 13, w(i) = 600;
        elseif n <= 25, w(i) = 350;
        else,           w(i) = 180;
        end
    end
end

% =========================================================================
% HARMONIC MODEL  (batched: TH is K x N)
% =========================================================================
function v = b1_of(TH, P)
%B1_OF  Fundamental amplitude b_1 of every row, K x 1.
    v = (4/pi) * (cos(TH) * P.s(:));
end

function B = bn_of(TH, P)
%BN_OF  Counted harmonic amplitudes b_n, n in P.ORD, K x M.
    [K, N] = size(TH);  M = P.M;
    A = reshape(P.ORD, 1, M, 1) .* reshape(TH, K, 1, N);
    B = reshape(reshape(cos(A), [], N) * P.s(:), K, M) .* (4./(P.ORD*pi));
end

function [B, JB, G1] = bnjac(TH, P)
%BNJAC  Harmonics and analytic Jacobians, Eq. (9).
%   B  (K x M)      b_n
%   JB (K x M x N)  d b_n / d theta_i = -(4/pi) * s_i * sin(n*theta_i)
%   G1 (K x N)      d b_1 / d theta_i
    [K, N] = size(TH);  M = P.M;
    A  = reshape(P.ORD, 1, M, 1) .* reshape(TH, K, 1, N);
    CA = cos(A);  SA = sin(A);
    B  = reshape(reshape(CA, [], N) * P.s(:), K, M) .* (4./(P.ORD*pi));
    JB = -(4/pi) * SA .* reshape(P.s, 1, 1, N);
    G1 = -(4/pi) * sin(TH) .* P.s;
end

function R = report(TH, Ma, P)
%REPORT  Reported quantities of every row (rows are sorted first).
%   R.THD49, R.b1, R.b1err = |b1 - Ma|, R.viol (worst violation of (4)),
%   R.MAXQ = max_n q_n with q_n = |b_n| / (L_n |b_1|), R.nviol (orders with
%   q_n > 1), R.ratio = |b_n|/|b_1|, R.compliant = (MAXQ <= 1 + COMP_TOL).
    TH = sort(TH, 2);
    B1 = b1_of(TH, P);  B = bn_of(TH, P);
    ratio = abs(B) ./ abs(B1);
    q     = ratio ./ P.LIM;
    gap_viol   = max([zeros(size(TH,1),1), P.minGap - diff(TH,1,2)], [], 2);
    bound_viol = max(max(0, P.lb - TH(:,1)), max(0, TH(:,end) - P.ub));
    R.THD49 = sqrt(sum(B.^2, 2)) ./ abs(B1);
    R.b1    = B1;
    R.b1err = abs(B1 - Ma);
    R.viol  = max(gap_viol, bound_viol);
    R.MAXQ  = max(q, [], 2);
    R.nviol = sum(q > 1, 2);
    R.ratio = ratio;
    R.compliant = R.MAXQ <= 1 + P.COMP_TOL;
end

function S = spectrum(TH, P, nmax)
%SPECTRUM  Extended-band diagnostics of final solutions.
%   S = SPECTRUM(TH, P)         uses nmax = 999
%   S = SPECTRUM(TH, P, nmax)
%
%   All bands are odd and non-triplen (triplen orders cancel in the
%   line-to-line voltage); the third harmonic is reported separately as
%   b3_ratio. WTHD = sqrt(sum (b_n/n)^2) / |b_1|. Charges no NFE: it is
%   evaluated only on returned solutions, never inside the search.
    if nargin < 3 || isempty(nmax), nmax = 999; end
    TH  = sort(TH, 2);
    a1  = abs(b1_of(TH, P));
    o   = odd_nontriplen(nmax);
    B   = bn_ord(TH, o, P.s);
    r   = abs(B) ./ a1;
    w   = 1 ./ o;
    i49 = o <= 49;   i99 = o <= 99;
    S.THD49    = sqrt(sum(B(:,i49).^2, 2)) ./ a1;
    S.THD99    = sqrt(sum(B(:,i99).^2, 2)) ./ a1;
    S.THDfull  = sqrt(sum(B.^2, 2))        ./ a1;
    S.WTHD49   = sqrt(sum((B(:,i49).*w(i49)).^2, 2)) ./ a1;
    S.WTHDfull = sqrt(sum((B.*w).^2, 2))             ./ a1;
    oob = ~i49;
    if any(oob)
        [mx, j] = max(r(:,oob), [], 2);
        oo = o(oob);
        S.OOB_max = mx;  S.OOB_order = reshape(oo(j), [], 1);
    else
        S.OOB_max = zeros(size(a1));  S.OOB_order = nan(size(a1));
    end
    S.b3_ratio = abs(bn_ord(TH, 3, P.s)) ./ a1;
    S.orders = o;  S.ratio = r;  S.nmax = nmax;
end

function B = bn_ord(TH, ord, s)
%BN_ORD  b_n for an arbitrary list of orders.
    [K, N] = size(TH);  ord = ord(:).';  M = numel(ord);
    A = reshape(ord, 1, M, 1) .* reshape(TH, K, 1, N);
    B = reshape(reshape(cos(A), [], N) * s(:), K, M) .* (4./(ord*pi));
end

% =========================================================================
% FEASIBLE SET  (Section 3.1)
% =========================================================================
function A = project(TH, P)
%PROJECT  Euclidean projection of every row onto the ordered set (4).
%   Isotonic regression of TH - shift via the min-max formula
%       u_i = max_{j<=i} min_{k>=i} mean(y_j..y_k),
%   O(N^2) per row and vectorised over the K rows.
    [K, N] = size(TH);
    Y  = TH - P.shift;
    csum = [zeros(K,1), cumsum(Y, 2)];
    j_idx = reshape(0:N-1, 1, N, 1);
    i_idx = reshape(1:N,  1, 1, N);
    len = i_idx - j_idx;                                  % 1 x N x N
    block_mean = (reshape(csum(:,2:end), K, 1, N) - reshape(csum(:,1:end-1), K, N, 1)) ...
         ./ max(len, 1);
    block_mean(repmat(len <= 0, [K 1 1])) = Inf;
    min_suffix = flip(cummin(flip(block_mean, 3), 3), 3);   % min over k >= i
    max_prefix = cummax(min_suffix, 2);                     % max over j <= i
    U = zeros(K, N);
    for i = 1:N, U(:,i) = max_prefix(:,i,i); end
    U = min(max(U, P.lb), P.hi);
    A = U + P.shift;
end

function V = th2v(TH, P)
%TH2V  Angles -> gap coordinates, Eq. (6).
    V = [TH(:,1) - P.lb, diff(TH, 1, 2) - P.minGap];
end

function TH = v2th(V, P)
%V2TH  Gap coordinates -> angles (inverse of Eq. (6)).
    TH = P.lb + P.shift + cumsum(V, 2);
end

function W = proj_simplex_le(V, S)
%PROJ_SIMPLEX_LE  Euclidean projection of every row onto {v >= 0, sum(v) <= S},
%   the simplex (7): v -> max(v - tau, 0), tau from one sort [22].
    W = max(V, 0);
    over = sum(W, 2) > S;
    if ~any(over), return; end
    X = V(over, :);  [m, n] = size(X);
    Xs   = sort(X, 2, 'descend');
    css  = cumsum(Xs, 2) - S;
    cond = (Xs - css ./ (1:n)) > 0;
    rho = ones(m, 1);
    for k = 1:n, rho(cond(:,k)) = k; end      % last index where cond holds
    idx = sub2ind([m n], (1:m).', rho);
    tau = css(idx) ./ rho;
    W(over, :) = max(X - tau, 0);
end

% =========================================================================
% SAMPLING  (Section 3.1: starts are drawn directly from the feasible set)
% =========================================================================
function TH = z2theta(Z, P)
%Z2THETA  Map Z in (0,1)^N (K x N) to feasible ordered angles.
%   Order-statistic (stick-breaking) map: a uniform Z gives a uniform point
%   of the feasible set, so random starts need no rejection or repair.
    Z = min(max(Z, 1e-12), 1-1e-12);
    [K, N] = size(Z);
    U = zeros(K, N);  prev = zeros(K, 1);
    for i = 1:N
        t = 1 - (1 - Z(:,i)).^(1/(N-i+1));
        prev = prev + (1 - prev).*t;
        U(:,i) = prev;
    end
    TH = P.lb + P.shift + P.span * U;
end

function Z = theta2z(TH, P)
%THETA2Z  Inverse of Z2THETA.
    [K, N] = size(TH);
    U = min(max((TH - P.lb - P.shift)/P.span, 0), 1-1e-12);
    Z = zeros(K, N);  prev = zeros(K, 1);
    for i = 1:N
        remain = max(1e-15, 1 - prev);
        Z(:,i) = 1 - (1 - (U(:,i) - prev)./remain).^(N-i+1);
        prev = U(:,i);
    end
    Z = min(max(Z, 0), 1);
end

% =========================================================================
% EVALUATION COUNTER
% =========================================================================
function n = nfe_get()
%NFE_GET  Evaluations charged since the last reset.
    n = nfe_bump(0);
end

function nfe_reset()
%NFE_RESET  Reset the evaluation counter to zero.
    nfe_bump(-inf);
end

function n = nfe_bump(k)
%NFE_BUMP  Charge k evaluations (k = -inf resets). External engines such as
%   fmincon must charge their work through this function (F.nfe_add). The
%   count is mirrored in the global NFE_COUNT used by baseline_common.
    persistent c
    if isempty(c), c = 0; end
    if isinf(k) && k < 0, c = 0; else, c = c + k; end
    n = c;
    global NFE_COUNT %#ok<GVMIS>
    NFE_COUNT = c;
end

% =========================================================================
% BATCHED CORRECTOR  (Sections 3.2 and 3.4)
% =========================================================================
function [TH, I] = solve(TH0, Ma, P, opts)
%SOLVE  Batched two-metric projected Levenberg-Marquardt in gap coordinates.
%
%   [TH, I] = SOLVE(TH0, Ma, P)
%   [TH, I] = SOLVE(TH0, Ma, P, opts)
%
%   TH0 is K x N (any points; they are projected onto (4) first). All K rows
%   are advanced together. Each solve ends with a minimum-norm Newton
%   correction onto {b_1 = Ma} followed by a projection.
%
%   opts fields (defaults in brackets)
%     mode      'thd'     min sum b_n^2,          b_1 by weighted residual  ['thd']
%               'minimax' min max_n q_n by Lawson reweighting (Section 3.4)
%               'comply'  min sum b_n^2 + rho*sum max(0, q_n - 1)^2
%     schedule  [w_g, iterations] per phase of the fundamental weight
%               [30 30; 1e3 12; 1e5 8]
%     narc      projected arc-search steps                               [2]
%     gtol      step-size stopping tolerance                             [1e-13]
%     rho       penalty weights of 'comply' mode                         [1e3 3e4]
%     lawson    Lawson reweightings of 'minimax' mode                    [4]
%
%   I fields: ok (accepted: |b1-Ma| < FEAS_B1 and violation < 1e-9), thd,
%   maxq (Inf where not ok), b1err, nviol, compliant.
    if nargin < 4, opts = struct(); end
    opts = set_default(opts, 'mode',     'thd');
    opts = set_default(opts, 'schedule', [3e1 30; 1e3 12; 1e5 8]);
    opts = set_default(opts, 'narc',     2);
    opts = set_default(opts, 'gtol',     1e-13);
    opts = set_default(opts, 'rho',      [1e3 3e4]);
    opts = set_default(opts, 'lawson',   4);

    TH0 = project(TH0, P);
    V   = proj_simplex_le(th2v(TH0, P), P.span);
    [K, ~] = size(V);  M = P.M;
    w_unit = ones(K, M);

    switch lower(opts.mode)
      case 'thd'
        for ph = 1:size(opts.schedule,1)
            V = lm_phase(V, Ma, P, opts.schedule(ph,1), w_unit, 0, opts.schedule(ph,2), opts);
        end
      case 'minimax'
        % Lawson iteratively reweighted least squares, per row:
        % w <- w.*q ./ sum(w.*q)  (Eq. (10), Section 3.4).
        w_lawson = ones(K, M) / M;
        q_scale  = 1 ./ (P.LIM * abs(Ma));
        for it_lawson = 1:opts.lawson
            w_res = sqrt(max(w_lawson, 1e-10)) .* q_scale;
            if it_lawson == 1, sched = opts.schedule; else, sched = [1e3 14; 1e5 8]; end
            for ph = 1:size(sched,1)
                V = lm_phase(V, Ma, P, sched(ph,1), w_res, 0, sched(ph,2), opts);
            end
            q  = abs(bn_of(v2th(V,P), P)) .* q_scale;
            wq = w_lawson .* q;
            w_lawson = wq ./ max(sum(wq,2), 1e-300);
        end
      case 'comply'
        for ir = 1:numel(opts.rho)
            if ir == 1, sched = [3e1 25; 1e3 12; 1e5 8]; else, sched = [1e3 14; 1e5 8]; end
            for ph = 1:size(sched,1)
                V = lm_phase(V, Ma, P, sched(ph,1), w_unit, opts.rho(ir), sched(ph,2), opts);
            end
        end
      otherwise
        error('plm_core:mode', 'mode must be thd | minimax | comply');
    end

    TH = newton_b1(v2th(V, P), Ma, P);
    R  = report(TH, Ma, P);
    I.ok    = (R.b1err < P.FEAS_B1) & (R.viol < 1e-9) & isfinite(R.THD49);
    I.thd   = R.THD49;  I.thd(~I.ok) = Inf;
    I.maxq  = R.MAXQ;   I.maxq(~I.ok) = Inf;
    I.b1err = R.b1err;  I.nviol = R.nviol;
    I.compliant = R.compliant & I.ok;
end

function s = set_default(s, field, value)
%SET_DEFAULT  Set s.(field) = value when the field is absent or empty.
    if ~isfield(s, field) || isempty(s.(field)), s.(field) = value; end
end

% -------------------------------------------------------------------------
function V = lm_phase(V, Ma, P, wg, w_res, rho, iters, opts)
%LM_PHASE  One projected-LM phase at a fixed fundamental weight wg.
%   Two-metric step: Newton metric on the free gaps, identity on the gaps
%   pinned at zero [24, 25]. The gain ratio is evaluated on the PROJECTED
%   step proj(v + d) - v, followed by a short projected arc search
%   (Section 3.2). Damping: lambda/6 if gain > 0.75, 7*lambda if gain < 0.10.
    [K, N] = size(V);
    I_N  = eye(N);
    span = P.span;
    lambda = 1e-3 * ones(K,1);
    active = true(K,1);
    phi = sum(resid(v2th(V,P), Ma, P, wg, w_res, rho, false).^2, 2);

    for it = 1:iters
        rows = find(active);
        if isempty(rows), break; end
        v  = V(rows,:);
        [R, J_theta] = resid(v2th(v,P), Ma, P, wg, w_res(rows,:), rho, true);
        J_gap = flip(cumsum(flip(J_theta, 3), 3), 3);       % chain rule through Eq. (6)
        [H, b] = gram(J_gap, R);

        pinned = (v <= 1e-11) & (b < 0);                    % gaps held at zero
        hdiag  = max(diag3(H), 1e-14);
        H_mod  = H;
        mask   = pinned;
        H_mod  = zero_rc(H_mod, mask);
        H_mod  = H_mod + pagediag(pinned + (~pinned) .* (lambda(rows) .* hdiag), I_N);
        b_mod  = b;  b_mod(pinned) = 0;
        e_free = ones(size(v)); e_free(pinned) = 0;
        sol = solve2(H_mod, b_mod, e_free);
        d = sol(:,:,1); w = sol(:,:,2);
        d(pinned) = 0; w(pinned) = 0;
        % keep the face sum(v) = span (theta_N = theta_max) if the step leaves it
        on_face = (sum(v,2) >= span - 1e-11) & (sum(d,2) > 0);
        den  = sum(w,2);  den(abs(den) < 1e-300) = 1e-300;
        d = d - (on_face .* (sum(d,2)./den)) .* w;

        alpha = ones(numel(rows),1);
        v_acc = v; phi_acc = phi(rows); gain = zeros(numel(rows),1);
        pending = true(numel(rows),1);
        for k_arc = 1:opts.narc
            trial = find(pending);
            if isempty(trial), break; end
            v_trial   = proj_simplex_le(v(trial,:) + alpha(trial).*d(trial,:), span);
            phi_trial = sum(resid(v2th(v_trial,P), Ma, P, wg, w_res(rows(trial),:), rho, false).^2, 2);
            dv   = v_trial - v(trial,:);                    % projected step
            pred = 2*sum(b(trial,:).*dv, 2) - quad3(H(:,:,trial), dv);
            pred(abs(pred) < 1e-300) = 1e-300;
            ratio = (phi(rows(trial)) - phi_trial) ./ pred;
            decreased = phi_trial < phi(rows(trial)) - 1e-20;
            acc = trial(decreased);
            v_acc(acc,:) = v_trial(decreased,:); phi_acc(acc) = phi_trial(decreased); gain(acc) = ratio(decreased);
            pending(acc) = false;
            alpha(trial(~decreased)) = alpha(trial(~decreased)) * 0.2;
        end
        moved = ~pending;
        rows_moved = rows(moved);
        V(rows_moved,:) = v_acc(moved,:);  phi(rows_moved) = phi_acc(moved);
        very_good  = moved & (gain > 0.75) & (alpha > 0.9);
        acceptable = moved & (gain > 0.10) & ~very_good;
        lambda(rows(very_good)) = max(lambda(rows(very_good))/6, 1e-13);
        poor = ~(very_good | acceptable);
        lambda(rows(poor)) = min(lambda(rows(poor))*7, 1e12);
        step = max(abs(v_acc - v), [], 2);
        active(rows((moved & step < opts.gtol) | (~moved & lambda(rows) > 1e9))) = false;
    end
end

% -------------------------------------------------------------------------
function [R, J] = resid(TH, Ma, P, wg, w_res, rho, needJ)
%RESID  Residual vector (and Jacobian in angle coordinates) of objective (8),
%   optionally with the compliance penalty of 'comply' mode. Charges K NFE.
    K = size(TH,1);
    nfe_bump(K);
    if needJ
        [B, JB, G1] = bnjac(TH, P);
    else
        B = bn_of(TH, P);
    end
    g = (4/pi)*(cos(TH)*P.s(:)) - Ma;
    R = [w_res .* B, wg * g];
    if rho > 0
        q_scale = 1 ./ (P.LIM * abs(Ma));
        q  = abs(B) .* q_scale;
        R  = [R, sqrt(rho) * max(0, q - 1)];
    end
    if ~needJ, J = []; return; end
    J = cat(2, w_res .* JB, reshape(wg*G1, K, 1, []));
    if rho > 0
        dq = ((q > 1) .* sign(B) .* q_scale);
        J = cat(2, J, sqrt(rho) * dq .* JB);
    end
end

function TH = newton_b1(TH, Ma, P)
%NEWTON_B1  Minimum-norm Newton correction onto {b_1 = Ma}, five steps,
%   each followed by the projection onto (4). Charges K NFE per step.
    for it = 1:5
        g = (4/pi)*(cos(TH)*P.s(:)) - Ma;
        nfe_bump(size(TH,1));
        if max(abs(g)) < 1e-15, break; end
        G1 = -(4/pi) * sin(TH) .* P.s;
        nr = max(sum(G1.^2, 2), 1e-300);
        TH = project(TH - (g./nr).*G1, P);
    end
end

% =========================================================================
% BATCHED LINEAR ALGEBRA  (uses pagemtimes / pagemldivide when available)
% =========================================================================
function [H, b] = gram(J, R)
%GRAM  H(:,:,k) = J_k' * J_k  and  b(k,:) = -(J_k' * R_k)'   (J is K x M x N).
    [K, M, N] = size(J);
    Jp = permute(J, [2 3 1]);                  % M x N x K
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
%DIAG3  Diagonals of the pages of H, K x N.
    N = size(H,1); K = size(H,3);
    d = zeros(K, N);
    for i = 1:N, d(:,i) = reshape(H(i,i,:), K, 1); end
end

function H = zero_rc(H, A)
%ZERO_RC  Zero the rows and columns of H(:,:,k) flagged in A(k,:).
    [K, N] = size(A);
    Am = permute(A, [2 3 1]);                  % N x 1 x K
    keep = ~Am & ~permute(Am, [2 1 3]);
    H = H .* keep;
end

function D = pagediag(d, I)
%PAGEDIAG  D(:,:,k) = diag(d(k,:)).
    D = permute(d, [2 3 1]) .* I;
end

function q = quad3(H, d)
%QUAD3  q(k) = d(k,:) * H(:,:,k) * d(k,:)'.
    K = size(H,3);
    q = zeros(K,1);
    dp = permute(d, [2 3 1]);
    if exist('pagemtimes','builtin') == 5 || exist('pagemtimes','file') == 2
        q = reshape(pagemtimes(dp, 'transpose', pagemtimes(H, dp), 'none'), K, 1);
    else
        for k = 1:K, q(k) = d(k,:) * H(:,:,k) * d(k,:).'; end
    end
end

function X = solve2(H, b1v, b2v)
%SOLVE2  X(k,:,1) = H_k \ b1v(k,:)'  and  X(k,:,2) = H_k \ b2v(k,:)'.
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

% =========================================================================
% SELF TEST
% =========================================================================
function pass = selftest(verbose)
%SELFTEST  Sampling, projection, Jacobian and solver checks.
%   pass = F.selftest()   prints one line per check and returns true if all pass.
    if nargin < 1, verbose = true; end
    P = defaults(); N = P.N; pass = true;
    rs = rand('state'); rand('state', 7); %#ok<RAND>
    chk = @(nm, c, extra) report_line(nm, c, extra, verbose);

    Z  = rand(200, N);
    TH = z2theta(Z, P);
    pass = chk('z2theta feasible', all(all(diff(TH,1,2) >= P.minGap - 1e-12)) && ...
        all(TH(:,1) >= P.lb-1e-12) && all(TH(:,end) <= P.ub+1e-12), '') && pass;
    pass = chk('z round-trip', max(max(abs(theta2z(TH,P) - Z))) < 1e-9, ...
        sprintf('err %.2e', max(max(abs(theta2z(TH,P)-Z))))) && pass;

    X  = P.lb + rand(200,N)*(P.ub-P.lb);
    Pr = project(X, P);
    pass = chk('project feasible', all(all(diff(Pr,1,2) >= P.minGap-1e-12)), '') && pass;
    pass = chk('project idempotent', max(max(abs(project(Pr,P)-Pr))) < 1e-12, '') && pass;

    % analytic Jacobian vs central differences
    Xs = z2theta(rand(4,N), P);
    [B, JB, G1] = bnjac(Xs, P);
    h = 1e-6; eJ = 0; eG = 0;
    for i = 1:N
        E = zeros(1,N); E(i) = h;
        eJ = max(eJ, max(max(abs(JB(:,:,i) - (bn_of(Xs+E,P)-bn_of(Xs-E,P))/(2*h)))));
        eG = max(eG, max(abs(G1(:,i) - (b1_of(Xs+E,P)-b1_of(Xs-E,P))/(2*h))));
    end
    pass = chk('d b_n/d th vs finite diff', eJ < 1e-6, sprintf('err %.2e', eJ)) && pass;
    pass = chk('d b_1/d th vs finite diff', eG < 1e-7, sprintf('err %.2e', eG)) && pass;
    pass = chk('bn vs bnjac agree', max(max(abs(B - bn_of(Xs,P)))) < 1e-14, '') && pass;

    % solver: fundamental met, ordering met, distortion plausible
    nfe_reset();
    [TH2, I] = solve(z2theta(rand(400,N),P), 0.60, P);
    n1 = nfe_get();
    pass = chk('solve: >=88% converged', mean(I.ok) >= 0.88, ...
        sprintf('%.1f%%', 100*mean(I.ok))) && pass;
    pass = chk('solve: b1 error < 1e-8', max(I.b1err(I.ok)) < 1e-8, ...
        sprintf('max %.2e', max(I.b1err(I.ok)))) && pass;
    pass = chk('solve: cost ~132 NFE/start', n1/400 < 200, ...
        sprintf('%.0f NFE/start', n1/400)) && pass;
    pass = chk('solve: finds a good minimum at Ma=0.60', min(I.thd) < 0.02, ...
        sprintf('best %.4f%%', 100*min(I.thd))) && pass;

    % spectrum() must agree with report() on the band they share
    Xs2 = z2theta(rand(8,N), P);
    Rr  = report(Xs2, 0.60, P);  Sp = spectrum(Xs2, P, 199);
    eS  = max(abs(Rr.THD49 - Sp.THD49));
    pass = chk('spectrum THD49 == report THD49', eS < 1e-12, ...
        sprintf('err %.2e', eS)) && pass;
    pass = chk('spectrum: WTHD <= THD on same band', ...
        all(Sp.WTHD49 <= Sp.THD49 + 1e-15), '') && pass;

    % minimax mode reaches a compliant point at Ma = 0.40
    [~, I2] = solve(z2theta(rand(96,N),P), 0.40, P, struct('mode','minimax'));
    pass = chk('minimax finds a grid-compliant point at Ma=0.40', ...
        min(I2.maxq) <= 1, sprintf('min MAXQ %.3f', min(I2.maxq))) && pass;

    rand('state', rs); %#ok<RAND>
    if verbose
        fprintf('\n  SELFTEST %s\n\n', tern(pass, 'PASSED', 'FAILED'));
    end
end

function ok = report_line(nm, c, extra, verbose)
    ok = logical(c);
    if verbose
        fprintf('  [%s] %-38s %s\n', tern(ok,'PASS','FAIL'), nm, extra);
    end
end

function s = tern(c,a,b), if c, s = a; else, s = b; end, end
