function F = baseline_common()
%BASELINE_COMMON  Objective, constraints and fmincon refinement shared by the
%   two toolbox-based baselines, GWO-SQP (gwo_sqp_sweep.m) and MS-TB
%   (mstb_sweep.m).
%
%   F = baseline_common() returns a struct of function handles. Angle vectors
%   are 1 x N rows in radians unless stated otherwise.
%
%   Objective of both baselines (Table 2):
%       J(theta) = THD49 + RHO * sum_n WN_n * max(0, q_n - 1)^2,
%       q_n = |b_n| / (L_n |b_1|),  RHO = 30,
%   minimised subject to b_1 = Ma (nonlinear equality) and the ordering and
%   minimum-gap conditions of Eq. (4) as LINEAR inequalities A*theta <= b.
%
%   Handles
%     P = F.defaults()                  problem constants (N = 15, Section 2.1)
%     F.z2theta(Z, P)                   uniform map (0,1)^N -> feasible set
%     F.project_feasible(x, P)          Euclidean projection onto (4) (PAVA)
%     [A, b] = F.gap_A(N, P)            linear gap constraints
%     [J, g] = F.obj(theta, Ma, P)      penalised THD49 and its gradient (1 NFE)
%     J = F.obj_batch(TH, Ma, P)        vectorised objective for a population
%                                       (K NFE), with fundamental penalty 5e3*(b1-Ma)^2
%     F.nonlcon(theta, Ma, P, hard)     b_1 = Ma; hard = true adds q_n <= 1
%     R = F.report(theta, Ma, P)        THD49, b1 error, violations, MAXQ, compliance
%     S = F.spectrum(theta, nmax)       extended-band diagnostics (no NFE)
%     [a, ok, R, cost] = F.refine(x0, Ma, P, opts)
%                                       fmincon SQP refinement + Newton restoration
%     [a, ok, R, cost] = F.refine_compliant(x0, Ma, P, opts)
%                                       min THD49 with q_n <= 1 as hard constraints
%     o = F.cap_options(o, ceiling)     bound MaxFunctionEvaluations by the
%                                       evaluations left below a ceiling
%
%   Evaluation counter: the global variable NFE_COUNT, one per evaluation of
%   the harmonic vector (Section 4). F.obj, F.obj_batch and the restoration
%   charge it; the constraint function does not.
%
%   Requires the Optimization Toolbox (fmincon).

    F.defaults          = @defaults;
    F.z2theta           = @z2theta;
    F.project_feasible  = @project_feasible;
    F.gap_A             = @gap_A;
    F.obj               = @obj;
    F.obj_batch         = @obj_batch;
    F.nonlcon           = @nonlcon;
    F.report            = @report;
    F.spectrum          = @spectrum;
    F.refine            = @refine;
    F.refine_compliant  = @refine_compliant;
    F.cap_options       = @cap_options;
end

% =========================================================================
function P = defaults()
%DEFAULTS  Problem constants (Section 2.1) and the baselines' acceptance rule.
    P.N        = 15;
    P.minGap   = deg2rad(0.25);
    P.lb       = deg2rad(0.5);
    P.ub       = deg2rad(89.0);
    P.ORD      = odd_nontriplen(49);
    P.LIM      = grid_limits(P.ORD);
    P.W        = harm_weights(P.ORD);
    P.WN       = P.W ./ max(P.W);      % relative order weights w_n of the penalty
    P.RHO      = 30;                   % penalty weight
    P.FEAS_VIOL= 1e-6;                 % internal acceptance of a refinement
    P.FEAS_B1  = 1e-4;                 %   (the reported acceptance is external,
                                       %    see evaluate_candidates.m)
    P.MAXIT    = 200;
    P.COMP_TOL = 1e-4;                 % compliance tolerance on MAXQ
end

% =========================================================================
% FEASIBLE SET
% =========================================================================
function TH = z2theta(Z, P)
%Z2THETA  (0,1)^N -> feasible ordered angles (order-statistic map), K x N.
    if nargin < 2, P = defaults(); end
    Z = min(max(Z, 1e-12), 1 - 1e-12);
    [K, N] = size(Z);
    Sspan  = P.ub - P.lb - (N-1)*P.minGap;
    U    = zeros(K, N);
    prev = zeros(K, 1);
    for i = 1:N
        p    = 1/(N - i + 1);
        om   = (1 - Z(:,i)).^p;
        t    = 1 - om;
        prev   = prev + (1 - prev).*t;
        U(:,i) = prev;
    end
    TH = P.lb + P.minGap*(0:N-1) + Sspan*U;
end

function a = project_feasible(x, P)
%PROJECT_FEASIBLE  Euclidean projection onto
%   {lb <= a_1, a_{i+1} - a_i >= minGap, a_N <= ub}  (PAVA on the shifted vector).
    if nargin < 2, P = defaults(); end
    x     = x(:).';
    N     = numel(x);
    shift = P.minGap * (0:N-1);
    hi    = P.ub - P.minGap*(N-1);
    y     = pava(x - shift);
    y     = min(max(y, P.lb), hi);   % clipping preserves monotonicity
    a     = y + shift;
end

function v = pava(y)
%PAVA  Pool-adjacent-violators isotonic regression (non-decreasing fit).
    n = numel(y);
    val = zeros(1,n); wt = zeros(1,n); k = 0;
    for i = 1:n
        k = k + 1; val(k) = y(i); wt(k) = 1;
        while k > 1 && val(k) < val(k-1)
            val(k-1) = (val(k)*wt(k) + val(k-1)*wt(k-1)) / (wt(k) + wt(k-1));
            wt(k-1)  = wt(k) + wt(k-1);
            k = k - 1;
        end
    end
    v = zeros(1,n); j = 1;
    for b = 1:k
        m = wt(b); v(j:j+m-1) = val(b); j = j + m;
    end
end

function [A, b] = gap_A(N, P)
%GAP_A  A*theta <= b  encoding  theta_{i+1} - theta_i >= minGap.
    if nargin < 2, P = defaults(); end
    A = zeros(N-1, N);
    for i = 1:N-1
        A(i,i) = 1; A(i,i+1) = -1;
    end
    b = -P.minGap * ones(N-1,1);
end

% =========================================================================
% OBJECTIVE AND CONSTRAINTS
% =========================================================================
function [J, g] = obj(a, Ma, P)
%OBJ  J = THD49 + RHO * sum_n WN_n * max(0, q_n - 1)^2 and its analytic
%   gradient. `a` must be ordered (guaranteed by the linear constraints).
%   Charges one NFE per call.
    global NFE_COUNT %#ok<GVMIS>
    if isempty(NFE_COUNT), NFE_COUNT = 0; end
    NFE_COUNT = NFE_COUNT + 1;

    a  = a(:).';
    Nv = numel(a);
    s  = (-1).^(0:Nv-1);
    ORD = P.ORD(:); LIM = P.LIM(:).'; WN = P.WN(:).';

    b1  = (4/pi) * sum(s .* cos(a));
    bn  = ((4./(ORD*pi)) .* (cos(ORD*a) * s.')).';       % 1 x M
    ab1 = abs(b1); sg1 = sign(b1);
    sq  = sqrt(sum(bn.^2));
    thd = sq / ab1;

    r = abs(bn) / ab1;
    q = r ./ LIM;
    e = max(0, q - 1);
    J = thd + P.RHO * sum(WN .* e.^2);

    if nargout > 1
        Gn = -(4/pi) * (sin(ORD*a) .* s);                % M x N : dbn/da
        g1 = -(4/pi) * (s .* sin(a));                    % 1 x N : db1/da
        gTHD = (bn * Gn) / (sq*ab1) - (sq * sg1 / b1^2) * g1;
        coef = 2 * P.RHO * (WN .* e) ./ LIM;
        gPEN = ((coef .* sign(bn)) * Gn) / ab1 - (sum(coef .* r) * sg1 / ab1) * g1;
        g = (gTHD + gPEN).';
    end
end

function J = obj_batch(A, Ma, P, fundW)
%OBJ_BATCH  Vectorised objective for K x N ordered angle sets, with the
%   fundamental as a penalty fundW*(b1 - Ma)^2 (default 5e3). Charges K NFE.
    global NFE_COUNT %#ok<GVMIS>
    if isempty(NFE_COUNT), NFE_COUNT = 0; end
    NFE_COUNT = NFE_COUNT + size(A,1);
    if nargin < 4 || isempty(fundW), fundW = 5e3; end

    [K, Nv] = size(A);
    s   = (-1).^(0:Nv-1);
    b1  = (4/pi) * (cos(A) * s.');
    ab1 = abs(b1);
    num = zeros(K,1); pen = zeros(K,1);
    for i = 1:numel(P.ORD)
        n  = P.ORD(i);
        bn = (4/(n*pi)) * (cos(n*A) * s.');
        num = num + bn.^2;
        e   = max(0, abs(bn)./ab1./P.LIM(i) - 1);
        pen = pen + P.WN(i) * e.^2;
    end
    J = sqrt(num)./ab1 + P.RHO*pen + fundW*(b1 - Ma).^2;
end

function [c, ceq, gc, gceq] = nonlcon(a, Ma, P, hardLimits)
%NONLCON  Fundamental equality b_1 = Ma (the gaps are linear, in A*theta <= b).
%   hardLimits = true adds -L_n <= b_n/|b_1| <= L_n as inequality constraints.
    if nargin < 4, hardLimits = false; end
    a  = a(:).';
    Nv = numel(a);
    s  = (-1).^(0:Nv-1);
    ORD = P.ORD(:);

    b1  = (4/pi) * sum(s .* cos(a));
    ceq = b1 - Ma;
    c   = [];
    if hardLimits
        bn  = ((4./(ORD*pi)) .* (cos(ORD*a) * s.')).';
        r   = bn / abs(b1);
        c   = [ r - P.LIM(:).' , -r - P.LIM(:).' ].';
    end
    if nargout > 2
        g1   = -(4/pi) * (s .* sin(a));
        gceq = g1(:);
        gc   = [];
        if hardLimits
            Gn  = -(4/pi) * (sin(ORD*a) .* s);
            bn  = ((4./(ORD*pi)) .* (cos(ORD*a) * s.')).';
            ab1 = abs(b1); sg1 = sign(b1);
            Jr  = Gn/ab1 - (bn.'*g1)*sg1/b1^2;           % M x N  d(r)/da
            gc  = [Jr; -Jr].';                            % N x 2M
        end
    end
end

% =========================================================================
% REPORTING
% =========================================================================
function R = report(a, Ma, P)
%REPORT  THD49, b1, |b1 - Ma|, worst gap violation, ratios, MAXQ, compliance.
    a  = sort(a(:).');
    Nv = numel(a);
    s  = (-1).^(0:Nv-1);
    ORD = P.ORD(:);
    b1 = (4/pi) * sum(s .* cos(a));
    bn = ((4./(ORD*pi)) .* (cos(ORD*a) * s.')).';
    R.THD49   = sqrt(sum(bn.^2)) / abs(b1);
    R.b1      = b1;
    R.b1err   = abs(b1 - Ma);
    R.gapviol = max([0, P.minGap - diff(a)]);
    R.ratio   = abs(bn) / abs(b1);
    R.MAXQ    = max(R.ratio ./ P.LIM(:).');
    R.nviol   = sum(R.ratio > P.LIM(:).');
    R.viol    = max([R.gapviol, R.b1err]);
    if ~isfield(P,'COMP_TOL'), P.COMP_TOL = 1e-4; end
    R.compliant = R.MAXQ <= 1 + P.COMP_TOL;
end

function S = spectrum(a, nmax)
%SPECTRUM  Extended-band diagnostics of one angle set (odd non-triplen
%   orders): THD49, THD99, WTHD49, WTHD999, largest order above 49, b3 ratio.
%   Charges no NFE.
    if nargin < 2 || isempty(nmax), nmax = 999; end
    a  = sort(a(:)).';
    Nv = numel(a);
    sg = (-1).^(0:Nv-1);
    bnf = @(n) (4/(n*pi)) * sum(sg .* cos(n*a));
    a1  = abs(bnf(1));
    o   = 5:2:nmax;  o = o(mod(o,3) ~= 0);
    B   = zeros(1, numel(o));
    for i = 1:numel(o), B(i) = bnf(o(i)); end
    i49 = o <= 49;  i99 = o <= 99;
    S.THD49   = sqrt(sum(B(i49).^2)) / a1;
    S.THD99   = sqrt(sum(B(i99).^2)) / a1;
    S.WTHD49  = sqrt(sum((B(i49)./o(i49)).^2)) / a1;
    S.WTHD999 = sqrt(sum((B./o).^2)) / a1;
    oob = ~i49;
    if any(oob)
        [S.OOB_max, j] = max(abs(B(oob)) / a1);
        oo = o(oob);  S.OOB_order = oo(j);
    else
        S.OOB_max = 0;  S.OOB_order = NaN;
    end
    S.b3_ratio = abs(bnf(3)) / a1;
end

% =========================================================================
% LOCAL REFINEMENT
% =========================================================================
function [a, ok, R, cost] = refine(x0, Ma, P, opts)
%REFINE  fmincon SQP on F.obj with the linear gap constraints and b_1 = Ma,
%   followed by the Newton restoration onto {b_1 = Ma} (charged).
%   ok uses the internal rule viol < FEAS_VIOL and |b1 - Ma| < FEAS_B1;
%   cost is the NFE charged by this call.
    global NFE_COUNT %#ok<GVMIS>
    n0 = NFE_COUNT;

    N = P.N;
    x0 = project_feasible(x0, P);
    [A, b] = gap_A(N, P);
    if nargin < 4 || isempty(opts)
        opts = default_options(P);
    end
    fun = @(z) obj(z, Ma, P);
    con = @(z) nonlcon(z, Ma, P, false);
    a = x0;
    try
        a = fmincon(fun, x0, A, b, [], [], ...
                    P.lb*ones(1,N), P.ub*ones(1,N), con, opts);
    catch
    end
    a  = restore_b1(a, Ma, P);
    R  = report(a, Ma, P);
    ok = (R.viol < P.FEAS_VIOL) && (R.b1err < P.FEAS_B1);
    cost = NFE_COUNT - n0;
end

function [a, ok, R, cost] = refine_compliant(x0, Ma, P, opts)
%REFINE_COMPLIANT  min THD49 s.t. b_1 = Ma, linear gaps and every q_n <= 1
%   as hard constraints, followed by the Newton restoration (charged).
    global NFE_COUNT %#ok<GVMIS>
    n0 = NFE_COUNT;
    N  = P.N;
    x0 = project_feasible(x0, P);
    [A, b] = gap_A(N, P);
    if nargin < 4 || isempty(opts)
        opts = default_options(P, 'MaxIterations', 400, 'MaxFunctionEvaluations', 6e4);
    end
    P0 = P; P0.RHO = 0;                      % pure THD49; limits are constraints
    fun = @(z) obj(z, Ma, P0);
    con = @(z) nonlcon(z, Ma, P, true);
    a = x0;
    try
        a = fmincon(fun, x0, A, b, [], [], ...
                    P.lb*ones(1,N), P.ub*ones(1,N), con, opts);
    catch
    end
    a  = restore_b1(a, Ma, P);
    R  = report(a, Ma, P);
    ok = (R.viol < P.FEAS_VIOL) && (R.b1err < P.FEAS_B1) && (R.MAXQ <= 1 + P.COMP_TOL);
    cost = NFE_COUNT - n0;
end

function a = restore_b1(a, Ma, P)
%RESTORE_B1  Minimum-norm Newton correction onto {b_1 = Ma}, five steps, each
%   followed by the projection onto (4). Charges one NFE per step.
    global NFE_COUNT %#ok<GVMIS>
    if isempty(NFE_COUNT), NFE_COUNT = 0; end
    a  = project_feasible(a, P);
    Nv = numel(a);
    s  = (-1).^(0:Nv-1);
    for it = 1:5
        g = (4/pi)*sum(s .* cos(a)) - Ma;
        NFE_COUNT = NFE_COUNT + 1;
        if abs(g) < 1e-15, break; end
        G1 = -(4/pi) * (s .* sin(a));
        a  = project_feasible(a - (g/max(sum(G1.^2),1e-300))*G1, P);
    end
end

function o = default_options(P, varargin)
%DEFAULT_OPTIONS  fmincon SQP options with analytic gradients.
    o = optimoptions('fmincon','Algorithm','sqp','Display','none', ...
        'MaxIterations', P.MAXIT, 'MaxFunctionEvaluations', 3e4, ...
        'OptimalityTolerance',1e-12,'StepTolerance',1e-12, ...
        'ConstraintTolerance',1e-10, ...
        'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);
    if ~isempty(varargin), o = optimoptions(o, varargin{:}); end
end

function o = cap_options(o, ceilNow)
%CAP_OPTIONS  Bound one fmincon call by the evaluations left below ceilNow
%   (MaxFunctionEvaluations = max(20, ceilNow - NFE_COUNT) when that is lower).
    global NFE_COUNT %#ok<GVMIS>
    if ~isfinite(ceilNow), return; end
    remain = ceilNow - NFE_COUNT;
    cur = o.MaxFunctionEvaluations;
    if ~isnumeric(cur) || ~isscalar(cur), cur = Inf; end
    if remain < cur
        o = optimoptions(o, 'MaxFunctionEvaluations', max(20, floor(remain)));
    end
end

% =========================================================================
% PROBLEM TABLES
% =========================================================================
function o = odd_nontriplen(nMax)
    a = 5:2:nMax;  o = a(mod(a,3) ~= 0);
end

function L = grid_limits(orders)
%GRID_LIMITS  Table 1 limits as fractions of the fundamental.
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
%HARM_WEIGHTS  Order weights of the penalty (normalised to WN = W/max(W)).
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
