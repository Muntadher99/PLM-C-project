function T = gauss_newton_rate(csvfile, arms, budget, outcsv)
%GAUSS_NEWTON_RATE  Local convergence factor rho of undamped Gauss-Newton at
%   the reported solutions (paper, Section 2.2).
%
%   T = gauss_newton_rate('bench_main.csv')                 % PLM-C at 20000 NFE
%   T = gauss_newton_rate(csv, {'gn_hybrid'}, 20000, 'gauss_newton_rate.csv')
%
%   At a local minimiser theta* of  min 1/2||b||^2  s.t.  b_1 = Ma  and (4),
%   with active linear constraints A and multiplier mu of the fundamental,
%   let Z span null([grad b_1; A]). Undamped Gauss-Newton, linearised on the
%   active manifold, contracts with factor
%
%       rho = spectral radius of (Z'J'JZ) \ (Z'(S + mu*H_b1)Z),
%       S = sum_n b_n * Hess(b_n)                     [21, Sec. 10.3]
%
%   rho = 0: quadratic; 0 < rho < 1: linear at rate rho; rho >= 1: repelled.
%   All Hessians are diagonal in theta, so the computation is exact.
%   Also reported: the relative KKT residual, the smallest eigenvalue of the
%   reduced Hessian and the number of active constraints.
%
%   Output CSV columns: arm, seed, Ma, THD49_pct, rho, eigG_max, eigG_min,
%   lambda_min_Hred, kkt_rel, n_active, mu.

if nargin < 2 || isempty(arms),   arms = {'gn_hybrid'}; end
if nargin < 3 || isempty(budget), budget = 20000; end
if nargin < 4 || isempty(outcsv), outcsv = 'gauss_newton_rate.csv'; end
if ischar(arms), arms = {arms}; end

F = plm_core();
D = read_benchmark_csv(csvfile);
P = F.defaults(D.N, 49);

rows = [];  armcol = {};
for a = 1:numel(arms)
    k = find(strcmp(D.arm, arms{a}) & D.nfe_target == budget & D.accepted == 1);
    for r = k(:).'
        th = deg2rad(D.theta(r,:));
        s  = local_rate(th, D.Ma(r), P, F);
        rows(end+1,:) = [D.seed(r), D.Ma(r), D.THD49(r), s.rho, s.eigmax, s.eigmin, ...
                         s.lmin_red, s.kkt, s.nact, s.mu]; %#ok<AGROW>
        armcol{end+1,1} = arms{a}; %#ok<AGROW>
    end
end

hdr = {'arm','seed','Ma','THD49_pct','rho','eigG_max','eigG_min','lambda_min_Hred', ...
       'kkt_rel','n_active','mu'};
fid = fopen(outcsv, 'w');
fprintf(fid, '%s', hdr{1}); fprintf(fid, ',%s', hdr{2:end}); fprintf(fid, '\n');
for i = 1:size(rows,1)
    fprintf(fid, '%s', armcol{i}); fprintf(fid, ',%.10g', rows(i,:)); fprintf(fid, '\n');
end
fclose(fid);

fprintf('\n== gauss_newton_rate: undamped Gauss-Newton contraction factor at reported solutions ==\n');
fprintf('%-16s %5s %8s %8s %8s %8s %10s %10s\n', 'arm','n','med rho','p90 rho','max rho', ...
        'rho>=1','max KKT','min l(Hr)');
for a = 1:numel(arms)
    m = strcmp(armcol, arms{a});
    if ~any(m), continue; end
    rr = rows(m,4);
    fprintf('%-16s %5d %8.3f %8.3f %8.3f %7.1f%% %10.2e %10.2e\n', arms{a}, sum(m), ...
        median(rr), prctile_(rr,90), max(rr), 100*mean(rr >= 1), max(rows(m,8)), min(rows(m,7)));
end
fprintf('\nrho = 0: quadratic;  0<rho<1: linear at rate rho;  rho>=1: undamped GN repelled.\n');
fprintf('wrote %s\n\n', outcsv);
T = struct('hdr', {hdr}, 'arm', {armcol}, 'data', rows);
end

% -------------------------------------------------------------------------
function s = local_rate(th, Ma, P, F)
N = P.N;  th = sort(th(:).');
[B, JB, G1] = F.bnjac(th, P);                 % B 1xM, JB 1xMxN, G1 1xN
J   = reshape(JB, P.M, N);                    % M x N
b   = B(:);
ORD = P.ORD(:);  sg = P.s(:).';
% diagonal Hessians
Hn  = -(4/pi) * (ORD .* cos(ORD .* th)) .* sg;          % M x N, row n = diag of Hess(b_n)
Sd  = (b.' * Hn).';                                       % N x 1  (diag of S)
Hb1 = -(4/pi) * (sg .* cos(th)).';                        % N x 1  (diag of Hess(b1))
g   = J.' * b;                                            % grad of 1/2||b||^2
% Active linear constraints (rows of A in A*theta >= c): slack below tolA,
% then constraints with wrongly signed multipliers are released one by one.
tolA = 1e-5;                    % rad; the CSV angles carry 10 digits
Aall = zeros(0, N); slk = zeros(0,1);
e = zeros(1,N); e(1) = 1;  Aall(end+1,:) = e; slk(end+1) = th(1) - P.lb;
for i = 1:N-1
    e = zeros(1,N); e(i) = -1; e(i+1) = 1; Aall(end+1,:) = e; %#ok<AGROW>
    slk(end+1) = th(i+1) - th(i) - P.minGap; %#ok<AGROW>
end
e = zeros(1,N); e(N) = -1;  Aall(end+1,:) = e; slk(end+1) = P.ub - th(N);
act = find(slk < tolA);
for pass = 1:N
    A  = Aall(act,:);
    C  = [G1(:).'; A];
    ml = -(C.') \ g;                     % g + C' * ml = 0 (least squares)
    lam = -ml(2:end);                     % multipliers of A*theta >= c, want >= 0
    if isempty(lam) || all(lam >= -1e-10 * max(1, norm(g))), break; end
    [~, j] = min(lam);  act(j) = [];
end
mu = ml(1);
kkt = norm(g + C.'*ml) / max(norm(g), 1e-300);
Z  = null(C);
if isempty(Z)
    s = struct('rho',0,'eigmax',0,'eigmin',0,'lmin_red',Inf,'kkt',kkt,'nact',size(A,1),'mu',mu);
    return;
end
Hr = Z.' * (J.'*J) * Z;
Cr = Z.' * diag(Sd + mu*Hb1) * Z;
Hr = (Hr + Hr.')/2;  Cr = (Cr + Cr.')/2;
ev = real(eig(Hr \ Cr));
s.rho    = max(abs(ev));
s.eigmax = max(ev);
s.eigmin = min(ev);
s.lmin_red = min(eig(Hr + Cr));
s.kkt  = kkt;
s.nact = size(A,1);
s.mu   = mu;
end

function q = prctile_(x, p)
x = sort(x(:)); n = numel(x);
if n == 0, q = NaN; return; end
r = 1 + (n-1)*p/100; lo = floor(r); hi = ceil(r);
q = x(lo) + (r-lo)*(x(hi)-x(lo));
end
