function S = paper_tables(source)
%PAPER_TABLES  Regenerate the tables and statistics of Sections 4 to 6 of the
%   paper from the result files of the package.
%
%   paper_tables()               % the result files shipped in */results/
%   paper_tables('reproduced')   % the files written by the run_* scripts in
%                                % */reproduced/ (results/ where none exists)
%   S = paper_tables(files)      % struct overriding any of the paths below
%
%   Printed, in order:
%     Table 3   median THD49 (%) by method and evaluation budget
%     Sec. 4    matched-budget check (ratio of median achieved NFE to PLM-C)
%     Table 4   Hodges-Lehmann shift of the per-seed differences, exact
%               Wilcoxon signed-rank p, Holm adjustment, W-L-T at 20000
%     Sec. 5.2  acceptance by budget, achieved NFE, wall time, per-point wins
%     Table 5   continuation ablation and corrector-by-strategy factorial,
%               with the paired contrasts quoted in Section 5.4
%     Table 6   corrector ablation at matched starts
%     Sec. 5.5  eleven and nineteen angles
%     Sec. 5.6  seed-to-seed spread and compliance at 20000 NFE
%     Table 8   device-level simulation against the analytic model
%     Table 9   distortion, compliance, switching and estimated loss
%   Table 7 is 6_compliance_witnesses/results/unified_table.csv
%   (build_witness_tables).
%
%   Paired statistics (Section 4): for every seed, the median THD49 over the
%   operating points accepted for BOTH methods; the ten per-seed differences
%   (other - PLM-C, positive favours PLM-C) are summarised by the
%   Hodges-Lehmann shift (median of the Walsh averages) and tested with the
%   exact two-sided Wilcoxon signed-rank test (zero differences dropped,
%   average ranks for ties).
%
%   Requires only base MATLAB (or GNU Octave); run setup_paths first.

root = fileparts(mfilename('fullpath'));
P.main     = {'3_benchmark_protocol', 'bench_main.csv'};
P.gwo      = {'3_benchmark_protocol', 'bench_gwo_sqp.csv'};
P.abl      = {'4_ablations', 'bench_forward_continuation.csv'};
P.sqpc     = {'4_ablations', 'bench_sqp_continuation.csv'};
P.abl_prod = {'4_ablations', 'corrector_ablation_production.csv'};
P.abl_conv = {'4_ablations', 'corrector_ablation_converged.csv'};
P.N11      = {'5_angle_count', 'bench_N11.csv'};
P.N19      = {'5_angle_count', 'bench_N19.csv'};
P.nonideal = {'7_device_level_validation', 'nonideal_results.csv'};
P.loss     = {'8_loss_analysis', 'loss_table.csv'};
use_rep = nargin >= 1 && ischar(source) && strcmpi(source, 'reproduced');
f = fieldnames(P);
for i = 1:numel(f)
    c = P.(f{i});
    p = fullfile(root, c{1}, 'results', c{2});
    if use_rep
        q = fullfile(root, c{1}, 'reproduced', c{2});
        if exist(q, 'file') == 2, p = q; else, fprintf('(no reproduced %s: stored results used)\n', c{2}); end
    end
    P.(f{i}) = p;
end
if nargin >= 1 && isstruct(source)
    f = fieldnames(source);
    for i = 1:numel(f), P.(f{i}) = source.(f{i}); end
end

M = cat_bench({read_benchmark_csv(P.main), read_benchmark_csv(P.gwo)});
budgets = [500 1000 2000 4000 8000 20000];
arms  = {'gn_hybrid','gn_multistart','sqp_multistart','gwo_sqp_box','multistart_tb','fssa'};
names = {'PLM-C','PLM-MS','SQP-MS','GWO-SQP','MS-TB','FSSA'};
S = struct();

% ======================= Table 3 ==========================================
banner('Table 3. Median THD49 (%) by evaluation budget (accepted rows)');
fprintf('%-9s', 'Method'); fprintf('%9g', budgets); fprintf('\n');
S.table3 = nan(numel(arms), numel(budgets));
for a = 1:numel(arms)
    fprintf('%-9s', names{a});
    for b = 1:numel(budgets)
        v = M.THD49(arm_rows(M, arms{a}, budgets(b)) & M.accepted == 1);
        S.table3(a,b) = med(v);
        fprintf('%9s', num_or_dash(S.table3(a,b), 4));
    end
    fprintf('\n');
end

% ======================= matched budgets ==================================
banner('Section 4. Median achieved NFE relative to PLM-C (max/min; > 1.15 is excluded)');
fprintf('%-9s', 'Method'); fprintf('%9g', budgets); fprintf('\n');
for a = 2:numel(arms)
    fprintf('%-9s', names{a});
    for b = 1:numel(budgets)
        r0 = med(M.nfe_achieved(arm_rows(M, arms{1}, budgets(b))));
        r1 = med(M.nfe_achieved(arm_rows(M, arms{a}, budgets(b))));
        if min(r0, r1) > 0, fprintf('%9.3f', max(r0, r1) / min(r0, r1));
        else,               fprintf('%9s', '-'); end           % no design returned
    end
    fprintf('\n');
end

% ======================= Table 4 ==========================================
banner('Table 4. Advantage of PLM-C (pp): Hodges-Lehmann shift of per-seed differences');
B4 = [1000 2000 4000 8000 20000];
fprintf('%-9s', 'Baseline'); fprintf('%16g', B4); fprintf('%10s\n', 'W-L-T');
pv = [];  lab = {};
S.table4 = nan(numel(arms)-1, numel(B4));
for a = 2:numel(arms)
    fprintf('%-9s', names{a});
    for b = 1:numel(B4)
        [hl, p] = contrast(M, 'gn_hybrid', arms{a}, B4(b));
        S.table4(a-1,b) = hl;
        fprintf('%9.3f (%.3f)', hl, p);
        pv(end+1) = p; lab{end+1} = sprintf('%s @ %g', names{a}, B4(b)); %#ok<AGROW>
    end
    [~, ~, w, l, t] = contrast(M, 'gn_hybrid', arms{a}, 20000);
    fprintf('%6d-%d-%d\n', w, l, t);
end
padj = holm(pv);
fprintf('\n(p) is the exact two-sided Wilcoxon p over ten seeds; %d of %d reach p <= 0.002;\n', ...
        sum(pv <= 0.002 + 1e-9), numel(pv));
fprintf('%d of %d are significant at 0.05 after Holm adjustment over the %d comparisons.\n', ...
        sum(padj < 0.05), numel(pv), numel(pv));

% ======================= Section 5.2 ======================================
banner('Section 5.2. Accepted cells (%), mean achieved NFE, median sweep time');
fprintf('%-9s', 'Method'); fprintf('%8g', budgets); fprintf('%12s %12s %10s\n', 'NFE@500', 'NFE@20000', 'time@20k');
for a = 1:numel(arms)
    fprintf('%-9s', names{a});
    for b = 1:numel(budgets)
        m = arm_rows(M, arms{a}, budgets(b));
        fprintf('%8.1f', 100*sum(M.accepted(m) == 1)/max(sum(m),1));
    end
    fprintf('%12.0f %12.0f %9.1fs\n', mean(M.nfe_achieved(arm_rows(M, arms{a}, 500))), ...
        mean(M.nfe_achieved(arm_rows(M, arms{a}, 20000))), med(M.time_s(arm_rows(M, arms{a}, 20000))));
end
mas = unique(round(M.Ma*1e4))/1e4;  nwin = 0;  notes = {};
for k = 1:numel(mas)
    best = Inf; who = '';
    for a = 1:numel(arms)
        m = arm_rows(M, arms{a}, 2000) & M.accepted == 1 & abs(M.Ma - mas(k)) < 1e-6;
        v = med(M.THD49(m));
        if v < best, best = v; who = names{a}; end
        if a == 1, v1 = v; end
    end
    if strcmp(who, 'PLM-C'), nwin = nwin + 1;
    else, notes{end+1} = sprintf('Ma = %.2f: %s lower by %.2f pp', mas(k), who, v1 - best); end %#ok<AGROW>
end
fprintf('\nPer operating point at 2000 NFE, PLM-C has the lowest median THD49 at %d of %d.\n', nwin, numel(mas));
for i = 1:numel(notes), fprintf('  %s\n', notes{i}); end

% ======================= Table 5 ==========================================
banner('Table 5. Continuation ablation and corrector-by-strategy factorial: median THD49 (%)');
A = cat_bench({M, read_benchmark_csv(P.abl), read_benchmark_csv(P.sqpc)});
B5 = [1000 2000 8000 20000];
cfg5 = {'gn_multistart','Multi-start, no continuation'; ...
        'gn_hybrid_zeroth_fwd','Zeroth predictor, forward'; ...
        'gn_hybrid','Zeroth predictor, bidirectional'; ...
        'sqp_multistart','SQP, multi-start (SQP-MS)'; ...
        'sqp_hybrid','SQP, bidirectional (SQP-C)'};
fprintf('%-34s', 'Configuration'); fprintf('%9g', B5); fprintf('\n');
for c = 1:size(cfg5,1)
    fprintf('%-34s', cfg5{c,2});
    for b = 1:numel(B5)
        fprintf('%9.3f', med(A.THD49(arm_rows(A, cfg5{c,1}, B5(b)) & A.accepted == 1)));
    end
    fprintf('\n');
end
fprintf('\nSection 5.4 contrasts (HL shift in pp, positive favours the first; exact p; W-L-T):\n');
pairs = {'gn_hybrid_zeroth_fwd','gn_multistart','forward continuation [2] vs multi-start'; ...
         'gn_hybrid','gn_hybrid_zeroth_fwd','reverse pooled pass (bidirectional vs forward)'; ...
         'sqp_hybrid','sqp_multistart','continuation with the SQP corrector'; ...
         'gn_hybrid','sqp_hybrid','PLM vs SQP corrector inside continuation'; ...
         'gn_multistart','sqp_multistart','PLM vs SQP corrector inside multi-start'; ...
         'gn_hybrid','gn_multistart','continuation with the PLM corrector'};
for q = 1:size(pairs,1)
    fprintf('  %s\n', pairs{q,3});
    for B = [1000 2000 4000 8000 20000]
        if ~any(arm_rows(A, pairs{q,1}, B)) || ~any(arm_rows(A, pairs{q,2}, B)), continue; end
        [hl, p, w, l, t] = contrast(A, pairs{q,1}, pairs{q,2}, B);
        fprintf('    %6g NFE: %+7.3f pp   p = %.3f   %d-%d-%d\n', B, hl, p, w, l, t);
    end
end

% ======================= Table 6 ==========================================
banner('Table 6. Corrector ablation at matched starts (means over five operating points)');
Tp = read_simple_csv(P.abl_prod);  Tc = read_simple_csv(P.abl_conv);
vars = {'proposed','gn','lm1','theta','batch1'};
vlab = {'Proposed','No damping','One metric','Angle space','Unbatched'};
fprintf('%-12s %8s %8s %12s %12s %8s %10s\n', 'Variant','Acc.(%)','KKT(%)','Median THD','Best THD','NFE','time (s)');
tprop = sum(col_of(Tp, 'time_s', strcmp(Tp.variant, 'proposed')));
for i = 1:numel(vars)
    mp = strcmp(Tp.variant, vars{i});  mc = strcmp(Tc.variant, vars{i});
    kkt = NaN;  if any(mc), kkt = mean(col_of(Tc, 'kkt_pct', mc)); end
    fprintf('%-12s %8.1f %8s %12.2f %12.3f %8.0f %10.2f\n', vlab{i}, ...
        mean(col_of(Tp,'accepted_pct',mp)), num_or_dash(kkt, 1), ...
        mean(col_of(Tp,'median_THD49_pct',mp)), mean(col_of(Tp,'best_THD49_pct',mp)), ...
        mean(col_of(Tp,'nfe_per_start',mp)), sum(col_of(Tp,'time_s',mp)));
end
fprintf('KKT share from the run at ten times the iteration caps; unbatched\n');
fprintf('wall time / batched wall time = %.1f.\n', ...
    sum(col_of(Tp,'time_s',strcmp(Tp.variant,'batch1'))) / tprop);

% ======================= Section 5.5 ======================================
banner('Section 5.5. Eleven and nineteen angles');
for NN = [11 19]
    if NN == 11, Dn = read_benchmark_csv(P.N11); else, Dn = read_benchmark_csv(P.N19); end
    fprintf('N = %d\n', NN);
    for a = {'gn_hybrid','gn_multistart','sqp_multistart'}
        fprintf('  %-15s', a{1});
        for B = [2000 20000]
            m = arm_rows(Dn, a{1}, B);
            fprintf('  median THD49 %8.4f%% (compliant %4.1f%%) @%g', med(Dn.THD49(m & Dn.accepted == 1)), ...
                100*sum(Dn.compliant(m) == 1 & Dn.accepted(m) == 1)/max(sum(m),1), B);
        end
        fprintf('\n');
    end
    for b = {'gn_multistart','sqp_multistart'}
        for B = [2000 20000]
            [hl, p, w, l, t] = contrast(Dn, 'gn_hybrid', b{1}, B);
            fprintf('  PLM-C vs %-15s @%5g: %+8.4f pp  p = %.3f  %d-%d-%d\n', b{1}, B, hl, p, w, l, t);
        end
    end
end

% ======================= Section 5.6 ======================================
banner('Section 5.6. Seed-to-seed spread (pp) and compliance at 20000 NFE');
fprintf('%-9s %14s %12s %14s\n', 'Method', 'median spread', 'max spread', 'compliant (%)');
for a = 1:numel(arms)
    m = arm_rows(M, arms{a}, 20000);
    sp = [];
    for k = 1:numel(mas)
        v = M.THD49(m & M.accepted == 1 & abs(M.Ma - mas(k)) < 1e-6);
        if ~isempty(v), sp(end+1) = max(v) - min(v); end %#ok<AGROW>
    end
    fprintf('%-9s %14.2f %12.2f %14.1f\n', names{a}, med(sp), max(sp), ...
        100*sum(M.accepted(m) == 1 & M.compliant(m) == 1)/max(sum(m),1));
end

% ======================= Table 8 ==========================================
banner('Table 8. Device-level simulation against the analytic model, ten angle sets');
Tn = read_simple_csv(P.nonideal);
cases = unique(Tn.case, 'stable');
fprintf('%-26s %8s %8s %8s %9s %9s %8s\n', 'Case', '|dTHD| med', 'max', 'q_max', 'b1 med', 'b1 worst', 'even');
for i = 1:numel(cases)
    m  = strcmp(Tn.case, cases{i});
    mc = m & strcmp(Tn.design, 'compliant');
    d  = abs(col_of(Tn, 'dTHD49_pp', m));  b = col_of(Tn, 'b1_err_pct_sim', m);
    fprintf('%-26s %8.3f %8.3f %8.3f %9.2f %9.2f %8.2f\n', cases{i}, med(d), max(d), ...
        max(col_of(Tn, 'maxq_sim', mc)), med(b), min(b), max(col_of(Tn, 'even_max_pct', m)));
end
fprintf('switching-function model vs simulation: median |dTHD49| = %.3f pp\n', ...
    med(abs(col_of(Tn,'THD49_swfn',true(size(Tn.case))) - col_of(Tn,'THD49_sim',true(size(Tn.case))))));

% ======================= Table 9 ==========================================
banner('Table 9. Distortion, compliance, switching and estimated loss on the 21-point grid');
Tl = read_simple_csv(P.loss);
des = unique(Tl.design, 'stable');
fprintf('%-18s %8s %9s %7s %8s %9s %8s\n', 'Design', 'fsw Hz', 'THD49 %', 'Compl.', 'q_max', 'Loss W', 'eta %');
for i = 1:numel(des)
    print_loss_row(Tl, strcmp(Tl.design, des{i}), des{i});
end
% deployed design: THD-min where it complies, otherwise the compliant witness
mt = strcmp(Tl.design, 'PLM-C THD-min');  mcmp = strcmp(Tl.design, 'PLM-C compliant');
ma_t = Tl.Ma(mt);  ok_t = col_of(Tl, 'compliant', mt) == 1;
dep = false(size(Tl.Ma));
it = find(mt);  ic = find(mcmp);
for k = 1:numel(it)
    if ok_t(k), dep(it(k)) = true;
    else, dep(ic(abs(Tl.Ma(ic) - ma_t(k)) < 1e-9)) = true; end
end
print_loss_row(Tl, dep, 'Deployed');
fprintf('Medians of THD49 and q_max; loss and efficiency averaged over the grid.\n\n');
end

% =========================================================================
% PAIRED STATISTICS
% =========================================================================
function [hl, p, w, l, t] = contrast(D, a, b, B)
%CONTRAST  Per-seed medians over the points both arms solve; d = b - a.
seeds = unique(D.seed(arm_rows(D, a, B)));
d = [];
for s = seeds(:).'
    ma = arm_rows(D, a, B) & D.seed == s & D.accepted == 1;
    mb = arm_rows(D, b, B) & D.seed == s & D.accepted == 1;
    [common, ia, ib] = intersect(round(D.Ma(ma)*1e4), round(D.Ma(mb)*1e4));
    if isempty(common), continue; end
    va = D.THD49(ma);  vb = D.THD49(mb);
    d(end+1) = med(vb(ib)) - med(va(ia)); %#ok<AGROW>
end
hl = hodges_lehmann(d);
p  = signrank_exact(d);
w  = sum(d > 1e-9);  l = sum(d < -1e-9);  t = numel(d) - w - l;
end

function h = hodges_lehmann(d)
%HODGES_LEHMANN  Median of the Walsh averages (d_i + d_j)/2, i <= j.
n = numel(d);  W = [];
for i = 1:n
    for j = i:n
        W(end+1) = (d(i) + d(j))/2; %#ok<AGROW>
    end
end
h = med(W);
end

function p = signrank_exact(d)
%SIGNRANK_EXACT  Exact two-sided Wilcoxon signed-rank p-value by enumerating
%   all 2^n sign assignments (zero differences dropped, average ranks).
d = d(abs(d) > 1e-12);  n = numel(d);
if n == 0, p = NaN; return; end
r  = tied_rank(abs(d));
Wp = sum(r(d > 0));  mu = sum(r)/2;
S  = dec2bin(0:2^n-1, n) - '0';            % all sign patterns
W  = S * r(:);
p  = mean(abs(W - mu) >= abs(Wp - mu) - 1e-9);
end

function r = tied_rank(x)
[xs, ix] = sort(x(:));  m = numel(xs);  rr = zeros(m,1);  i = 1;
while i <= m
    j = i;
    while j < m && xs(j+1) == xs(i), j = j + 1; end
    rr(i:j) = (i+j)/2;
    i = j + 1;
end
r = zeros(m,1);  r(ix) = rr;
end

function q = holm(p)
%HOLM  Holm step-down adjusted p-values.
[ps, o] = sort(p(:));  m = numel(ps);
adj = min(1, cummax((m - (1:m)' + 1) .* ps));
q = zeros(size(p));  q(o) = adj;
end

% =========================================================================
% DATA HELPERS
% =========================================================================
function m = arm_rows(D, arm, B)
m = strcmp(D.arm, arm) & D.nfe_target == B;
end

function D = cat_bench(list)
%CAT_BENCH  Concatenate read_benchmark_csv structs (same N).
D = list{1};
f = {'arm','nfe_target','nfe_achieved','seed','Ma','accepted','THD49','MAXQ', ...
     'compliant','time_s','theta'};
for k = 2:numel(list)
    for i = 1:numel(f), D.(f{i}) = [D.(f{i}); list{k}.(f{i})]; end
end
end

function T = read_simple_csv(fn)
%READ_SIMPLE_CSV  Header-driven reader. A column is numeric when every
%   non-empty cell parses as a number (or NaN); otherwise it is a cellstr.
txt = fileread(fn);  L = regexp(txt, '\r?\n', 'split');  L = L(~cellfun(@isempty, strtrim(L)));
H = regexp(strtrim(L{1}), ',', 'split');  n = numel(L) - 1;
C = repmat({''}, n, numel(H));
for r = 1:n
    v = regexp(L{r+1}, ',', 'split');  k = min(numel(v), numel(H));
    C(r, 1:k) = v(1:k);
end
T = struct();
for c = 1:numel(H)
    x = str2double(C(:,c));
    txtcell = isnan(x) & ~cellfun(@isempty, strtrim(C(:,c))) & ~strcmpi(strtrim(C(:,c)), 'nan');
    nm = regexprep(strtrim(H{c}), '\W', '_');
    if any(txtcell), T.(nm) = C(:,c); else, T.(nm) = x; end
end
end

function v = col_of(T, name, m)
x = T.(name);  v = x(m);
end

function print_loss_row(T, m, label)
fprintf('%-18s %8.0f %9.3f %4d/%-2d %8.3f %9.2f %8.2f\n', label, ...
    med(col_of(T,'fsw_dev_Hz',m)), med(col_of(T,'THD49_pct',m)), sum(col_of(T,'compliant',m)), ...
    sum(m), med(col_of(T,'maxq',m)), mean(col_of(T,'P_loss_W',m)), mean(col_of(T,'eta_pct',m)));
end

function y = med(x)
x = sort(x(isfinite(x)));  n = numel(x);
if n == 0, y = NaN; elseif mod(n,2), y = x((n+1)/2); else, y = 0.5*(x(n/2)+x(n/2+1)); end
end

function s = num_or_dash(v, d)
if ~isfinite(v), s = '-'; else, s = sprintf('%.*f', d, v); end
end

function banner(s)
fprintf('\n%s\n%s\n', s, repmat('-', 1, numel(s)));
end
