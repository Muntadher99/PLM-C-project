function OUT = summarize_benchmark(csvfile, OPT)
%SUMMARIZE_BENCHMARK  Aggregate a run_benchmark CSV into summary tables.
%
%   S = summarize_benchmark('bench_main.csv')
%   S = summarize_benchmark(csvfile, struct('ref','gn_hybrid','print',false))
%
%   Aggregates only; every THD49, NFE and acceptance flag is read from the CSV
%   as written by evaluate_candidates.
%
%   OPT fields (defaults in brackets)
%     ref     reference arm of the paired comparisons ('' = gn_hybrid) ['']
%     print   print the tables                                   [true]
%     out     output base name ('' = <csv base>_summary)          ['']
%
%   Output <out>.xlsx (and one <out>_<sheet>.csv per sheet, always written):
%     summary    one row per (arm, budget): acceptance, mean / median / p90 /
%                worst / sd THD49 over accepted rows, median WTHD999,
%                compliant share, mean achieved NFE, achieved/target ratio,
%                median wall time, share of points where all seeds agree,
%                and a status flag (see notes)
%     by_budget  arm x budget pivots of mean THD49, median THD49, acceptance
%                and achieved/target NFE
%     by_arm     one row per arm, pooled over budgets
%     by_Ma      per (arm, Ma) at the largest budget
%     head2head  Wilcoxon signed-rank on per-(Ma, seed) THD49 differences
%                against the reference arm, with both arms' achieved NFE
%     notes      definitions and flagged cells
%   The per-seed paired statistics of the paper (Tables 4 and 5) are computed
%   by paper_tables.m.

if nargin < 1 || isempty(csvfile)
    error('summarize_benchmark:arg', 'usage: summarize_benchmark(''file.csv'')');
end
D.ref     = '';
D.xlsx    = true;
D.print   = true;
D.out     = '';
if nargin < 2, OPT = struct(); end
fn = fieldnames(OPT); C = D;
for i = 1:numel(fn), C.(fn{i}) = OPT.(fn{i}); end
if isempty(C.out)
    [pp, bb] = fileparts(csvfile);
    C.out = fullfile(pp, [bb '_summary']);
end

T = read_long_csv(csvfile);
n = numel(T.arm);
if n == 0, error('summarize_benchmark:empty', '%s has no data rows', csvfile); end

has_spwm = isfield(T, 'use_spwm');
budgets  = unique(T.nfe_target(isfinite(T.nfe_target)));
arms     = unique(T.arm);

% Order arms by median THD49 at the largest budget.
Bmax = max(budgets);
key  = nan(numel(arms),1);
for a = 1:numel(arms)
    m = strcmp(T.arm, arms{a}) & T.nfe_target == Bmax & T.accepted == 1;
    key(a) = nmedian(T.THD49_pct(m));
end
key(~isfinite(key)) = inf;
[~, ord] = sort(key);  arms = arms(ord);

% =========================================================================
% SHEET summary: one row per (arm, budget)
% =========================================================================
SH = {'name','arm_key','family','engine','strategy','spwm','nfe_target','n_points', ...
      'accept_pct','mean_THD49_pct','median_THD49_pct','p90_THD49_pct', ...
      'worst_THD49_pct','sd_THD49_pct','median_WTHD999_pct','compliant_pct', ...
      'mean_NFE_achieved','nfe_ratio','median_time_s','seed_identical_pct','status'};
SR = cell(0, numel(SH));
flagged = {};
for a = 1:numel(arms)
    for b = 1:numel(budgets)
        m = strcmp(T.arm, arms{a}) & T.nfe_target == budgets(b);
        if ~any(m), continue; end
        s = cell_stats(T, m);
        st = status_of(s);
        % A cell that did not run keeps its counts and cost but not its THD.
        if ~isempty(strfind(st, 'DID-NOT-RUN')) %#ok<STREMP>
            s.mean=NaN; s.med=NaN; s.p90=NaN; s.worst=NaN; s.sd=NaN;
            s.wthd=NaN; s.comp=NaN;
        end
        sp = 0;
        if has_spwm, sp = nmean(T.use_spwm(m)); end
        SR(end+1,:) = {paper_name(arms{a}), arms{a}, first_of(T.family,m), ...
            engine_name(first_of(T.engine,m)), ...
            first_of(T.strategy,m), sp, budgets(b), s.n, s.acc, s.mean, s.med, ...
            s.p90, s.worst, s.sd, s.wthd, s.comp, s.nfe, s.ratio, s.t, s.sid, st}; %#ok<AGROW>
        if ~strcmp(st,'ok')
            flagged{end+1,1} = sprintf('%-12s B=%-7g  %-18s (accept %.1f%%, %.2fx target, achieved NFE %.0f)', ...
                paper_name(arms{a}), budgets(b), st, s.acc, s.ratio, s.nfe); %#ok<AGROW>
        end
    end
end

% Seed check at the smallest budget (see seed_identical).
ci0 = @(nm) find(strcmp(SH, nm), 1);
seedbug = {};
Bmin = min(budgets);
for a = 1:numel(arms)
    r = find(strcmp(SR(:,2), arms{a}) & cell2mat(SR(:,ci0('nfe_target'))) == Bmin, 1);
    if isempty(r), continue; end
    v = SR{r, ci0('seed_identical_pct')};
    if isfinite(v) && v >= 95
        seedbug{end+1,1} = paper_name(arms{a}); %#ok<AGROW>
    end
end

% =========================================================================
% SHEET by_budget
% =========================================================================
PV = cell(0, numel(budgets)+1);
PV(end+1,:) = pad_row({'MEAN THD49 (%)  -- accepted points only'}, numel(budgets)+1);
PV(end+1,:) = [{'method'}, arrayfun(@(x) sprintf('%g',x), budgets(:).', 'UniformOutput', false)];
PV = [PV; pivot_block(T, arms, budgets, 'mean')];
PV(end+1,:) = pad_row({''}, numel(budgets)+1);
PV(end+1,:) = pad_row({'MEDIAN THD49 (%)  -- accepted points only'}, numel(budgets)+1);
PV(end+1,:) = [{'method'}, arrayfun(@(x) sprintf('%g',x), budgets(:).', 'UniformOutput', false)];
PV = [PV; pivot_block(T, arms, budgets, 'med')];
PV(end+1,:) = pad_row({''}, numel(budgets)+1);
PV(end+1,:) = pad_row({'ACCEPTED (%)'}, numel(budgets)+1);
PV(end+1,:) = [{'method'}, arrayfun(@(x) sprintf('%g',x), budgets(:).', 'UniformOutput', false)];
PV = [PV; pivot_block(T, arms, budgets, 'acc')];
PV(end+1,:) = pad_row({''}, numel(budgets)+1);
PV(end+1,:) = pad_row({'ACHIEVED NFE / TARGET'}, numel(budgets)+1);
PV(end+1,:) = [{'method'}, arrayfun(@(x) sprintf('%g',x), budgets(:).', 'UniformOutput', false)];
PV = [PV; pivot_block(T, arms, budgets, 'ratio')];

% =========================================================================
% SHEET by_arm: pooled over budgets
% =========================================================================
AH = {'name','arm_key','family','engine','strategy','n_points','accept_pct', ...
      'mean_THD49_pct','median_THD49_pct','p90_THD49_pct','worst_THD49_pct', ...
      'median_WTHD999_pct','compliant_pct','mean_nfe_ratio','median_time_s','p90_time_s'};
AR = cell(numel(arms), numel(AH));
for a = 1:numel(arms)
    m = strcmp(T.arm, arms{a});
    s = cell_stats(T, m);
    AR(a,:) = {paper_name(arms{a}), arms{a}, first_of(T.family,m), ...
        engine_name(first_of(T.engine,m)), ...
        first_of(T.strategy,m), s.n, s.acc, s.mean, s.med, s.p90, s.worst, ...
        s.wthd, s.comp, s.ratio, s.t, nquantile(T.time_s(m), 0.90)};
end

% =========================================================================
% SHEET by_Ma: per operating point at the largest budget
% =========================================================================
MH = {'name','nfe_target','Ma','n_seeds','accept_pct','mean_THD49_pct', ...
      'median_THD49_pct','worst_THD49_pct','median_WTHD999_pct'};
MR = cell(0, numel(MH));
MAS = unique(T.Ma(isfinite(T.Ma)));
for a = 1:numel(arms)
    for k = 1:numel(MAS)
        m = strcmp(T.arm, arms{a}) & T.nfe_target == Bmax & abs(T.Ma - MAS(k)) < 1e-9;
        if ~any(m), continue; end
        s = cell_stats(T, m);
        MR(end+1,:) = {paper_name(arms{a}), Bmax, MAS(k), s.n, s.acc, s.mean, s.med, ...
                       s.worst, s.wthd}; %#ok<AGROW>
    end
end

% =========================================================================
% SHEET head2head: paired against the reference arm
% =========================================================================
ref = C.ref;
if isempty(ref)
    if any(strcmp(arms,'gn_hybrid')), ref = 'gn_hybrid'; else, ref = arms{1}; end
end
% A row whose two arms differ by more than 15% in achieved NFE compares
% budgets, not methods, and is labelled as such (Section 4).
HH = {'budget','reference','other','n_pairs','n_possible', ...
      'ref_advantage_pp','ref_wins','ref_losses','ties','p_signrank', ...
      'significant','ref_NFE_spent','other_NFE_spent','budget_mismatch_x','comparison'};
HR = cell(0, numel(HH));
for b = 1:numel(budgets)
    for a = 1:numel(arms)
        if strcmp(arms{a}, ref), continue; end
        r = paired_vs(T, ref, arms{a}, budgets(b));
        if isempty(r), continue; end
        mm = max(r.nfeA, r.nfeB) / max(min(r.nfeA, r.nfeB), eps);
        if mm > 1.15
            verdict = sprintf('UNFAIR -- %s had %.1fx the evaluations', ...
                              paper_name(tern(r.nfeA > r.nfeB, ref, arms{a})), mm);
        else
            verdict = 'matched';
        end
        HR(end+1,:) = {budgets(b), paper_name(ref), paper_name(arms{a}), r.n, r.npos, r.adv, r.win, ...
                       r.lose, r.tie, r.p, tern(isfinite(r.p) && r.p < 0.05,'yes','no'), ...
                       r.nfeA, r.nfeB, mm, verdict}; %#ok<AGROW>
    end
end

% =========================================================================
% SHEET notes
% =========================================================================
NH = {'item','text'};
NR = {
 'source file',              csvfile
 'rows aggregated',          sprintf('%d', n)
 'arms',                     strjoin(cellfun(@paper_name, arms(:).', 'UniformOutput', false), ', ')
 'budgets (NFE per Ma)',     strjoin(arrayfun(@(x) sprintf('%g',x), budgets(:).', 'UniformOutput', false), ', ')
 'THD49',                    'RMS of odd non-triplen orders 5..49 relative to b1, in percent'
 'WTHD999',                  'sqrt(sum((b_n/n)^2))/|b1| to order 999'
 'acceptance rule',          'One rule, applied by evaluate_candidates to every arm: |b1-Ma| < 1e-8 and constraint violation < 1e-9'
 'averaging',                'THD statistics use accepted points only; accept_pct is on the same row.'
 'nfe_ratio',                'achieved / target evaluations per operating point.'
 'head2head sign',           sprintf('Positive ref_advantage_pp means %s has the lower THD49.', paper_name(ref))
 'head2head fairness',       'budget_mismatch_x is max/min of the two arms'' achieved NFE; above 1.15 the row is labelled UNFAIR.'
 'paired test',              'Wilcoxon signed-rank on per-(Ma, seed) differences, normal approximation with tie and continuity correction; only pairs where both arms were accepted.'
 'seed_identical_pct',       'Share of operating points at which every seed returned the identical THD49.'
 'status DID-NOT-RUN',       'accept < 5% or achieved NFE below 10% of target; THD statistics blanked.'
 'status OVERSPENT',         'achieved NFE above 1.10x target.'
 'status UNDERSPENT',        'achieved NFE below 0.85x target.'
 'status PARTIAL',           'accept below 80%.'
 };
NR(end+1,:) = {'--- GLOSSARY ---', 'short name  ->  full name  (arm key in the CSV)'};
for a = 1:numel(arms)
    [sh, fl] = paper_name(arms{a});
    NR(end+1,:) = {sh, sprintf('%s   [CSV arm = %s]', fl, arms{a})}; %#ok<AGROW>
end
if isempty(seedbug)
    NR(end+1,:) = {'SEED CHECK', 'passed -- every arm produced seed-to-seed variation at the smallest budget'};
else
    NR(end+1,:) = {'SEED CHECK FAILED', sprintf('%s returned the identical result for every seed at %g NFE/Ma.', ...
        strjoin(seedbug(:).', ', '), Bmin)};
end
if isempty(flagged)
    NR(end+1,:) = {'flagged cells', 'none'};
else
    for i = 1:numel(flagged)
        NR(end+1,:) = {sprintf('flagged %d', i), flagged{i}}; %#ok<AGROW>
    end
end

% =========================================================================
% WRITE
% =========================================================================
sheets = { 'summary',   SH, SR
           'by_budget', {}, PV
           'by_arm',    AH, AR
           'by_Ma',     MH, MR
           'head2head', HH, HR
           'notes',     NH, NR };
wrote = write_book(C.out, sheets);

OUT = struct('summary',{[SH; SR]}, 'by_budget',{PV}, 'by_arm',{[AH; AR]}, ...
             'by_Ma',{[MH; MR]}, 'head2head',{[HH; HR]}, 'notes',{[NH; NR]}, ...
             'file', C.out, 'ref', ref, 'flagged', {flagged}, 'seedbug', {seedbug});

if C.print
    print_summary(SH, SR, PV, HR, arms, budgets, ref, flagged, seedbug, csvfile, n, wrote, C.out);
end
end

% =========================================================================
% CONSOLE REPORT
% =========================================================================
function print_summary(SH, SR, PV, HR, arms, budgets, ref, flagged, seedbug, csvfile, n, wrote, base)
W = 118;
fprintf('\n%s\n', repmat('=',1,W));
fprintf('BENCHMARK SUMMARY   %s   (%d rows, %d arms, %d budgets)\n', csvfile, n, numel(arms), numel(budgets));
fprintf('%s\n', repmat('=',1,W));

ci = @(nm) find(strcmp(SH, nm), 1);
fprintf('\n%s\n', repmat('-',1,W));
fprintf('PER ARM AND BUDGET   (arms ordered best-first at %g NFE/Ma)\n', max(budgets));
fprintf('%s\n', repmat('-',1,W));
fprintf('%-12s %8s %6s %9s %9s %9s %9s %9s %6s %9s %6s %7s %9s  %s\n', ...
        'method','budget','acc%','meanTHD','medTHD','p90','worst','WTHD999','comp%','NFE','xtgt','time','seed==','status');
last = '';
for r = 1:size(SR,1)
    if ~strcmp(SR{r,1}, last) && ~isempty(last), fprintf('\n'); end
    last = SR{r,1};
    fprintf('%-12s %8g %5.1f%% %s %s %s %s %s %s %9.0f %6.2f %6.1fs %s  %s\n', ...
        SR{r,ci('name')}, SR{r,ci('nfe_target')}, SR{r,ci('accept_pct')}, ...
        pct(SR{r,ci('mean_THD49_pct')}), pct(SR{r,ci('median_THD49_pct')}), ...
        pct(SR{r,ci('p90_THD49_pct')}), pct(SR{r,ci('worst_THD49_pct')}), ...
        pct(SR{r,ci('median_WTHD999_pct')}), pct(SR{r,ci('compliant_pct')}, 0), ...
        SR{r,ci('mean_NFE_achieved')}, SR{r,ci('nfe_ratio')}, SR{r,ci('median_time_s')}, ...
        pct(SR{r,ci('seed_identical_pct')}, 0), SR{r,ci('status')});
end

fprintf('\n%s\n', repmat('-',1,W));
fprintf('MEAN THD49 (%%) BY BUDGET\n');
fprintf('%s\n', repmat('-',1,W));
print_block(PV, 1, numel(budgets));
fprintf('\nACCEPTED (%%) BY BUDGET\n');
k = find(strncmp(PV(:,1), 'ACCEPTED', 8), 1);
print_block(PV, k, numel(budgets));

if ~isempty(HR)
    fprintf('\n%s\n', repmat('-',1,W));
    fprintf('PAIRED vs %s   (+ve pp = %s is better; paired on (Ma, seed); * = p<0.05)\n', ...
            paper_name(ref), paper_name(ref));
    fprintf('%s\n', repmat('-',1,W));
    fprintf('%8s  %-12s %8s %10s %6s %6s %6s %11s   %8s %8s  %s\n', ...
            'budget','vs','pairs','adv (pp)','win','lose','tie','p','refNFE','othNFE','comparison');
    for r = 1:size(HR,1)
        star = ' '; if isfinite(HR{r,10}) && HR{r,10} < 0.05, star = '*'; end
        fprintf('%8g  %-12s %4d/%-3d %+9.3f %6d %6d %6d %10.1e%s  %8.0f %8.0f  %s\n', ...
            HR{r,1}, HR{r,3}, HR{r,4}, HR{r,5}, HR{r,6}, HR{r,7}, HR{r,8}, HR{r,9}, ...
            HR{r,10}, star, HR{r,12}, HR{r,13}, HR{r,15});
    end
end

if ~isempty(seedbug)
    fprintf('\nSEED CHECK FAILED: %s returned the identical result for all seeds at %g NFE/Ma.\n', ...
            strjoin(seedbug(:).', ', '), min(budgets));
end

if ~isempty(flagged)
    fprintf('\n%s\n', repmat('-',1,W));
    fprintf('FLAGGED CELLS  (%d)\n', numel(flagged));
    fprintf('%s\n', repmat('-',1,W));
    for i = 1:numel(flagged), fprintf('  %s\n', flagged{i}); end
end

if wrote
    fprintf('\nwrote %s.xlsx  (sheets: summary | by_budget | by_arm | by_Ma | head2head | notes)\n', base);
else
    fprintf('\nwrote %s_*.csv\n', base);
end
fprintf('%s\n\n', repmat('=',1,W));
end

function print_block(PV, k, nb)
% k points at the block title row; k+1 is the header row.
fprintf('%-12s', 'method');
for c = 2:nb+1, fprintf(' %10s', PV{k+1,c}); end
fprintf('\n');
r = k + 2;
while r <= size(PV,1) && ~isempty(PV{r,1}) && ~any(strncmp(PV{r,1}, {'MEAN','MEDIAN','ACCEPTED','ACHIEVED'}, 4))
    fprintf('%-12s', PV{r,1});
    for c = 2:nb+1
        v = PV{r,c};
        if ischar(v), fprintf(' %10s', v);
        elseif ~isfinite(v), fprintf(' %10s', '-');
        else, fprintf(' %10.4f', v);
        end
    end
    fprintf('\n');
    r = r + 1;
end
end

% =========================================================================
% NAMES USED IN THE PAPER
% =========================================================================
function [short, full] = paper_name(key)
%PAPER_NAME  Arm key in the CSV -> short and full name used in the paper.
switch lower(key)
  case 'gn_hybrid',            short='PLM-C';     full='Projected Levenberg-Marquardt with bidirectional continuation (proposed)';
  case 'gn_hybrid_zeroth_fwd', short='PLM-C/Z1P'; full='PLM corrector with the forward-only continuation of [2] (Table 5, row 2)';
  case 'gn_multistart',        short='PLM-MS';    full='Projected Levenberg-Marquardt, independent multi-start';
  case 'sqp_hybrid',           short='SQP-C';     full='fmincon SQP with bidirectional continuation';
  case 'sqp_multistart',       short='SQP-MS';    full='fmincon SQP, independent multi-start (same formulation and starts)';
  case 'fssa',                 short='FSSA';      full='Simulated annealing with continuation, reimplementation of [2]';
  case 'fssa_grid',            short='FSSA-grid'; full='FSSA with the weighted grid-code objective of [2] (Section 5.6)';
  case 'gwo_sqp_box',          short='GWO-SQP';   full='Grey wolf optimiser with SQP refinement [7, 27]';
  case 'multistart_tb',        short='MS-TB';     full='MultiStart, MATLAB Global Optimization Toolbox [28]';
  otherwise,                   short=key;         full=key;
end
end

function e = engine_name(e)
%ENGINE_NAME  Display label of the local engine.
switch lower(strtrim(e))
  case {'gauss-newton','gn'},  e = 'projected-LM';
  case {'fmincon-sqp','sqp'},  e = 'fmincon-SQP';
  case 'simulated-annealing',  e = 'simulated-annealing';
end
end

% =========================================================================
% AGGREGATION
% =========================================================================
function s = cell_stats(T, m)
%CELL_STATS  Statistics of one selection of rows. THD statistics are over
%   accepted rows; counts, NFE and time are over all rows.
ok  = m & T.accepted == 1 & isfinite(T.THD49_pct);
thd = T.THD49_pct(ok);
s.n     = sum(m);
s.acc   = 100 * sum(ok) / max(sum(m),1);
s.mean  = nmean(thd);
s.med   = nmedian(thd);
s.p90   = nquantile(thd, 0.90);
s.worst = nmax(thd);
s.sd    = nstd(thd);
s.wthd  = nmedian(T.WTHD999_pct(ok));
if isfield(T,'MAXQ'), s.comp = 100*nmean(double(T.MAXQ(ok) <= 1.0001)); else, s.comp = NaN; end
s.nfe   = nmean(T.nfe_achieved(m));
s.ratio = nmean(T.nfe_ratio(m));
s.t     = nmedian(T.time_s(m));
s.sid   = seed_identical(T, m);
end

function pc = seed_identical(T, m)
%SEED_IDENTICAL  Percentage of operating points at which every seed returned
%   the identical THD49. High at every budget would indicate that the seed
%   does not reach the arm; high only at the largest budget is convergence.
mas = unique(T.Ma(m));
if numel(mas) == 0, pc = NaN; return; end
same = 0;  used = 0;
for k = 1:numel(mas)
    mk = m & T.Ma == mas(k);
    v  = T.THD49_pct(mk);
    v  = v(isfinite(v));
    if numel(v) < 2, continue; end
    used = used + 1;
    if max(v) - min(v) == 0, same = same + 1; end
end
if used == 0, pc = NaN; else, pc = 100 * same / used; end
end

function st = status_of(s)
%STATUS_OF  Flags for cells that are not a matched, fully accepted comparison.
f = {};
if s.acc < 5 || (isfinite(s.ratio) && s.ratio < 0.10)
    f{end+1} = 'DID-NOT-RUN';
elseif s.acc < 80
    f{end+1} = 'PARTIAL';
end
if isfinite(s.ratio) && s.ratio > 1.10, f{end+1} = 'OVERSPENT';
elseif isfinite(s.ratio) && s.ratio < 0.85 && s.ratio >= 0.10, f{end+1} = 'UNDERSPENT';
end
if isempty(f), st = 'ok'; else, st = strjoin(f, '+'); end
end

function B = pivot_block(T, arms, budgets, what)
B = cell(numel(arms), numel(budgets)+1);
for a = 1:numel(arms)
    B{a,1} = paper_name(arms{a});
    for b = 1:numel(budgets)
        m = strcmp(T.arm, arms{a}) & T.nfe_target == budgets(b);
        if ~any(m), B{a,b+1} = NaN; continue; end
        s = cell_stats(T, m);
        switch what
          case 'mean',  B{a,b+1} = s.mean;
          case 'med',   B{a,b+1} = s.med;
          case 'acc',   B{a,b+1} = s.acc;
          case 'ratio', B{a,b+1} = s.ratio;
        end
        if any(strcmp(what,{'mean','med'})) && strncmp(status_of(s),'DID-NOT-RUN',11)
            B{a,b+1} = NaN;
        end
    end
end
end

function r = paired_vs(T, refarm, other, B)
%PAIRED_VS  Wilcoxon signed-rank on the per-(Ma, seed) THD49 differences.
r = [];
mR = strcmp(T.arm, refarm) & T.nfe_target == B;
mO = strcmp(T.arm, other)  & T.nfe_target == B;
if ~any(mR) || ~any(mO), return; end
kR = T.Ma(mR)*1e6 + T.seed(mR);   aR = T.THD49_pct(mR);  okR = T.accepted(mR)==1;
kO = T.Ma(mO)*1e6 + T.seed(mO);   aO = T.THD49_pct(mO);  okO = T.accepted(mO)==1;
[cm, iR, iO] = intersect(round(kR), round(kO));
npos = numel(cm);
good = okR(iR) & okO(iO) & isfinite(aR(iR)) & isfinite(aO(iO));
x = aR(iR(good));  y = aO(iO(good));
if numel(x) < 6, return; end
d = x - y;                                   % negative -> reference is lower
r.n = numel(d);  r.npos = npos;
r.adv  = -nmean(d);                          % positive -> reference is better
r.win  = sum(d < -1e-9);
r.lose = sum(d >  1e-9);
r.tie  = sum(abs(d) <= 1e-9);
r.p    = signrank_p(d);
r.nfeA = nmean(T.nfe_achieved(mR));
r.nfeB = nmean(T.nfe_achieved(mO));
end

function p = signrank_p(d)
%SIGNRANK_P  Two-sided Wilcoxon signed-rank test, normal approximation with
%   tie and continuity correction (no toolbox required).
d = d(isfinite(d) & abs(d) > 1e-12);
n = numel(d);
if n < 6, p = NaN; return; end
a = abs(d);
r = tiedrank_local(a);
W = sum(r(d > 0));
mu = n*(n+1)/4;
[~, ~, g] = unique(a);
tc = accumarray(g(:), 1);
v  = n*(n+1)*(2*n+1)/24 - sum(tc.^3 - tc)/48;
if v <= 0, p = NaN; return; end
z = (W - mu - 0.5*sign(W - mu)) / sqrt(v);
p = erfc(abs(z)/sqrt(2));
end

function r = tiedrank_local(x)
[xs, ix] = sort(x(:));  m = numel(xs);  rr = zeros(m,1);  i = 1;
while i <= m
    j = i;
    while j < m && xs(j+1) == xs(i), j = j + 1; end
    rr(i:j) = (i+j)/2;
    i = j + 1;
end
r = zeros(m,1);  r(ix) = rr;
end

% =========================================================================
% INPUT / OUTPUT
% =========================================================================
function T = read_long_csv(fn)
%READ_LONG_CSV  Header-driven CSV reader; text columns are detected from the
%   first data row.
fid = fopen(fn, 'r');
if fid < 0, error('summarize_benchmark:open', 'cannot open %s', fn); end
hdr = strtrim(fgetl(fid));
row1 = fgetl(fid);
if ~ischar(row1), fclose(fid); error('summarize_benchmark:empty','%s has no data rows', fn); end
H = strsplit(hdr, ',');
V = strsplit(row1, ',');
fmt = '';  istxt = false(1, numel(H));
for i = 1:numel(H)
    v = '';  if i <= numel(V), v = strtrim(V{i}); end
    if ~isempty(v) && isnan(str2double(v)), istxt(i) = true; end
end
for i = 1:numel(H)
    if istxt(i), fmt = [fmt '%s']; else, fmt = [fmt '%f']; end %#ok<AGROW>
end
frewind(fid);  fgetl(fid);                                  % skip header
Craw = textscan(fid, fmt, 'Delimiter', ',', 'EmptyValue', NaN, ...
                'CollectOutput', false, 'ReturnOnError', false);
fclose(fid);
T = struct();
for i = 1:numel(H)
    nm = valid_name(H{i});
    T.(nm) = Craw{i};
end
end

function s = valid_name(s)
%VALID_NAME  Header text -> struct field name.
s = regexprep(strtrim(s), '\W', '_');
if isempty(s) || ~isletter(s(1)), s = ['v' s]; end
end

function ok = write_book(base, sheets)
%WRITE_BOOK  One CSV per sheet always; an .xlsx when the host can write one.
for i = 1:size(sheets,1)
    wcell([base '_' sheets{i,1} '.csv'], [reshape_hdr(sheets{i,2}, sheets{i,3}); sheets{i,3}]);
end
ok = false;
fx = [base '.xlsx'];
if exist('writecell','file') == 2 || exist('writetable','file') == 2
    try
        if exist(fx,'file') == 2, delete(fx); end
        for i = 1:size(sheets,1)
            M = blankify([reshape_hdr(sheets{i,2}, sheets{i,3}); sheets{i,3}]);
            if exist('writecell','file') == 2
                writecell(M, fx, 'Sheet', sheets{i,1});
            else
                vn = cellfun(@valid_name, M(1,:), 'UniformOutput', false);
                writetable(cell2table(M(2:end,:), 'VariableNames', vn), ...
                           fx, 'Sheet', sheets{i,1});
            end
        end
        ok = true;
    catch e
        fprintf('  (xlsx skipped: %s -- the CSVs were still written)\n', e.message);
    end
end
end

function M = blankify(M)
%BLANKIFY  NaN -> empty cell for the spreadsheet.
for i = 1:numel(M)
    v = M{i};
    if isnumeric(v) && (isempty(v) || ~all(isfinite(v(:)))), M{i} = ''; end
end
end

function h = reshape_hdr(hdr, body)
if isempty(hdr), h = cell(0, size(body,2)); else, h = hdr(:).'; end
end

function wcell(fn, M)
fid = fopen(fn, 'w');
for r = 1:size(M,1)
    for c = 1:size(M,2)
        if c > 1, fprintf(fid, ','); end
        v = M{r,c};
        if ischar(v)
            if any(v == ','), fprintf(fid, '"%s"', v); else, fprintf(fid, '%s', v); end
        elseif isempty(v) || (isnumeric(v) && ~isfinite(v))
            % blank
        else
            fprintf(fid, '%.6g', double(v));
        end
    end
    fprintf(fid, '\n');
end
fclose(fid);
end

% =========================================================================
% NaN-SAFE STATISTICS
% =========================================================================
function y = nmean(x)
x = x(isfinite(x));
if isempty(x), y = NaN; else, y = sum(x)/numel(x); end
end

function y = nmax(x)
x = x(isfinite(x));
if isempty(x), y = NaN; else, y = max(x); end
end

function y = nstd(x)
x = x(isfinite(x));
if numel(x) < 2, y = NaN; else, y = std(x); end
end

function y = nmedian(x)
x = sort(x(isfinite(x)));  m = numel(x);
if m == 0, y = NaN; elseif mod(m,2), y = x((m+1)/2); else, y = 0.5*(x(m/2)+x(m/2+1)); end
end

function y = nquantile(x, p)
x = sort(x(isfinite(x)));  m = numel(x);
if m == 0, y = NaN; return; end
y = x(max(1, min(m, ceil(p*m))));
end

function r = pad_row(c, w)
r = cell(1, w);  r(:) = {''};  r(1:numel(c)) = c;
end

function v = first_of(c, m)
i = find(m, 1);
if isempty(i), v = ''; else, v = c{i}; end
end

function y = tern(c, a, b)
if c, y = a; else, y = b; end
end

function s = pct(v, d)
%PCT  Fixed-width percentage, or '--' for a blank cell.
if nargin < 2, d = 4; end
if isempty(v) || ~isnumeric(v) || ~isfinite(v)
    s = sprintf('%9s', '--');
else
    s = sprintf('%8.*f%%', d, v);
end
end
