function T = minimax_witness_search(CFG, varargin)
%MINIMAX_WITNESS_SEARCH  Minimax compliance stage of Section 3.4 at every
%   operating point (paper, Section 5.7 and Table 7).
%
%   csvfile = minimax_witness_search()
%   csvfile = minimax_witness_search(struct('seeds',1:5,'NFE_TARGET',20000))
%   minimax_witness_search('-summary', 'minimax_search.csv')  % summary only
%
%   For every (Ma, seed) three objectives are run on the same engine
%   (plm_core solve), with the same allowance and the same random starts:
%     'thd'      min sum b_n^2                           s.t. b_1 = Ma
%     'comply'   min sum b_n^2 + rho*sum max(0, q_n - 1)^2   (weighted penalty)
%     'minimax'  min max_n q_n by four Lawson reweightings   s.t. b_1 = Ma
%   q_n = |b_n| / (L_n |b_1|). Each mode keeps the angle set that is best by
%   its own criterion (lowest THD49, lowest compliant THD49, lowest max q_n).
%   The minimax angle set is a compliance witness when max_n q_n <= 1 + 1e-4;
%   its loading Q^ = max_n q_n is an upper bound on the global minimax value
%   Q*, not Q* itself (Section 3.4).
%
%   CFG fields (defaults in brackets)
%     Ma_list      operating points                   [1.20:-0.05:0.20]
%     NFE_TARGET   allowance per Ma, per mode           [20000]
%     seeds        seed list; point i of seed s uses rng(1000*s + i) [1:5]
%     N, nmax      angles and highest counted order     [15, 49]
%     CHUNK        starts advanced per batched solve    [24]
%     out          output CSV ('' = minimax_search_<date>.csv) ['']
%     verbose      progress lines                       [true]
%
%   Output CSV, one row per (Ma, seed, mode): Ma, seed, mode, nfe_target,
%   nfe_achieved, starts, accepted, best_THD49_pct, best_MAXQ,
%   compliant_found, THD49_of_compliant_pct, MAXQ_of_best_THD,
%   THD49_at_bestMAXQ_pct, time_s, theta01_deg ... thetaNN_deg (the angle set
%   kept by that mode). A per-Ma summary (medians over the seeds) is written
%   to <out>_summary.xlsx (sheets by_Ma | headline | notes, CSV copies
%   always). build_witness_tables audits the angle sets and builds Table 7
%   from the best seed at each point.

% Summary of an existing run.
if nargin == 2 && ischar(CFG) && strcmpi(CFG, '-summary')
    F = plm_core();  summarise(varargin{1}, F.defaults());
    T = varargin{1};
    return;
end

% ------------------------------------------------------------------ config
D.Ma_list    = 1.20:-0.05:0.20;
D.NFE_TARGET = 20000;
D.seeds      = 1:5;
D.N          = 15;
D.nmax       = 49;
D.CHUNK      = 24;
D.out        = '';
D.verbose    = true;
if nargin < 1, CFG = struct(); end
fn = fieldnames(CFG); C = D;
for i = 1:numel(fn), C.(fn{i}) = CFG.(fn{i}); end
if isempty(C.out)
    c = clock(); %#ok<CLOCK>
    C.out = sprintf('minimax_search_%04d%02d%02d_%02d%02d.csv', c(1),c(2),c(3),c(4),c(5));
end

F   = plm_core();
P   = F.defaults(C.N, C.nmax);
MAS = C.Ma_list(:);
nM  = numel(MAS);
MODES = {'thd','comply','minimax'};

hdr = {'Ma','seed','mode','nfe_target','nfe_achieved','starts','accepted', ...
       'best_THD49_pct','best_MAXQ','compliant_found','THD49_of_compliant_pct', ...
       'MAXQ_of_best_THD','THD49_at_bestMAXQ_pct','time_s'};
hdr = [hdr, arrayfun(@(k) sprintf('theta%02d_deg',k), 1:C.N, 'UniformOutput', false)];
fid = fopen(C.out,'w');
fprintf(fid,'%s',hdr{1});
for i=2:numel(hdr), fprintf(fid,',%s',hdr{i}); end
fprintf(fid,'\n'); fclose(fid);

if C.verbose
    fprintf('\n=================== minimax_witness_search ===================\n');
    fprintf('Ma      : %g : %g : %g  (%d points)\n', MAS(1), MAS(2)-MAS(1), MAS(end), nM);
    fprintf('seeds   : %s\n', mat2str(C.seeds));
    fprintf('budget  : %g NFE per Ma per mode\n', C.NFE_TARGET);
    fprintf('modes   : %s\n', strjoin(MODES,', '));
    fprintf('compliance: max_n q_n <= 1 + %g\n', P.COMP_TOL);
    fprintf('out     : %s\n\n', C.out);
end

for mi = 1:nM
    Ma = MAS(mi);
    for sd = C.seeds(:).'
        for k = 1:numel(MODES)
            mode = MODES{k};
            rng_set(sd*1000 + mi);         % same random starts for every mode
            R = run_one(F, P, Ma, mode, C.NFE_TARGET, C.CHUNK);
            row = [Ma, sd, NaN, C.NFE_TARGET, R.nfe, R.starts, R.accepted, ...
                   100*R.bestTHD, R.bestMAXQ, double(R.comp_found), ...
                   100*R.thd_comp, R.maxq_at_thdmin, 100*R.thd_at_bestmaxq, ...
                   R.time_s, rad2deg(R.x)];
            append_row(C.out, mode, row);
        end
        if C.verbose
            fprintf('  Ma=%.2f seed=%d done\n', Ma, sd);
        end
    end
end

if C.verbose, fprintf('\nwrote %s\n', C.out); end
T = C.out;
summarise(C.out, P);
end

% =========================================================================
function R = run_one(F, P, Ma, mode, budget, chunk)
%RUN_ONE  Multi-start in one objective mode under a fixed NFE budget.
%   Batches are sized to the remaining budget from a per-mode cost estimate
%   (initial values 132 / 183 / 269 NFE per start, then measured within this
%   call only).
t0 = tic;
F.nfe_reset();
switch lower(mode)
  case 'thd',     cost = 132;
  case 'comply',  cost = 183;
  case 'minimax', cost = 269;
  otherwise,      cost = 200;
end
R = struct('bestTHD',Inf,'bestMAXQ',Inf,'comp_found',false,'thd_comp',NaN, ...
           'maxq_at_thdmin',NaN,'thd_at_bestmaxq',NaN,'x',nan(1,P.N), ...
           'starts',0,'accepted',0,'nfe',0);
xbest = nan(1,P.N);

while true
    remaining = budget - F.nfe();
    K   = min(chunk, floor(remaining / max(cost,1)));
    if K < 1, break; end                         % cannot fund another start
    n0  = F.nfe();
    TH0 = F.zsample(rand(K, P.N), P);            % uniform starts on the feasible set
    [TH, I] = F.solve(TH0, Ma, P, struct('mode', mode));
    cost = max(1, (F.nfe() - n0) / K);           % measured cost per start
    R.starts   = R.starts + K;
    R.accepted = R.accepted + sum(I.ok);

    ok = find(I.ok);
    for j = ok(:).'
        % lowest THD49 and its loading
        if I.thd(j) < R.bestTHD
            R.bestTHD        = I.thd(j);
            R.maxq_at_thdmin = I.maxq(j);
            if strcmpi(mode,'thd'), xbest = TH(j,:); end
        end
        % lowest loading max_n q_n and its THD49
        if I.maxq(j) < R.bestMAXQ
            R.bestMAXQ        = I.maxq(j);
            R.thd_at_bestmaxq = I.thd(j);
            if strcmpi(mode,'minimax'), xbest = TH(j,:); end
        end
        % lowest THD49 among compliant solutions
        if I.compliant(j) && (isnan(R.thd_comp) || I.thd(j) < R.thd_comp)
            R.comp_found = true;
            R.thd_comp   = I.thd(j);
            if strcmpi(mode,'comply'), xbest = TH(j,:); end
        end
    end
end

R.nfe    = F.nfe();
R.x      = xbest;
R.time_s = toc(t0);
if isinf(R.bestTHD),  R.bestTHD  = NaN; end
if isinf(R.bestMAXQ), R.bestMAXQ = NaN; end
end

% =========================================================================
function summarise(fn, P)
%SUMMARISE  Per-Ma summary (medians over seeds), printed and written.
rows = read_csv(fn);
if isempty(rows), fprintf('no rows to summarise\n'); return; end
MAS = unique_sorted([rows.Ma]);

fprintf('\n%s\n', repmat('=',1,104));
fprintf('GRID-CODE COMPLIANCE BY OPERATING POINT  (medians over seeds)\n');
fprintf('%s\n', repmat('=',1,104));
fprintf(['  Q^ any  : lowest max_n q_n found by any mode\n' ...
         '  Q^ mmax : lowest max_n q_n found by the minimax mode\n' ...
         '  Q^ <= 1 exhibits a compliant angle set; Q^ is an upper bound on Q*.\n\n']);
fprintf('%6s %9s %9s %9s %9s %11s %11s %10s %10s\n', ...
        'Ma','Q^ any','Q^ mmax','witness','THD@minTHD','MAXQ@minTHD','THD@Q^','price pp','penalty ok');
nCert = 0; nThdViol = 0; nGap = 0; nPenFail = 0; nMmaxLoses = 0;
outR = {};
for i = 1:numel(MAS)
    ma = MAS(i);
    qm = med_of(rows, ma, 'minimax', 'best_MAXQ');
    qd = med_of(rows, ma, 'thd',     'best_MAXQ');
    qc = med_of(rows, ma, 'comply',  'best_MAXQ');
    q  = min([qm qd qc]);
    td = med_of(rows, ma, 'thd',     'best_THD49_pct');
    qt = med_of(rows, ma, 'thd',     'MAXQ_of_best_THD');
    tc = med_of(rows, ma, 'comply',  'THD49_of_compliant_pct');
    pf = med_of(rows, ma, 'comply',  'compliant_found');
    cert = q <= 1 + P.COMP_TOL;
    if isfinite(qm) && isfinite(qd) && qm > qd + 1e-9, nMmaxLoses = nMmaxLoses + 1; end
    thdviol = qt > 1 + P.COMP_TOL;
    penok   = pf > 0.5;
    nCert = nCert + cert;  nThdViol = nThdViol + thdviol;
    nGap  = nGap  + (cert && thdviol);
    nPenFail = nPenFail + (cert && ~penok);
    tq = med_of(rows, ma, 'minimax', 'THD49_at_bestMAXQ_pct');
    fprintf('%6.2f %9.4f %9.4f %9s %10.4f%% %11.4f %10.4f%% %+9.3f %10s\n', ...
        ma, q, qm, tern(cert,'YES','no'), td, qt, tq, tq - td, tern(penok,'yes','NO'));
    outR(end+1,:) = {ma, q, qm, qd, double(cert), td, qt, tc, double(penok), tq, tq - td}; %#ok<AGROW>
end

fprintf('\n%s\n', repmat('-',1,104));
fprintf('  operating points                                              : %d\n', numel(MAS));
fprintf('  compliant witness exhibited (Q^ <= 1)                         : %d\n', nCert);
fprintf('  distortion-minimising solution violates a limit (MAXQ > 1)    : %d\n', nThdViol);
fprintf('  witness exhibited where the distortion minimum violates       : %d\n', nGap);
fprintf('  witness exhibited where the weighted penalty found none       : %d\n', nPenFail);
fprintf('  points where THD mode reached a lower MAXQ than minimax mode  : %d\n', nMmaxLoses);

write_outputs(fn, outR, struct('nMa',numel(MAS),'nCert',nCert,'nThdViol',nThdViol, ...
                               'nGap',nGap,'nPenFail',nPenFail,'nMmaxLoses',nMmaxLoses), P);
end

% =========================================================================
function write_outputs(fn, outR, K, P)
%WRITE_OUTPUTS  <fn>_summary.xlsx (by_Ma | headline | notes) and one
%   CSV per sheet (always written).
base = strrep(fn, '.csv', '_summary');

CH = {'Ma','Qhat_any_mode','Qhat_minimax_mode','Qhat_thd_mode','witness_found', ...
      'THD_at_minTHD_pct','MAXQ_at_minTHD','THD_at_comply_pct','penalty_found_compliant', ...
      'THD_at_Qhat_pct','price_of_compliance_pp', ...
      'minTHD_violates_limit','witness_but_minTHD_violates'};
CR = cell(size(outR,1), numel(CH));
for i = 1:size(outR,1)
    ma=outR{i,1}; q=outR{i,2}; qm=outR{i,3}; qd=outR{i,4}; ce=outR{i,5};
    td=outR{i,6}; qt=outR{i,7}; tc=outR{i,8}; pk=outR{i,9};
    tq=outR{i,10}; pr=outR{i,11};
    viol = double(qt > 1 + P.COMP_TOL);
    CR(i,:) = {ma, q, qm, qd, ce, td, qt, tc, pk, tq, pr, viol, double(ce>0.5 && viol>0.5)};
end

HH = {'quantity','value','meaning'};
HR = { 'operating points total',                K.nMa,        'points on the modulation-index grid'
       'compliance witness exhibited',          K.nCert,      'median over the seeds of Qhat <= 1 (Table 7 uses the best seed)'
       'min-THD optimum violates a limit',      K.nThdViol,   'MAXQ > 1 at the minimum-THD solution'
       'witness but min-THD violates',          K.nGap,       'compliance achieved by the minimax stage where the distortion minimum violates'
       'witness but penalty missed it',         K.nPenFail,   'the weighted-penalty formulation found no compliant point'
       'minimax beaten by THD mode on MAXQ',    K.nMmaxLoses, 'THD mode buys about twice the starts at matched NFE' };

NH = {'item','text'};
NR = {
 'source file',        fn
 'Q definition',       'Qhat = max_n q_n of the returned angle set, q_n = |b_n|/(L_n |b1|); the global minimax value Q* = min over feasible theta of max_n q_n s.t. b1 = Ma (Eq. (10)); Qhat is an upper bound on Q*'
 'Qhat <= 1',          'An explicit feasible angle set meets every limit: compliance is achievable at this operating point (constructive).'
 'Qhat > 1',           'No compliant set was found from the starts used; this is not a proof that none exists.'
 'Qhat any vs mmax',   'Qhat any is the lowest loading from any objective mode; Qhat mmax that of the minimax mode alone.'
 'matched budget',     'All three modes get the same NFE budget and the same random starts; per-start cost differs (thd ~132, comply ~183, minimax ~269 NFE).'
 'compliance test',    sprintf('MAXQ <= 1 + %g', P.COMP_TOL)
 'price of compliance','THD_at_Qhat_pct minus THD_at_minTHD_pct.'
 };

sheets = { 'by_Ma',       CH, CR
           'headline',    HH, HR
           'notes',       NH, NR };

for i = 1:size(sheets,1)
    wcell([base '_' sheets{i,1} '.csv'], [sheets{i,2}(:).'; sheets{i,3}]);
end

wrote = false;
fx = [base '.xlsx'];
if exist('writecell','file') == 2 || exist('writetable','file') == 2
    try
        if exist(fx,'file') == 2, delete(fx); end
        for i = 1:size(sheets,1)
            M = blankify([sheets{i,2}(:).'; sheets{i,3}]);
            if exist('writecell','file') == 2
                writecell(M, fx, 'Sheet', sheets{i,1});
            else
                vn = cellfun(@matlab_name, M(1,:), 'UniformOutput', false);
                writetable(cell2table(M(2:end,:), 'VariableNames', vn), fx, 'Sheet', sheets{i,1});
            end
        end
        wrote = true;
    catch e
        fprintf('  (xlsx skipped: %s -- the CSVs were still written)\n', e.message);
    end
end

if wrote
    fprintf('\nwrote %s   (sheets: by_Ma | headline | notes)\n', fx);
else
    fprintf('\nwrote %s_*.csv\n', base);
end
fprintf('%s\n\n', repmat('=',1,104));
end

function M = blankify(M)
%BLANKIFY  NaN -> empty cell for the spreadsheet.
for i = 1:numel(M)
    v = M{i};
    if isnumeric(v) && (isempty(v) || ~all(isfinite(v(:)))), M{i} = ''; end
end
end

function wcell(fn, M)
fid = fopen(fn, 'w');
for r = 1:size(M,1)
    for c = 1:size(M,2)
        if c > 1, fprintf(fid, ','); end
        v = M{r,c};
        if ischar(v)
            if any(v == ',') || any(v == '"')
                fprintf(fid, '"%s"', strrep(v, '"', '""'));
            else
                fprintf(fid, '%s', v);
            end
        elseif isempty(v) || (isnumeric(v) && ~isfinite(v))
            % blank
        else
            d = double(v);
            if d == fix(d) && abs(d) < 1e15
                fprintf(fid, '%d', d);          % integers as integers
            else
                fprintf(fid, '%.6f', d);
            end
        end
    end
    fprintf(fid, '\n');
end
fclose(fid);
end

% =========================================================================
% HELPERS
% =========================================================================
function append_row(fn, mode, row)
fid = fopen(fn,'a');
fprintf(fid, '%.6g,%.6g,%s', row(1), row(2), mode);
for c = 4:numel(row)
    if isnan(row(c)), fprintf(fid, ',');
    else,             fprintf(fid, ',%.10g', row(c));
    end
end
fprintf(fid, '\n');
fclose(fid);
end

function rows = read_csv(fn)
fid = fopen(fn,'r');  h = strsplit(strtrim(fgetl(fid)), ',');
rows = struct('Ma',{},'mode',{},'vals',{});
while true
    l = fgetl(fid);
    if ~ischar(l), break; end
    p = strsplit(strtrim(l), ',', 'CollapseDelimiters', false);
    r.Ma = str2double(p{1});  r.mode = p{3};
    r.vals = struct();
    for i = 1:min(numel(p), numel(h))
        nm = matlab_name(h{i});
        if i ~= 3, r.vals.(nm) = str2double(p{i}); end
    end
    rows(end+1) = r; %#ok<AGROW>
end
fclose(fid);
end

function v = med_of(rows, ma, mode, field)
x = [];
for i = 1:numel(rows)
    if abs(rows(i).Ma - ma) < 1e-9 && strcmpi(rows(i).mode, mode)
        f = matlab_name(field);
        if isfield(rows(i).vals, f), x(end+1) = rows(i).vals.(f); end %#ok<AGROW>
    end
end
x = sort(x(isfinite(x)));  n = numel(x);
if n == 0, v = NaN; elseif mod(n,2), v = x((n+1)/2); else, v = 0.5*(x(n/2)+x(n/2+1)); end
end

function s = matlab_name(s)
s = regexprep(strtrim(s), '\W', '_');
if isempty(s) || ~isletter(s(1)), s = ['v' s]; end
end

function u = unique_sorted(x)
u = sort(unique(x(isfinite(x))), 'descend');
end

function y = tern(c,a,b), if c, y = a; else, y = b; end, end

function rng_set(s)
try, rng(s,'twister'); catch, rand('state',s); randn('state',s); end %#ok<RAND>
end
