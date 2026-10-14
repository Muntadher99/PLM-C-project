function OUT = export_validation_cases(csvfile, CFG)
%EXPORT_VALIDATION_CASES  Extract the angle sets of the ideal-switch Simulink
%   check from a benchmark CSV (paper, Section 6.2: "on the angles of all six
%   methods at thirty cases").
%
%   export_validation_cases('bench_all.csv')
%   export_validation_cases(csvfile, struct('Ma_list',[1.20 0.80 0.40]))
%
%   The angles are taken from the benchmark rows themselves (no re-solving),
%   so the continuation arms keep the path of the full sweep. Build the input
%   with all six arms first:
%       merge_benchmark_csv({'bench_main.csv','bench_gwo_sqp.csv'}, 'bench_all.csv')
%
%   CFG fields (defaults in brackets)
%     budget         nfe_target to select                   [20000]
%     seed           seed to select                          [1]
%     Ma_list        operating points                        [1.20 1.00 0.80 0.60 0.40]
%     arms           arm keys            [gn_hybrid gn_multistart sqp_multistart
%                                         gwo_sqp_box multistart_tb fssa]
%     out            output workbook                         ['validation_cases.xlsx']
%     N              number of angles                        [15]
%     only_accepted  export accepted rows only               [true]
%     verbose        print the exported table                [true]
%
%   Output workbook (CSV files per sheet if no spreadsheet writer exists)
%     'cases'      arm | Ma | THD49_pct | theta01_deg..theta15_deg
%                  (first sheet: input of ideal_switch_validation)
%     'predicted'  analytic metrics of every case
%     'notes'      provenance

D.budget   = 20000;
D.seed     = 1;
D.Ma_list  = [1.20 1.00 0.80 0.60 0.40];
D.arms     = {'gn_hybrid','gn_multistart','sqp_multistart','gwo_sqp_box','multistart_tb','fssa'};
D.out      = 'validation_cases.xlsx';
D.N        = 15;
D.only_accepted = true;
D.verbose  = true;

if nargin < 2, CFG = struct(); end
fn = fieldnames(CFG);  C = D;
for i = 1:numel(fn), C.(fn{i}) = CFG.(fn{i}); end
if ischar(C.arms), C.arms = {C.arms}; end

% ---------------- read the benchmark CSV -----------------------------------
fid = fopen(csvfile, 'r');
if fid < 0, error('export_validation_cases:open','cannot open %s', csvfile); end
H  = strsplit(strtrim(fgetl(fid)), ',');
nC = numel(H);
Craw = textscan(fid, [repmat('%s',1,4) repmat('%f',1,nC-4)], ...
                'Delimiter', ',', 'EmptyValue', NaN, 'ReturnOnError', false);
fclose(fid);

col = @(nm) find(strcmp(H, nm), 1);
need = {'arm','nfe_target','seed','Ma','accepted','THD49_pct'};
for i = 1:numel(need)
    if isempty(col(need{i}))
        error('export_validation_cases:schema','%s has no column "%s"', csvfile, need{i});
    end
end
num = @(nm) Craw{col(nm)};
armcol = Craw{col('arm')};
thcol  = arrayfun(@(k) col(sprintf('theta%02d_deg',k)), 1:C.N);
if any(thcol == 0) || numel(thcol) ~= C.N
    error('export_validation_cases:angles','%s does not carry theta01_deg..theta%02d_deg', csvfile, C.N);
end
TH = zeros(numel(armcol), C.N);
for k = 1:C.N, TH(:,k) = Craw{thcol(k)}; end

B  = num('nfe_target');  SD = num('seed');  MA = num('Ma');
AC = num('accepted');    T49 = num('THD49_pct');

if C.verbose
    fprintf('\n=========== export_validation_cases: %s ===========\n', csvfile);
    fprintf('filter : nfe_target = %g, seed = %d\n', C.budget, C.seed);
    fprintf('arms   : %s\n', strjoin(C.arms, ', '));
    fprintf('Ma     : %s\n\n', mat2str(C.Ma_list(:).',3));
end

% ---------------- select ---------------------------------------------------
cases = {};  pred = {};  miss = {};
base  = (abs(B - C.budget) < 1e-9) & (abs(SD - C.seed) < 1e-9);
for ai = 1:numel(C.arms)
    key = C.arms{ai};
    nm  = paper_label(key);
    for q = 1:numel(C.Ma_list)
        r = find(base & strcmp(armcol, key) & abs(MA - C.Ma_list(q)) < 1e-9, 1);
        if isempty(r)
            miss{end+1,1} = sprintf('%-16s Ma=%.2f  not in file at budget %g seed %d', ...
                                    key, C.Ma_list(q), C.budget, C.seed); %#ok<AGROW>
            continue;
        end
        if C.only_accepted && ~(AC(r) == 1)
            miss{end+1,1} = sprintf('%-16s Ma=%.2f  not accepted -- not exported', ...
                                    key, C.Ma_list(q)); %#ok<AGROW>
            continue;
        end
        cases(end+1,:) = [{nm, MA(r), T49(r)}, num2cell(TH(r,:))]; %#ok<AGROW>
        row = {nm, key, MA(r), AC(r), T49(r)};
        for f = {'THD99_pct','WTHD49_pct','WTHD999_pct','worst_harm_above49_pct', ...
                 'worst_harm_order','b3_pct','MAXQ','b1','b1_err','viol', ...
                 'nfe_achieved','time_s'}
            c = col(f{1});
            if isempty(c), row{end+1} = NaN; else, row{end+1} = Craw{c}(r); end %#ok<AGROW>
        end
        pred(end+1,:) = row; %#ok<AGROW>
    end
end

if isempty(cases)
    error('export_validation_cases:empty', ['nothing matched; check that %s contains ' ...
          'nfe_target = %g and seed = %d.'], csvfile, C.budget, C.seed);
end
for i = 1:numel(miss), fprintf('  -- %s\n', miss{i}); end

% ---------------- write ----------------------------------------------------
case_hdr = [{'arm','Ma','THD49_pct'}, ...
            arrayfun(@(k) sprintf('theta%02d_deg',k), 1:C.N, 'UniformOutput', false)];
pred_hdr = {'method','arm_key','Ma','accepted','THD49_pct','THD99_pct','WTHD49_pct', ...
            'WTHD999_pct','worst_harm_above49_pct','worst_harm_order','b3_pct','MAXQ', ...
            'b1','b1_err','viol','nfe_achieved','time_s'};
N = {};
N{end+1,1} = 'Ideal-switch Simulink check (Section 6.2): input cases.';
N{end+1,1} = sprintf('Extracted from : %s', csvfile);
N{end+1,1} = sprintf('Filter         : nfe_target = %g, seed = %d', C.budget, C.seed);
N{end+1,1} = sprintf('Arms           : %s', strjoin(C.arms, ', '));
N{end+1,1} = sprintf('Ma             : %s', mat2str(C.Ma_list(:).',3));
N{end+1,1} = sprintf('Cases          : %d', size(cases,1));
N{end+1,1} = 'Angles are the benchmark rows themselves (not re-solved).';
N{end+1,1} = 'Acceptance: |b1-Ma| < 1e-8 and violation < 1e-9 (evaluate_candidates).';
N{end+1,1} = 'Counted band: odd non-triplen orders 5..49 (16 orders).';
for i = 1:numel(miss), N{end+1,1} = ['SKIPPED: ' miss{i}]; end %#ok<AGROW>

[ok, what] = write_book(C.out, {'cases',     case_hdr, cases; ...
                                'predicted', pred_hdr, pred; ...
                                'notes',     {'note'}, N});

% ---------------- report ---------------------------------------------------
if C.verbose
    u = unique(cases(:,1), 'stable');
    fprintf('\n--- analytic THD49 (%%) of the exported cases ---\n%-8s', 'Ma');
    for i = 1:numel(u), fprintf('%12s', u{i}); end
    fprintf('\n');
    for q = 1:numel(C.Ma_list)
        fprintf('%-8.2f', C.Ma_list(q));
        for i = 1:numel(u)
            r = find(strcmp(cases(:,1), u{i}) & ...
                     abs(cell2mat(cases(:,2)) - C.Ma_list(q)) < 1e-9, 1);
            if isempty(r), fprintf('%12s','--'); else, fprintf('%12.4f', cases{r,3}); end
        end
        fprintf('\n');
    end
    fprintf('\nwrote %s  (%d cases)%s\n', what, size(cases,1), ...
            tern(ok,'',' -- no spreadsheet writer on this host'));
end

OUT = struct('file', what, 'cases', {cases}, 'predicted', {pred}, 'skipped', {miss});
end

% =========================================================================
function nm = paper_label(key)
%PAPER_LABEL  Arm key -> name used in the paper.
switch lower(key)
  case 'gn_hybrid',      nm = 'PLM-C';
  case 'gn_multistart',  nm = 'PLM-MS';
  case 'sqp_hybrid',     nm = 'SQP-C';
  case 'sqp_multistart', nm = 'SQP-MS';
  case 'fssa',           nm = 'FSSA';
  case 'gwo_sqp_box',    nm = 'GWO-SQP';
  case 'multistart_tb',  nm = 'MS-TB';
  otherwise,             nm = key;
end
end

function [ok, what] = write_book(fx, sheets)
%WRITE_BOOK  sheets = {name, header, body} per row; .xlsx first, CSV per
%   sheet if the host has no spreadsheet writer.
ok = false;  what = fx;
[pth, base] = fileparts(fx);
have = (exist('writecell','file') == 2) || (exist('writetable','file') == 2) || ...
       (exist('xlswrite','file') == 2);
if have
    try
        if exist(fx,'file') == 2, delete(fx); end
        for i = 1:size(sheets,1)
            M = blankify([sheets{i,2}(:).'; sheets{i,3}]);
            if exist('writecell','file') == 2
                writecell(M, fx, 'Sheet', sheets{i,1});
            elseif exist('writetable','file') == 2
                vn = cellfun(@valid_name, M(1,:), 'UniformOutput', false);
                writetable(cell2table(M(2:end,:), 'VariableNames', vn), fx, ...
                           'Sheet', sheets{i,1});
            else
                xlswrite(fx, M, sheets{i,1});
            end
        end
        ok = true;
        what = fx;
        return;
    catch e
        fprintf('  (xlsx failed: %s -- writing CSV sheets instead)\n', e.message);
    end
end
names = cell(size(sheets,1),1);
for i = 1:size(sheets,1)
    fn = fullfile(pth, sprintf('%s_%s.csv', base, sheets{i,1}));
    wcell(fn, [sheets{i,2}(:).'; sheets{i,3}]);
    names{i} = fn;
end
what = names{1};
end

function M = blankify(M)
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
            if any(v == ','), fprintf(fid, '"%s"', v); else, fprintf(fid, '%s', v); end
        elseif isempty(v) || (isnumeric(v) && ~isfinite(v))
            % blank
        else
            fprintf(fid, '%.10g', double(v));
        end
    end
    fprintf(fid, '\n');
end
fclose(fid);
end

function s = valid_name(s)
s = regexprep(strtrim(s), '\W', '_');
if isempty(s) || ~isletter(s(1)), s = ['v' s]; end
end

function s = tern(c,a,b), if c, s = a; else, s = b; end, end
