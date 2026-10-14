function out = merge_benchmark_csv(files, outfile, keep)
%MERGE_BENCHMARK_CSV  Combine run_benchmark CSVs; a later file supersedes
%   every row of an arm it contains.
%
%   merge_benchmark_csv({'a.csv','b.csv'}, 'merged.csv')
%   merge_benchmark_csv(files, outfile, {'gn_hybrid','fssa'})   % keep only these arms
%
%   If a later file contains any row of arm X, all rows of X in the earlier
%   files are dropped, so one arm never mixes two runs. All files must have
%   identical headers. The merged CSV is also written as .xlsx when possible.

if ischar(files), files = {files}; end
if nargin < 2 || isempty(outfile), outfile = 'bench_merged.csv'; end
if nargin < 3, keep = {}; end
if ischar(keep) && ~isempty(keep), keep = {keep}; end

% ---- read every file, header-checked --------------------------------------
H0 = [];  raw = cell(numel(files),1);
for i = 1:numel(files)
    [H, L] = read_lines(files{i});
    if isempty(H0)
        H0 = H;
    elseif ~isequal(H, H0)
        error('merge_benchmark_csv:schema', ...
              '%s has a different header from %s', files{i}, files{1});
    end
    raw{i} = L;
end

% ---- which arms does each file supply? -----------------------------------
arms = cell(numel(files),1);
for i = 1:numel(files)
    a = cell(numel(raw{i}),1);
    for r = 1:numel(raw{i}), a{r} = first_field(raw{i}{r}); end
    arms{i} = a;
end

owner = containers.Map();     % arm -> index of the LAST file that supplies it
for i = 1:numel(files)
    u = unique(arms{i});
    for j = 1:numel(u), owner(u{j}) = i; end
end

% ---- write ----------------------------------------------------------------
fid = fopen(outfile, 'w');
fprintf(fid, '%s\n', H0);
kept = containers.Map();  nrow = 0;
for i = 1:numel(files)
    for r = 1:numel(raw{i})
        a = arms{i}{r};
        if owner(a) ~= i, continue; end                       % superseded
        if ~isempty(keep) && ~any(strcmp(keep, a)), continue; end
        fprintf(fid, '%s\n', raw{i}{r});
        nrow = nrow + 1;
        if kept.isKey(a), kept(a) = kept(a) + 1; else, kept(a) = 1; end
    end
end
fclose(fid);

% ---- report ---------------------------------------------------------------
fprintf('\n=================== merge_benchmark_csv ===================\n');
for i = 1:numel(files)
    u = unique(arms{i});
    own = u(cellfun(@(a) owner(a) == i, u));
    sup = u(cellfun(@(a) owner(a) ~= i, u));
    fprintf('  %s\n', files{i});
    if ~isempty(own), fprintf('      supplies  : %s\n', strjoin(own(:).', ', ')); end
    if ~isempty(sup), fprintf('      superseded: %s\n', strjoin(sup(:).', ', ')); end
end
if ~isempty(keep)
    missing = keep(~cellfun(@(a) kept.isKey(a), keep));
    if ~isempty(missing)
        fprintf('  requested but not present in any file: %s\n', strjoin(missing, ', '));
    end
end
k = kept.keys();
fprintf('\n  wrote %s  (%d rows)\n', outfile, nrow);
for i = 1:numel(k), fprintf('      %-20s %6d rows\n', k{i}, kept(k{i})); end

try, csv_to_xlsx(outfile, 'runs'); catch, end
out = outfile;
end

% =========================================================================
function [H, L] = read_lines(fn)
fid = fopen(fn, 'r');
if fid < 0, error('merge_benchmark_csv:open','cannot open %s', fn); end
H = strtrim(fgetl(fid));
L = {};
while true
    s = fgetl(fid);
    if ~ischar(s), break; end
    s = strtrim(s);
    if ~isempty(s), L{end+1,1} = s; end %#ok<AGROW>
end
fclose(fid);
end

function a = first_field(s)
i = find(s == ',', 1);
if isempty(i), a = s; else, a = s(1:i-1); end
end
