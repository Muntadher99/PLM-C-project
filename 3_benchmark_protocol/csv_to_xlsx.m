function fx = csv_to_xlsx(csvfile, sheetname)
%CSV_TO_XLSX  Write a run_benchmark CSV as an .xlsx workbook (one sheet).
%
%   csv_to_xlsx('bench_main.csv')
%   csv_to_xlsx(csvfile, 'runs')
%
%   The first four columns are text labels, the rest numeric; NaN is written
%   as an empty cell. Spreadsheet programs that split .csv files on a regional
%   list separator other than ',' read the .xlsx unchanged.

if nargin < 2 || isempty(sheetname), sheetname = 'runs'; end
[pth, base] = fileparts(csvfile);
fx = fullfile(pth, [base '.xlsx']);

% ---- read: first four columns are labels, the rest are numeric ------------
fid = fopen(csvfile, 'r');
if fid < 0, error('csv_to_xlsx:open','cannot open %s', csvfile); end
H = strsplit(strtrim(fgetl(fid)), ',');
nC = numel(H);
fmt = [repmat('%s', 1, 4) repmat('%f', 1, nC-4)];
Craw = textscan(fid, fmt, 'Delimiter', ',', 'EmptyValue', NaN, 'ReturnOnError', false);
fclose(fid);

nR = numel(Craw{1});
if nR == 0, error('csv_to_xlsx:empty','%s has a header but no rows', csvfile); end

M = cell(nR+1, nC);
M(1,:) = H;
for c = 1:nC
    v = Craw{c};
    if c <= 4
        M(2:end,c) = v;
    else
        M(2:end,c) = num2cell(v);
    end
end
for i = 1:numel(M)
    v = M{i};
    if isnumeric(v) && (isempty(v) || ~all(isfinite(v(:)))), M{i} = ''; end
end

% ---- write ---------------------------------------------------------------
if nR + 1 > 1048576
    error('csv_to_xlsx:rows','%d rows exceeds the xlsx limit', nR);
end
try
    if exist(fx,'file') == 2, delete(fx); end
    if exist('writecell','file') == 2
        writecell(M, fx, 'Sheet', sheetname);
    elseif exist('writetable','file') == 2
        vn = cellfun(@(s) valid_name(s), M(1,:), 'UniformOutput', false);
        writetable(cell2table(M(2:end,:), 'VariableNames', vn), fx, 'Sheet', sheetname);
    elseif exist('xlswrite','file') == 2
        xlswrite(fx, M, sheetname);
    else
        fx = '';
        fprintf('  (no spreadsheet writer on this host -- %s left as the record)\n', csvfile);
        return;
    end
    fprintf('  wrote %s  (%d rows x %d columns, sheet ''%s'')\n', fx, nR, nC, sheetname);
catch e
    fprintf('  (xlsx skipped: %s -- %s is unchanged)\n', e.message, csvfile);
    fx = '';
end
end

function s = valid_name(s)
s = regexprep(strtrim(s), '\W', '_');
if isempty(s) || ~isletter(s(1)), s = ['v' s]; end
end
