function D = read_benchmark_csv(csvfile)
%READ_BENCHMARK_CSV  Read a run_benchmark CSV into plain arrays (MATLAB or Octave).
%
%   D = read_benchmark_csv('bench_main.csv')
%
%   D fields: arm (cellstr), nfe_target, nfe_achieved, seed, Ma, accepted,
%   THD49 (%), MAXQ, compliant, time_s, theta (degrees, rows x N), N.
%   Lines with the wrong number of fields and rows with an empty arm field
%   are skipped and reported.

txt = fileread(csvfile);
L = regexp(txt, '\r?\n', 'split');
L = L(~cellfun(@isempty, strtrim(L)));
H = regexp(strtrim(L{1}), ',', 'split');
nc = numel(H);
nr = numel(L) - 1;
C = repmat({cell(nr,1)}, 1, nc);
bad = false(nr,1);
for r = 1:nr
    f = regexp(L{r+1}, ',', 'split');
    if numel(f) ~= nc, bad(r) = true; continue; end   % malformed line: skip
    for c = 1:nc, C{c}{r} = f{c}; end
end
if any(bad)
    fprintf('read_benchmark_csv: skipped %d malformed line(s) in %s\n', sum(bad), csvfile);
    for c = 1:nc, C{c} = C{c}(~bad); end
end
col = @(name) C{find(strcmp(H, name), 1)};
num = @(name) str2double(col(name));

arm = col('arm');
keep = ~cellfun(@isempty, arm);
D.arm          = arm(keep);
D.nfe_target   = pick(num('nfe_target'), keep);
D.nfe_achieved = pick(num('nfe_achieved'), keep);
D.seed         = pick(num('seed'), keep);
D.Ma           = pick(num('Ma'), keep);
D.accepted     = pick(num('accepted'), keep);
D.THD49        = pick(num('THD49_pct'), keep);
D.MAXQ         = pick(num('MAXQ'), keep);
D.compliant    = pick(num('compliant'), keep);
D.time_s       = pick(num('time_s'), keep);
tc = find(strncmp(H, 'theta', 5));
D.N = numel(tc);
D.theta = zeros(sum(keep), D.N);
for j = 1:D.N
    v = str2double(C{tc(j)});
    D.theta(:,j) = v(keep);
end
end

function y = pick(x, k)
y = x(k);  y = y(:);
end
