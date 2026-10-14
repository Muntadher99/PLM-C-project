function run_angle_count(outdir)
%RUN_ANGLE_COUNT  Rerun the benchmark at eleven and nineteen angles (paper,
%   Section 5.5; Reviewer 1 comment 1, Reviewer 2 comment 7).
%
%   run_angle_count()            % writes to 5_angle_count/reproduced
%   run_angle_count(outdir)
%
%   Same protocol as the main benchmark (21 points, seeds 1-10, orders up to
%   49, minimum gap 0.25 deg) for PLM-C, PLM-MS and SQP-MS at 2000 and 20000
%   NFE per point. With 16 counted orders, N = 11 is more over-determined
%   than N = 15, while N = 19 admits exact solutions.
%
%   Outputs (each with .xlsx and _summary.xlsx)
%     bench_N11.csv, bench_N19.csv

if nargin < 1 || isempty(outdir)
    outdir = fullfile(fileparts(mfilename('fullpath')), 'reproduced');
end
if ~exist(outdir, 'dir'), mkdir(outdir); end

arms = {'gn_hybrid','gn_multistart','sqp_multistart'};
for N = [11 19]
    run_benchmark(struct('arms', {arms}, 'budgets', [2000 20000], 'N', N, ...
                         'out', fullfile(outdir, sprintf('bench_N%d.csv', N))));
end
end
