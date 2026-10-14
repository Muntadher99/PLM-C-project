function run_main_benchmark(outdir)
%RUN_MAIN_BENCHMARK  Rerun the matched-budget benchmark of the six methods
%   (paper, Sections 4, 5.1-5.3 and 5.6; Tables 2, 3 and 4; Reviewer 1
%   comments 7, 8, 9 and 11, Reviewer 2 comment 6).
%
%   run_main_benchmark()          % writes to 3_benchmark_protocol/reproduced
%   run_main_benchmark(outdir)
%
%   Configuration of the paper (run_benchmark defaults): 21 operating points
%   Ma = 1.20:-0.05:0.20, budgets 500 1000 2000 4000 8000 20000 NFE per
%   point, seeds 1-10, N = 15 angles, orders up to 49, minimum gap 0.25 deg.
%
%   Outputs (each with .xlsx and _summary.xlsx)
%     bench_main.csv     PLM-C, PLM-MS, SQP-MS, MS-TB, FSSA
%     bench_gwo_sqp.csv  GWO-SQP
%   Stored copies of both are in 3_benchmark_protocol/results.
%
%   Requires the Optimization Toolbox (SQP-MS, GWO-SQP) and the Global
%   Optimization Toolbox (MS-TB).

if nargin < 1 || isempty(outdir)
    outdir = fullfile(fileparts(mfilename('fullpath')), 'reproduced');
end
if ~exist(outdir, 'dir'), mkdir(outdir); end

run_benchmark(struct('arms', {{'gn_hybrid','gn_multistart','sqp_multistart','multistart_tb','fssa'}}, ...
                     'out', fullfile(outdir, 'bench_main.csv')));
run_benchmark(struct('arms', {{'gwo_sqp_box'}}, ...
                     'out', fullfile(outdir, 'bench_gwo_sqp.csv')));
end
