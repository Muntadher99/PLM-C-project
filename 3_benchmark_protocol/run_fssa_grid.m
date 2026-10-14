function run_fssa_grid(outdir)
%RUN_FSSA_GRID  Rerun FSSA with the weighted grid-code objective of [2]
%   (paper, Section 5.6; response to Reviewer 1, comment 7).
%
%   run_fssa_grid()          % writes to 3_benchmark_protocol/reproduced
%   run_fssa_grid(outdir)
%
%   Same protocol as the main benchmark (21 operating points
%   Ma = 1.20:-0.05:0.20, seeds 1-10, N = 15 angles, orders up to 49, minimum
%   gap 0.25 deg) at the largest budget, 20000 NFE per point. The arm
%   fssa_grid is the FSSA arm of the main benchmark with OBJ = 'grid' in place
%   of OBJ = 'ls' (arm registry of run_benchmark); cold restarts are off, as
%   in [2].
%
%   Output (with .xlsx and _summary.xlsx)
%     bench_fssa_grid.csv
%   The stored copy is in 3_benchmark_protocol/results.

if nargin < 1 || isempty(outdir)
    outdir = fullfile(fileparts(mfilename('fullpath')), 'reproduced');
end
if ~exist(outdir, 'dir'), mkdir(outdir); end

run_benchmark(struct('arms', {{'fssa_grid'}}, 'budgets', 20000, ...
                     'out', fullfile(outdir, 'bench_fssa_grid.csv')));
end
