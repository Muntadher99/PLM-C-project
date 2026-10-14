function run_sensitivity(main_csv, outdir)
%RUN_SENSITIVITY  Rerun the implementation-sensitivity study of Section 6.1
%   (Reviewer 1 comment 4).
%
%   run_sensitivity()                   % writes to 9_sensitivity/reproduced
%   run_sensitivity(main_csv, outdir)
%
%   main_csv : main benchmark CSV whose PLM-C angles at 20000 NFE are
%              perturbed [3_benchmark_protocol/results/bench_main.csv]
%
%   Outputs
%     sensitivity.xlsx and sensitivity_<study>.csv for the studies deadtime,
%       npripple, quant, matol (PLM-C angles of seeds 1-10) and mingap
%       (re-optimised once); printed by sensitivity_summary
%     bench_mingap050.csv (+ .xlsx, _summary.xlsx): the main benchmark of the
%       five methods repeated with a minimum gap of 0.50 deg

here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(main_csv)
    main_csv = fullfile(fileparts(here), '3_benchmark_protocol', 'results', 'bench_main.csv');
end
if nargin < 2 || isempty(outdir), outdir = fullfile(here, 'reproduced'); end
if ~exist(outdir, 'dir'), mkdir(outdir); end

base = fullfile(outdir, 'sensitivity');
sensitivity_study(struct('csv', main_csv, 'seed', 1:10, 'out', base));
sensitivity_summary(base);

run_benchmark(struct('arms', {{'gn_hybrid','gn_multistart','sqp_multistart','multistart_tb','fssa'}}, ...
                     'minGap_deg', 0.50, 'out', fullfile(outdir, 'bench_mingap050.csv')));
end
