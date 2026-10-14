function run_ablations(main_csv, outdir)
%RUN_ABLATIONS  Rerun the component ablations (paper, Sections 2.2 and 5.4;
%   Tables 5 and 6; Reviewer 1 comment 6, Reviewer 2 comment 1).
%
%   run_ablations()                    % writes to 4_ablations/reproduced
%   run_ablations(main_csv, outdir)
%
%   main_csv : main benchmark CSV whose 20000-NFE solutions are analysed by
%              gauss_newton_rate [3_benchmark_protocol/results/bench_main.csv]
%
%   Outputs
%     bench_forward_continuation.csv   forward-only continuation of [2] with the
%                                      PLM corrector (Table 5, row 2),
%                                      budgets 1000 2000 8000 20000
%     bench_sqp_continuation.csv       SQP-C: fmincon SQP corrector inside the
%                                      bidirectional continuation (Table 5,
%                                      row 5), budgets 1000-20000
%     corrector_ablation_production.csv  same 200 starts per point, production
%                                      iteration caps (Table 6: acceptance,
%                                      distortion, NFE per start, batching)
%     corrector_ablation_converged.csv   the same starts at ten times the caps
%                                      (Table 6: share reaching a KKT point)
%     gauss_newton_rate.csv            local rate rho of undamped Gauss-Newton
%                                      at the 20000-NFE solutions (Section 2.2)
%   Table 5 also uses the PLM-MS, PLM-C and SQP-MS rows of bench_main.csv.

here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(main_csv)
    main_csv = fullfile(fileparts(here), '3_benchmark_protocol', 'results', 'bench_main.csv');
end
if nargin < 2 || isempty(outdir), outdir = fullfile(here, 'reproduced'); end
if ~exist(outdir, 'dir'), mkdir(outdir); end

% ---- Table 5: start strategy and corrector --------------------------------
run_benchmark(struct('arms', {{'gn_hybrid_zeroth_fwd'}}, 'budgets', [1000 2000 8000 20000], ...
                     'out', fullfile(outdir, 'bench_forward_continuation.csv')));
run_benchmark(struct('arms', {{'sqp_hybrid'}}, 'budgets', [1000 2000 4000 8000 20000], ...
                     'out', fullfile(outdir, 'bench_sqp_continuation.csv')));

% ---- Table 6: the local corrector at matched starts -----------------------
corrector_ablation(struct('out', fullfile(outdir, 'corrector_ablation_production.csv')));
corrector_ablation(struct('iter_scale', 10, 'variants', {{'proposed','gn','lm1','theta'}}, ...
                          'out', fullfile(outdir, 'corrector_ablation_converged.csv')));

% ---- Section 2.2: measured local rate of undamped Gauss-Newton ------------
gauss_newton_rate(main_csv, {'gn_hybrid','gn_multistart','sqp_multistart'}, 20000, ...
                  fullfile(outdir, 'gauss_newton_rate.csv'));
end
