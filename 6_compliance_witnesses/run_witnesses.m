function run_witnesses(main_csv, outdir)
%RUN_WITNESSES  Rerun the minimax compliance stage and rebuild the audited
%   witness tables (paper, Sections 3.4 and 5.7; Table 7; Reviewer 1
%   comments 2, 3 and 10, Reviewer 2 comments 2 and 5).
%
%   run_witnesses()                    % writes to 6_compliance_witnesses/reproduced
%   run_witnesses(main_csv, outdir)
%
%   main_csv : main benchmark CSV; its lowest-THD49 PLM-C solution at 20000
%              NFE is the distortion-minimising design of each point
%              [3_benchmark_protocol/results/bench_main.csv]
%
%   Outputs
%     minimax_search.csv (+ _summary.xlsx)  minimax_witness_search, seeds 1-5,
%                               20000 NFE per point and mode
%     witnesses_all_audit.csv   every minimax angle set with b1 error, gaps,
%                               all sixteen q_n and q_max
%     witnesses_selected.csv    certificate (lowest Q^) and compliant design
%                               (lowest THD49 with q_max <= 1) per point
%     thdmin_selected.csv       distortion-minimising design per point
%     unified_table.csv         Table 7
%     nonideal_cases.csv        the ten angle sets of the device-level study

here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(main_csv)
    main_csv = fullfile(fileparts(here), '3_benchmark_protocol', 'results', 'bench_main.csv');
end
if nargin < 2 || isempty(outdir), outdir = fullfile(here, 'reproduced'); end
if ~exist(outdir, 'dir'), mkdir(outdir); end

mm = fullfile(outdir, 'minimax_search.csv');
minimax_witness_search(struct('seeds', 1:5, 'NFE_TARGET', 20000, 'out', mm));
build_witness_tables(mm, main_csv, 20000, outdir);
end
