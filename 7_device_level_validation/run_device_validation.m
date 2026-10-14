function run_device_validation(cases_csv, outdir, ideal_check, bench_csvs)
%RUN_DEVICE_VALIDATION  Rerun the Simulink validation of Section 6.2
%   (Table 8; Reviewer 1 comment 4, Reviewer 2 comments 3 and 7).
%
%   run_device_validation()                            % both studies
%   run_device_validation(cases_csv, outdir, false)    % device-level study only
%   run_device_validation(cases_csv, outdir, true, {main_csv, gwo_csv})
%
%   cases_csv   : the ten angle sets of the device-level study, written by
%                 build_witness_tables [7_device_level_validation/results/
%                 nonideal_cases.csv]
%   outdir      : output folder [7_device_level_validation/reproduced]
%   ideal_check : also run the ideal-switch check [true]
%   bench_csvs  : benchmark CSVs holding the six methods [the stored
%                 bench_main.csv and bench_gwo_sqp.csv of 3_benchmark_protocol]
%
%   1. Device-level study (Table 8): nonideal_simulink_runs on NPC_SHM.slx,
%      six cases x ten angle sets -> nonideal_results.csv
%   2. Ideal-switch check (first sentence of Section 6.2): the six methods at
%      seed 1, 20000 NFE, Ma = 1.20 1.00 0.80 0.60 0.40 (thirty angle sets)
%      are taken from the stored benchmark rows (export_validation_cases) and
%      simulated in NPC_SHM_ideal.slx (ideal_switch_validation, run in the
%      base workspace because parsim transfers it to the workers)
%      -> ideal_switch_check/reproduced/validation_cases.xlsx and
%         validation_cases_sim.xlsx
%
%   Requires Simulink and Simscape Electrical; the ideal-switch check also
%   needs the Parallel Computing Toolbox (parsim). MATLAB R2024b.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
if nargin < 1 || isempty(cases_csv), cases_csv = fullfile(here, 'results', 'nonideal_cases.csv'); end
if nargin < 2 || isempty(outdir)
    outdir = fullfile(here, 'reproduced');
    idir   = fullfile(here, 'ideal_switch_check', 'reproduced');
else
    idir   = fullfile(outdir, 'ideal_switch_check');
end
if nargin < 3 || isempty(ideal_check), ideal_check = true; end
if ~exist(outdir, 'dir'), mkdir(outdir); end

% ---- 1. device-level study (Table 8) --------------------------------------
nonideal_simulink_runs(cases_csv, struct('out', fullfile(outdir, 'nonideal_results.csv')));

% ---- 2. ideal-switch check -------------------------------------------------
if ~ideal_check, return; end
if ~exist(idir, 'dir'), mkdir(idir); end
if nargin < 4 || isempty(bench_csvs)
    res = fullfile(root, '3_benchmark_protocol', 'results');
    bench_csvs = {fullfile(res, 'bench_main.csv'), fullfile(res, 'bench_gwo_sqp.csv')};
end
merged = fullfile(idir, 'bench_all.csv');
merge_benchmark_csv(bench_csvs, merged);
export_validation_cases(merged, struct('out', fullfile(idir, 'validation_cases.xlsx')));
old = cd(idir);
restore = onCleanup(@() cd(old));   %#ok<NASGU>
evalin('base', 'ideal_switch_validation');
end
