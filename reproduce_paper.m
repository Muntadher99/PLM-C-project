function reproduce_paper(steps)
%REPRODUCE_PAPER  Single entry point of the PLM-C reproducibility package
%   (paper IJIES 20266308).
%
%   reproduce_paper                     % print Tables 3-6, 8, 9 from the stored results
%   reproduce_paper('all')              % rerun every experiment, then print the
%                                       % tables from the reproduced files
%   reproduce_paper({'benchmark','witnesses'})   % rerun selected steps only
%
%   Steps, run in this order (logged solver time of the published runs):
%     'benchmark'    3_benchmark_protocol       run_main_benchmark       ~1.4 h
%     'ablations'    4_ablations                run_ablations            ~0.5 h
%     'angle_count'  5_angle_count              run_angle_count          ~0.6 h
%     'witnesses'    6_compliance_witnesses     run_witnesses            < 0.1 h
%     'device'       7_device_level_validation  run_device_validation    ~0.7 h + ideal check
%     'loss'         8_loss_analysis            run_loss_analysis        ~1 min
%     'sensitivity'  9_sensitivity              run_sensitivity          ~1.3 h
%   Each step writes to <folder>/reproduced/ under the file names of
%   <folder>/results/. A step reads the reproduced output of an earlier step
%   when that step is part of the same call, and the stored results
%   otherwise (ablations, witnesses and sensitivity read bench_main.csv;
%   device reads nonideal_cases.csv; loss reads the witness tables).
%
%   'device' needs Simulink and Simscape Electrical (and the Parallel
%   Computing Toolbox for the ideal-switch check); every other step needs
%   MATLAB with the Optimization and Global Optimization Toolboxes.

root = fileparts(mfilename('fullpath'));
setup_paths();
if nargin < 1 || isempty(steps)
    paper_tables();
    return;
end

order = {'benchmark','ablations','angle_count','witnesses','device','loss','sensitivity'};
if ischar(steps)
    if strcmpi(steps, 'all'), steps = order; else, steps = {steps}; end
end
unknown = setdiff(steps, order);
if ~isempty(unknown)
    error('reproduce_paper:step', 'unknown step(s): %s (valid: %s)', ...
          strjoin(unknown, ', '), strjoin(order, ', '));
end

res = @(folder, f) fullfile(root, folder, 'results', f);
rep = @(folder, f) fullfile(root, folder, 'reproduced', f);
main_csv  = res('3_benchmark_protocol', 'bench_main.csv');
gwo_csv   = res('3_benchmark_protocol', 'bench_gwo_sqp.csv');
wit_dir   = fullfile(root, '6_compliance_witnesses', 'results');
cases_csv = res('7_device_level_validation', 'nonideal_cases.csv');

t_all = tic;
for k = 1:numel(order)
    step = order{k};
    if ~any(strcmp(steps, step)), continue; end
    fprintf('\n######## reproduce_paper: %s ########\n', step);
    t0 = tic;
    switch step
        case 'benchmark'
            run_main_benchmark();
            main_csv = rep('3_benchmark_protocol', 'bench_main.csv');
            gwo_csv  = rep('3_benchmark_protocol', 'bench_gwo_sqp.csv');
        case 'ablations'
            run_ablations(main_csv);
        case 'angle_count'
            run_angle_count();
        case 'witnesses'
            run_witnesses(main_csv);
            wit_dir   = fullfile(root, '6_compliance_witnesses', 'reproduced');
            cases_csv = fullfile(wit_dir, 'nonideal_cases.csv');
        case 'device'
            run_device_validation(cases_csv, [], true, {main_csv, gwo_csv});
        case 'loss'
            run_loss_analysis(wit_dir);
        case 'sensitivity'
            run_sensitivity(main_csv);
    end
    fprintf('######## %s done in %.1f min ########\n', step, toc(t0)/60);
end
fprintf('\nall requested steps done in %.1f h\n', toc(t_all)/3600);
paper_tables('reproduced');
end
