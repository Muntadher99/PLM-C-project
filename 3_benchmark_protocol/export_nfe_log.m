function outfile = export_nfe_log(which, outfile)
%EXPORT_NFE_LOG  Achieved-NFE log of every benchmark run (Reviewer 1,
%   comment 8; paper, Section 4).
%
%   export_nfe_log()                 % stored results/ files
%                                    %   -> 3_benchmark_protocol/results/nfe_log.csv
%   export_nfe_log('reproduced')     % reproduced/ files of a rerun
%                                    %   -> 3_benchmark_protocol/reproduced/nfe_log.csv
%   export_nfe_log(which, outfile)   % explicit output file
%
%   Reads the nfe_target and nfe_achieved columns of every benchmark CSV of
%   the package and writes one row per sweep (file, method, budget, seed).
%   Nothing is rerun or recomputed: nfe_per_point is the nfe_achieved value
%   stored in the benchmark CSV, which run_benchmark writes on each of the
%   rows of a sweep.
%
%   Columns of nfe_log.csv
%     source           benchmark CSV the row is read from
%     arm, method      arm key and the name used in the paper
%     nfe_target       evaluation budget per operating point
%     seed             seed of the sweep
%     n_points         operating points in the sweep
%     nfe_per_point    achieved evaluations per operating point (nfe_achieved)
%     nfe_sweep_total  n_points * nfe_per_point when nfe_per_point is a mean
%                      over the sweep; empty otherwise
%     nfe_ratio        nfe_per_point / nfe_target
%     count_basis      what nfe_per_point measures:
%                      'mean_over_points'  PLM-C, PLM-MS, SQP-MS, SQP-C,
%                                          forward continuation, FSSA: the
%                                          evaluations of the sweep divided by
%                                          its operating points
%                      'last_point'        GWO-SQP, MS-TB: the evaluations
%                                          charged at the last operating point
%                                          of the sweep (their counter is reset
%                                          at every point, see run_benchmark)
%
%   One NFE is one evaluation of Eq. (2) for one angle vector, Jacobian
%   included; the external restoration of evaluate_candidates is not charged.
%   Runs in MATLAB and GNU Octave.

if nargin < 1 || isempty(which), which = 'results'; end
which = lower(char(which));
if ~any(strcmp(which, {'results', 'reproduced'}))
    error('export_nfe_log:which', 'which must be ''results'' or ''reproduced''');
end
root = fileparts(fileparts(mfilename('fullpath')));
if nargin < 2 || isempty(outfile)
    outdir = fullfile(root, '3_benchmark_protocol', which);
    if ~exist(outdir, 'dir'), mkdir(outdir); end
    outfile = fullfile(outdir, 'nfe_log.csv');
end

% benchmark CSVs of the package, in the order of the paper
files = { ...
    '3_benchmark_protocol', 'bench_main.csv'
    '3_benchmark_protocol', 'bench_gwo_sqp.csv'
    '3_benchmark_protocol', 'bench_fssa_grid.csv'
    '4_ablations',          'bench_forward_continuation.csv'
    '4_ablations',          'bench_sqp_continuation.csv'
    '5_angle_count',        'bench_N11.csv'
    '5_angle_count',        'bench_N19.csv'
    '9_sensitivity',        'bench_mingap050.csv'};

keys  = {'gn_hybrid','gn_multistart','sqp_multistart','sqp_hybrid', ...
         'gn_hybrid_zeroth_fwd','gwo_sqp_box','multistart_tb','fssa','fssa_grid'};
names = {'PLM-C','PLM-MS','SQP-MS','SQP-C', ...
         'Forward continuation (Table 5 row 2)','GWO-SQP','MS-TB','FSSA', ...
         'FSSA with the objective of [2]'};
last_point_arms = {'gwo_sqp_box', 'multistart_tb'};

fid = fopen(outfile, 'w');
if fid < 0, error('export_nfe_log:open', 'cannot write %s', outfile); end
fprintf(fid, ['source,arm,method,nfe_target,seed,n_points,nfe_per_point,', ...
              'nfe_sweep_total,nfe_ratio,count_basis\n']);
nrows = 0;
for f = 1:size(files, 1)
    rel = [files{f,1} '/' which '/' files{f,2}];
    fn  = fullfile(root, files{f,1}, which, files{f,2});
    if ~exist(fn, 'file')
        fprintf('export_nfe_log: %s not found, skipped\n', rel);
        continue;
    end
    D = read_benchmark_csv(fn);
    arms = unique_stable(D.arm);
    for a = 1:numel(arms)
        ia = strcmp(D.arm, arms{a});
        k  = find(strcmp(keys, arms{a}), 1);
        if isempty(k), method = arms{a}; else, method = names{k}; end
        if any(strcmp(last_point_arms, arms{a}))
            basis = 'last_point';
        else
            basis = 'mean_over_points';
        end
        budgets = unique(D.nfe_target(ia));
        for b = budgets(:).'
            ib = ia & D.nfe_target == b;
            seeds = unique(D.seed(ib));
            for s = seeds(:).'
                is = ib & D.seed == s;
                v  = D.nfe_achieved(is);
                if any(abs(v - v(1)) > 1e-9 * max(1, abs(v(1))))
                    error('export_nfe_log:value', ...
                          '%s: nfe_achieved differs within %s, budget %g, seed %g', ...
                          rel, arms{a}, b, s);
                end
                npt = sum(is);
                per = v(1);
                if strcmp(basis, 'mean_over_points')
                    tot = npt * per;
                    if abs(tot - round(tot)) < 1e-3, tot = round(tot); end
                    tot_s = sprintf('%.10g', tot);
                else
                    tot_s = '';
                end
                fprintf(fid, '%s,%s,%s,%.10g,%.10g,%d,%.10g,%s,%.10g,%s\n', ...
                        rel, arms{a}, method, b, s, npt, per, tot_s, per / b, basis);
                nrows = nrows + 1;
            end
        end
    end
end
fclose(fid);
fprintf('export_nfe_log: %d sweeps written to %s\n', nrows, outfile);
end

function u = unique_stable(c)
%UNIQUE_STABLE  Unique strings in order of first appearance (MATLAB and Octave).
u = {};
for i = 1:numel(c)
    if ~any(strcmp(u, c{i})), u{end+1} = c{i}; end %#ok<AGROW>
end
end
