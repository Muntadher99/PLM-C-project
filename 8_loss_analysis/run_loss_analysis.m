function run_loss_analysis(witness_dir, outdir)
%RUN_LOSS_ANALYSIS  Reproduce the loss estimate of Section 6.3 and Table 9
%   (Reviewer 1 comment 5, Reviewer 2 comments 4 and 5).
%
%   run_loss_analysis()                      % designs from 6_compliance_witnesses/results
%   run_loss_analysis(witness_dir, outdir)
%
%   The two programmed designs on the 21-point grid are read from the audited
%   witness tables written by build_witness_tables:
%     'PLM-C THD-min'    thdmin_selected.csv
%     'PLM-C compliant'  rows 'compliant_design' of witnesses_selected.csv
%   and compared with phase-disposition carrier PWM with min-max injection at
%   carrier ratios 31, 63 and 101 (loss_tradeoff_table, npc_loss_model).
%
%   Outputs in outdir [8_loss_analysis/reproduced]
%     loss_table.csv        device data at Tvj = 25 C (Table 9)
%     loss_table_175C.csv   device data at Tvj = 175 C (bracketing case quoted
%                           in Section 6.3)
%   paper_tables prints Table 9, including the deployed design, from
%   loss_table.csv.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
if nargin < 1 || isempty(witness_dir)
    witness_dir = fullfile(root, '6_compliance_witnesses', 'results');
end
if nargin < 2 || isempty(outdir), outdir = fullfile(here, 'reproduced'); end
if ~exist(outdir, 'dir'), mkdir(outdir); end

[grid, TH]  = read_designs(fullfile(witness_dir, 'thdmin_selected.csv'), 'thd_min', []);
[~, THc]    = read_designs(fullfile(witness_dir, 'witnesses_selected.csv'), 'compliant_design', grid);
des(1).label = 'PLM-C THD-min';   des(1).Ma = grid;  des(1).TH_deg = TH;
des(2).label = 'PLM-C compliant'; des(2).Ma = grid;  des(2).TH_deg = THc;

fprintf('\n== Table 9: device data at Tvj = 25 C ==\n');
loss_tradeoff_table(des, [31 63 101], struct(), struct(), fullfile(outdir, 'loss_table.csv'));
fprintf('\n== bracketing case: device data at Tvj = 175 C ==\n');
loss_tradeoff_table(des, [31 63 101], struct(), struct('Tvj', 175), ...
                    fullfile(outdir, 'loss_table_175C.csv'));
end

% =========================================================================
function [Ma, TH] = read_designs(csvfile, label, grid)
%READ_DESIGNS  Angle sets of one label from an audit_angle_sets CSV. With a
%   grid, rows are placed on that grid (NaN where the label has no row).
txt = fileread(csvfile);
L = regexp(txt, '\r?\n', 'split');  L = L(~cellfun(@isempty, strtrim(L)));
H = regexp(strtrim(L{1}), ',', 'split');
tc = find(strncmp(H, 'theta', 5));  mc = find(strcmp(H, 'Ma'));
Ma = [];  TH = zeros(0, numel(tc));
for r = 2:numel(L)
    v = regexp(L{r}, ',', 'split');
    if ~strcmp(v{1}, label), continue; end
    Ma(end+1,1) = str2double(v{mc});             %#ok<AGROW>
    TH(end+1,:) = str2double(v(tc));             %#ok<AGROW>
end
if ~isempty(grid)
    G = nan(numel(grid), numel(tc));
    for i = 1:numel(grid)
        k = find(abs(Ma - grid(i)) < 1e-6, 1);
        if ~isempty(k), G(i,:) = TH(k,:); end
    end
    Ma = grid;  TH = G;
end
end
