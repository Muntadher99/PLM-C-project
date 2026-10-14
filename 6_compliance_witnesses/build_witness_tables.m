function W = build_witness_tables(minimax_csv, bench_csv, budget, outdir)
%BUILD_WITNESS_TABLES  Audited compliance witnesses, Table 7 and the angle
%   sets of the device-level and loss studies (paper, Sections 5.7, 6.2 and
%   6.3; Reviewer 1 comments 2, 3 and 10, Reviewer 2 comments 2 and 5).
%
%   W = build_witness_tables('minimax_search.csv', 'bench_main.csv')
%   W = build_witness_tables(minimax_csv, bench_csv, 20000, 'reproduced')
%
%   Inputs
%     minimax_csv  output of minimax_witness_search (21 Ma x 5 seeds x 3 modes);
%                  the rows of mode 'minimax' are the witnesses
%     bench_csv    main benchmark CSV; the distortion-minimising design at each
%                  Ma is the lowest-THD49 accepted PLM-C (gn_hybrid) solution
%                  over the ten seeds at BUDGET
%     budget       NFE per operating point of that design [20000]
%     outdir       folder of the output files [current folder]
%
%   Selection per operating point
%     certificate       the witness with the lowest loading Q^ = max_n q_n
%     compliant design  the compliant witness (q_max <= 1 + 1e-4) of lowest THD49
%     deployed design   the distortion-minimising design where it complies,
%                       otherwise the compliant design (Section 3.4)
%
%   Outputs
%     witnesses_all_audit.csv   every minimax angle set, audited (audit_angle_sets)
%     witnesses_selected.csv    certificate and compliant design per Ma, audited
%     thdmin_selected.csv       distortion-minimising design per Ma, audited
%     unified_table.csv         Table 7, one row per Ma
%     nonideal_cases.csv        the ten angle sets of the device-level study
%                               (Ma = 0.40 0.60 0.80 1.00 1.20, both designs)
%   The loss study (Table 9) reads thdmin_selected.csv and
%   witnesses_selected.csv, see run_loss_analysis.

if nargin < 3 || isempty(budget), budget = 20000; end
if nargin < 4 || isempty(outdir), outdir = '.'; end
if ~exist(outdir, 'dir'), mkdir(outdir); end
out = @(f) fullfile(outdir, f);

% ---------------- minimax witnesses ----------------------------------------
[H, C] = read_long(minimax_csv);
mcol = find(ismember(lower(H), {'mode','arm','strategy','objective'}), 1);
if isempty(mcol), error('build_witness_tables:mode', 'no mode column in %s', minimax_csv); end
mode = C{mcol};
isMM = ~cellfun(@isempty, strfind(lower(mode), 'minimax')); %#ok<STRCLFH>
Ma   = str2double(C{strcmp(H,'Ma')});
seed = str2double(C{strcmp(H,'seed')});
tc   = find(strncmp(H, 'theta', 5));
TH   = zeros(numel(Ma), numel(tc));
for j = 1:numel(tc), TH(:,j) = str2double(C{tc(j)}); end
ok   = isMM & all(isfinite(TH), 2);
labs = arrayfun(@(s) sprintf('minimax_seed%d', s), seed(ok), 'UniformOutput', false);
Aall = audit_angle_sets(Ma(ok), TH(ok,:), labs, out('witnesses_all_audit.csv'));

grid = sort(unique(round(Ma(ok)*1e6)/1e6), 'descend');
nG = numel(grid);
selC = nan(nG,1);  selD = nan(nG,1);
for i = 1:nG
    k = find(abs(Aall.Ma - grid(i)) < 1e-6 & Aall.feasible);
    if isempty(k), continue; end
    [~, j] = min(Aall.qmax(k));  selC(i) = k(j);                       % certificate
    c = k(Aall.compliant(k));
    if ~isempty(c), [~, j] = min(Aall.THD49(c)); selD(i) = c(j); end   % compliant design
end
fprintf('\ncompliant witnesses (q_max <= 1) per operating point:\n');
for i = 1:nG
    k = abs(Aall.Ma - grid(i)) < 1e-6 & Aall.feasible;
    fprintf('  Ma = %.2f : %d of %d seeds\n', grid(i), sum(Aall.compliant(k)), sum(k));
end

% ---------------- distortion-minimising designs ----------------------------
D = read_benchmark_csv(bench_csv);
THm = nan(nG, D.N);  nComp = zeros(nG,1);
for i = 1:nG
    k = find(strcmp(D.arm,'gn_hybrid') & D.nfe_target == budget & ...
             abs(D.Ma - grid(i)) < 1e-6 & D.accepted == 1);
    if isempty(k), continue; end
    [~, j] = min(D.THD49(k));  THm(i,:) = D.theta(k(j),:);
    nComp(i) = sum(D.compliant(k) == 1);
end
Am = audit_angle_sets(grid, THm, 'thd_min', out('thdmin_selected.csv'));

% ---------------- selected witnesses (certificate and design) --------------
kc = selC(isfinite(selC));  kd = selD(isfinite(selD));
audit_angle_sets([Aall.Ma(kc); Aall.Ma(kd)], [Aall.theta(kc,:); Aall.theta(kd,:)], ...
    [repmat({'certificate'}, numel(kc), 1); repmat({'compliant_design'}, numel(kd), 1)], ...
    out('witnesses_selected.csv'));

% ---------------- Table 7 --------------------------------------------------
fid = fopen(out('unified_table.csv'), 'w');
fprintf(fid, ['Ma,THDmin_THD49_pct,THDmin_qmax,THDmin_compliant,seeds_compliant_of10,' ...
              'Qhat,Qhat_binding,design_THD49_pct,design_qmax,design_binding,' ...
              'dTHD49_pp,design_b1_err,design_min_gap_deg\n']);
fprintf('\n%5s | %7s %6s %3s | %6s %-9s | %7s %6s | %7s\n', 'Ma', 'THD49', 'qmax', 'ok', ...
        'Qhat', 'binding', 'THD49', 'qmax', 'dTHD pp');
for i = 1:nG
    c = selC(i);  d = selD(i);
    if isnan(c) || isnan(d), continue; end
    fprintf(fid, '%.2f,%.4f,%.4f,%d,%d,%.4f,%s,%.4f,%.4f,%s,%.4f,%.2e,%.4f\n', grid(i), ...
        Am.THD49(i), Am.qmax(i), Am.compliant(i), nComp(i), Aall.qmax(c), Aall.binding{c}, ...
        Aall.THD49(d), Aall.qmax(d), Aall.binding{d}, Aall.THD49(d) - Am.THD49(i), ...
        Aall.b1_err(d), Aall.min_gap_deg(d));
    fprintf('%5.2f | %7.3f %6.3f %3d | %6.3f %-9s | %7.3f %6.3f | %+7.3f\n', grid(i), ...
        Am.THD49(i), Am.qmax(i), Am.compliant(i), Aall.qmax(c), Aall.binding{c}, ...
        Aall.THD49(d), Aall.qmax(d), Aall.THD49(d) - Am.THD49(i));
end
fclose(fid);
fprintf('\nTHD-min design compliant at %d of %d points (some seed compliant at %d).\n', ...
        sum(Am.compliant), nG, sum(nComp > 0));
fprintf('witness with Q^ <= 1 at %d of %d points; Q^ range %.3f .. %.3f\n', sum(isfinite(selC) & ...
        Aall.qmax(max(selC,1)) <= 1+1e-4), nG, min(Aall.qmax(kc)), max(Aall.qmax(kc)));
fprintf('median THD49: THD-min %.3f%%, compliant design %.3f%%, all compliant minimax rows %.3f%%\n', ...
        median(Am.THD49), median(Aall.THD49(kd)), median(Aall.THD49(Aall.compliant)));

% ---------------- inputs of the device-level study -------------------------
pick = [0.40 0.60 0.80 1.00 1.20];
fid = fopen(out('nonideal_cases.csv'), 'w');
fprintf(fid, 'design,Ma'); fprintf(fid, ',theta%02d_deg', 1:D.N); fprintf(fid, '\n');
for p = pick
    i = find(abs(grid - p) < 1e-6, 1);
    if isempty(i) || isnan(selD(i)), continue; end
    fprintf(fid, 'thd_min,%.2f', p);  fprintf(fid, ',%.10f', Am.theta(i,:));  fprintf(fid, '\n');
    fprintf(fid, 'compliant,%.2f', p); fprintf(fid, ',%.10f', Aall.theta(selD(i),:)); fprintf(fid, '\n');
end
fclose(fid);
fprintf('wrote the witness tables to %s\n', outdir);

W = struct('grid', grid, 'thdmin', Am, 'witnesses', Aall, 'selC', selC, 'selD', selD);
end

% =========================================================================
function [H, C] = read_long(f)
txt = fileread(f);  L = regexp(txt, '\r?\n', 'split');  L = L(~cellfun(@isempty, strtrim(L)));
H = regexp(strtrim(L{1}), ',', 'split');  nc = numel(H);  nr = numel(L) - 1;
C = repmat({cell(nr,1)}, 1, nc);  bad = false(nr,1);
for r = 1:nr
    v = regexp(L{r+1}, ',', 'split');
    if numel(v) ~= nc, bad(r) = true; continue; end
    for c = 1:nc, C{c}{r} = v{c}; end
end
for c = 1:nc, C{c} = C{c}(~bad); end
end
