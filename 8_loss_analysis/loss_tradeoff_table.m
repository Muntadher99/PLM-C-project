function T = loss_tradeoff_table(designs, pwm_mf, cfg, dev, outcsv)
%LOSS_TRADEOFF_TABLE  Distortion, compliance, switching count and estimated
%   loss on the operating-point grid (paper, Table 9; Reviewer 1 comments 5
%   and 10, Reviewer 2 comments 4 and 5).
%
%   T = loss_tradeoff_table(designs, [31 63 101])
%   T = loss_tradeoff_table(designs, pwm_mf, cfg, dev, 'loss_table.csv')
%
%   designs(k).label   e.g. 'PLM-C THD-min' | 'PLM-C compliant'
%   designs(k).Ma      n x 1
%   designs(k).TH_deg  n x N angle sets in degrees (NaN rows are skipped)
%   pwm_mf             carrier ratios of the PD carrier PWM references with
%                      min-max injection, evaluated on the same Ma grid
%                      (linear to Ma = 1.155, over-modulating above)
%   cfg, dev           converter and device data, see npc_loss_model
%
%   Writes one row per (design, Ma): design, Ma, fsw_dev_Hz, THD49_pct,
%   THD50_all_pct, maxq, compliant, P_cond_W, P_sw_W, P_loss_W, P_out_W,
%   eta_pct, ITHD_pct, b1, and prints the per-design medians and means.
%   The deployed design of Table 9 is the THD-min design where it complies
%   and the compliant design elsewhere (Section 3.4).

if nargin < 2 || isempty(pwm_mf), pwm_mf = [31 63 101]; end
if nargin < 3, cfg = struct(); end
if nargin < 4, dev = struct(); end
if nargin < 5 || isempty(outcsv), outcsv = 'loss_tradeoff_table.csv'; end

rows = {};
for k = 1:numel(designs)
    d = designs(k);
    for i = 1:numel(d.Ma)
        if any(~isfinite(d.TH_deg(i,:))), continue; end
        L = npc_loss_model(struct('theta_deg', d.TH_deg(i,:)), d.Ma(i), cfg, dev);
        rows(end+1,:) = pack(d.label, d.Ma(i), L); %#ok<AGROW>
    end
end
Ma_grid = unique(vertcat(designs.Ma));
for mf = pwm_mf(:).'
    lab = sprintf('PD-PWM mf=%d', mf);
    for i = 1:numel(Ma_grid)
        L = npc_loss_model(struct('pwm','pd','mf',mf,'thi',true), Ma_grid(i), cfg, dev);
        rows(end+1,:) = pack(lab, Ma_grid(i), L); %#ok<AGROW>
    end
end

hdr = {'design','Ma','fsw_dev_Hz','THD49_pct','THD50_all_pct','maxq','compliant', ...
       'P_cond_W','P_sw_W','P_loss_W','P_out_W','eta_pct','ITHD_pct','b1'};
fid = fopen(outcsv, 'w');
fprintf(fid, '%s', hdr{1}); fprintf(fid, ',%s', hdr{2:end}); fprintf(fid, '\n');
for r = 1:size(rows,1)
    fprintf(fid, '%s', rows{r,1});
    fprintf(fid, ',%.6g', rows{r,2:end}); fprintf(fid, '\n');
end
fclose(fid);

labs = unique(rows(:,1), 'stable');
fprintf('\n%-22s %8s %9s %8s %9s %9s %9s %8s\n', 'design','fsw Hz','medTHD49','compl', ...
        'Pcond W','Psw W','Ploss W','eta %');
for j = 1:numel(labs)
    m = strcmp(rows(:,1), labs{j});
    X = cell2mat(rows(m, 2:end));
    fprintf('%-22s %8.0f %9.3f %4d/%-3d %9.2f %9.2f %9.2f %8.3f\n', labs{j}, ...
        median(X(:,2)), median(X(:,3)), sum(X(:,6)), size(X,1), mean(X(:,7)), ...
        mean(X(:,8)), mean(X(:,9)), mean(X(:,11)));
end
fprintf('wrote %s\n', outcsv);
T = struct('hdr', {hdr}, 'rows', {rows});
end

function r = pack(lab, Ma, L)
r = {lab, Ma, L.fsw_dev, L.THD49, L.THD50_all, L.maxq, double(L.maxq <= 1 + 1e-4), ...
     L.P_cond_T + L.P_cond_D, L.P_sw_T + L.P_rec_D, L.P_loss, L.P_out, L.eta, L.ITHD, L.b1};
end
