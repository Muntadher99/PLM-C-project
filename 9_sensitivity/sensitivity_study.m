function T = sensitivity_study(CFG)
%SENSITIVITY_STUDY  Sensitivity of the reported angle sets to implementation
%   effects (paper, Section 6.1).
%
%   T = sensitivity_study(struct('csv','bench_main.csv','seed',1:10))
%   T = sensitivity_study(struct('csv',f,'studies',{{'deadtime','quant'}}))
%
%   Studies (the first four perturb stored angle sets and re-evaluate them
%   without re-solving; the fifth re-optimises):
%     deadtime  every switching instant shifted by the commutation delay,
%               theta_i -> theta_i + s_i * Delta * sign(sin(theta_i - phi)),
%               Delta = 2*pi*f*t_dead, for dead times t_dead and load angles phi
%     npripple  neutral-point current ripple at the odd triplen orders
%               3, 9, 15, 21 for load angles phi (the dc component is
%               identically zero under half-wave symmetry)
%     quant     angles rounded to a timer of b bits over the quarter wave
%     matol     nearest-entry and linear interpolation between neighbouring
%               operating points of the stored grid
%     mingap    minimum gap delta swept and the design re-optimised at each
%               value (plm_core solve, cold starts plus the solution carried
%               from the previous, looser gap)
%   Perturbed angles are evaluated as they are (no restoration onto b1 = Ma),
%   so both THD49 and the fundamental error are reported.
%
%   CFG fields (defaults in brackets)
%     csv        benchmark CSV holding the angle sets      (required)
%     arm        arm whose angles are used                 ['gn_hybrid']
%     budget     nfe_target of those angles                [20000]
%     seed       seed(s) of those angles; a vector runs every study per seed [1]
%     studies    subset of the five studies                [all]
%     f_fund     fundamental frequency (Hz)                [50]
%     t_dead     dead times (s)                            [0 1 2 3 5 8]*1e-6
%     phi_deg    load angles (deg)                         [0 15 30 45 60]
%     res_bits   timer resolutions (bits per quarter wave) [8 9 10 11 12 14]
%     gap_deg    minimum gaps for the re-optimisation      [0.10 0.25 0.50 1.00 2.00]
%     gap_seeds  independent seeds per (Ma, gap) cell      [3]
%     gap_NFE    evaluations per (Ma, gap) cell            [20000]
%     gap_Ma     operating points of the gap sweep         [1.10 0.90 0.70 0.50 0.30]
%     N, nmax    angles and highest counted order          [15, 49]
%     out        output base name ('' = sensitivity_<date>) ['']
%
%   Output <out>.xlsx with one sheet per study (CSV per sheet always written);
%   summarise with sensitivity_summary(out).

D.csv      = '';
D.arm      = 'gn_hybrid';
D.budget   = 20000;
D.seed     = 1;
D.studies  = {'deadtime','npripple','quant','matol','mingap'};
D.f_fund   = 50;
D.t_dead   = [0 1 2 3 5 8] * 1e-6;
D.phi_deg  = [0 15 30 45 60];
D.res_bits = [8 9 10 11 12 14];
D.gap_deg  = [0.10 0.25 0.50 1.00 2.00];
D.gap_seeds= 3;
D.gap_NFE  = 20000;
D.gap_Ma   = [1.10 0.90 0.70 0.50 0.30];
D.N = 15;  D.nmax = 49;
D.out      = '';
D.verbose  = true;
if nargin < 1, CFG = struct(); end
fn = fieldnames(CFG); C = D;
for i = 1:numel(fn), C.(fn{i}) = CFG.(fn{i}); end
if ischar(C.studies), C.studies = {C.studies}; end
if isempty(C.csv), error('sensitivity_study:csv','set CFG.csv to a run_benchmark output file'); end
if isempty(C.out)
    c = clock(); %#ok<CLOCK>
    C.out = sprintf('sensitivity_%04d%02d%02d_%02d%02d', c(1),c(2),c(3),c(4),c(5));
end

F = plm_core();  P = F.defaults(C.N, C.nmax);
seeds = C.seed(:).';
if C.verbose
    fprintf('\n=================== sensitivity_study ===================\n');
    fprintf('angles  : %s  arm=%s  budget=%g  seeds=%s\n', ...
            C.csv, C.arm, C.budget, mat2str(seeds));
    fprintf('studies : %s\n\n', strjoin(C.studies, ', '));
end

acc = struct();                      % study name -> accumulated rows
hdrs = struct();
for si = 1:numel(seeds)
    sd = seeds(si);
    [MAS, TH] = load_angles(C.csv, C.arm, C.budget, sd, C.N);
    if isempty(MAS)
        warning('sensitivity_study:noangles', ...
            'no accepted angles for arm=%s budget=%g seed=%g -- skipped', C.arm, C.budget, sd);
        continue;
    end
    for k = 1:numel(C.studies)
        nm = lower(C.studies{k});
        switch nm
          case 'deadtime', [h, r] = deadtime_study(F, P, MAS, TH, C);
          case 'npripple', [h, r] = npripple_study(F, P, MAS, TH, C);
          case 'quant',    [h, r] = quant_study(F, P, MAS, TH, C);
          case 'matol',    [h, r] = matol_study(F, P, MAS, TH, C);
          case 'mingap'
            % re-solves from scratch, independent of the stored angles: run once
            if si > 1, continue; end
            [h, r] = mingap_study(C);
          otherwise, warning('sensitivity_study:study','unknown study "%s" skipped', nm); continue;
        end
        f = vname(nm);
        r = [repmat({sd}, size(r,1), 1), r];      % tag every row with its seed
        if ~isfield(acc, f), acc.(f) = {}; hdrs.(f) = [{'angle_seed'}, h(:).']; end
        acc.(f) = [acc.(f); r];
        if C.verbose && si == 1
            fprintf('  %-10s %d rows per seed\n', nm, size(r,1));
        end
    end
end

sheets = {};
fs = fieldnames(acc);
for k = 1:numel(fs)
    sheets(end+1,:) = {fs{k}, hdrs.(fs{k}), acc.(fs{k})}; %#ok<AGROW>
end
if C.verbose
    fprintf('\n  total rows: %d across %d seed(s)\n', ...
            sum(cellfun(@(x) size(x,1), sheets(:,3))), numel(seeds));
end

write_book(C.out, sheets);
T = C.out;
end

% =========================================================================
% DEAD TIME
% =========================================================================
function [hdr, out] = deadtime_study(F, P, MAS, TH, C)
%DEADTIME_STUDY  First-order commutation-delay model: rising and falling
%   edges shifted by Delta according to the load-current polarity at a
%   single sinusoidal load current of angle phi; ideal switches otherwise.
hdr = {'Ma','t_dead_us','phi_deg','delta_deg','THD49_pct','dTHD49_pp', ...
       'b1','b1_err','b1_err_pct_of_Ma','MAXQ','compliant','minGap_violated'};
out = {};
for m = 1:numel(MAS)
    Ma = MAS(m);  th0 = TH(m,:);
    R0 = F.report(th0, Ma, P);
    for td = C.t_dead(:).'
        Delta = 2*pi*C.f_fund*td;
        for ph = C.phi_deg(:).'
            phi = deg2rad(ph);
            sgn = sign(sin(th0 - phi));  sgn(sgn == 0) = 1;
            th  = th0 + Delta .* P.s .* sgn;
            R   = F.report(th, Ma, P);
            gapv = any(diff(sort(th)) < P.minGap - 1e-12) || th(1) < P.lb || th(end) > P.ub;
            out(end+1,:) = {Ma, td*1e6, ph, rad2deg(Delta), 100*R.THD49, ...
                            100*(R.THD49 - R0.THD49), R.b1, R.b1err, ...
                            100*R.b1err/max(abs(Ma),eps), R.MAXQ, ...
                            double(R.compliant), double(gapv)}; %#ok<AGROW>
        end
    end
end
end

% =========================================================================
% NEUTRAL-POINT RIPPLE
% =========================================================================
function [hdr, out] = npripple_study(F, P, MAS, TH, C) %#ok<INUSL>
%NPRIPPLE_STUDY  Neutral-point current at the odd triplen orders h:
%   |i_NP,h| = (3/2) sqrt((D_{h-1}-D_{h+1})^2 cos^2(phi) + (D_{h-1}+D_{h+1})^2 sin^2(phi)),
%   D_k = 4/(k*pi) * sum_i s_i sin(k*theta_i) (k even), per unit load current.
H = [3 9 15 21];
hdr = {'Ma','phi_deg','NPR_rms','NPR_h3','NPR_h9','NPR_h15','NPR_h21','THD49_pct'};
out = {};
for m = 1:numel(MAS)
    Ma = MAS(m);  th = TH(m,:);
    R  = F.report(th, Ma, P);
    for ph = C.phi_deg(:).'
        phi = deg2rad(ph);
        amp = zeros(1, numel(H));
        for j = 1:numel(H)
            h  = H(j);
            Dm = Dk(th, P, h-1);
            Dp = Dk(th, P, h+1);
            amp(j) = 1.5 * sqrt( (Dm-Dp)^2 * cos(phi)^2 + (Dm+Dp)^2 * sin(phi)^2 );
        end
        out(end+1,:) = {Ma, ph, sqrt(mean(amp.^2)), amp(1), amp(2), amp(3), ...
                        amp(4), 100*R.THD49}; %#ok<AGROW>
    end
end
end

function v = Dk(th, P, k)
%DK  Even-harmonic coefficient of the zero-state indicator (k even).
v = (4/(k*pi)) * sum(P.s(:).' .* sin(k*th(:).'));
end

% =========================================================================
% TIMER QUANTISATION
% =========================================================================
function [hdr, out] = quant_study(F, P, MAS, TH, C)
%QUANT_STUDY  Angles rounded to an LSB of 90/2^b degrees and evaluated as
%   rounded. gap_deficit_deg and deficit_over_LSB show whether a minimum-gap
%   violation is a sub-LSB rounding of an active constraint.
hdr = {'Ma','bits','LSB_deg','THD49_pct','dTHD49_pp','b1_err','MAXQ', ...
       'compliant','minGap_violated','min_gap_deg','gap_deficit_deg','deficit_over_LSB'};
out = {};
for m = 1:numel(MAS)
    Ma = MAS(m);  th0 = TH(m,:);
    R0 = F.report(th0, Ma, P);
    for b = C.res_bits(:).'
        lsb = (pi/2) / 2^b;
        th  = round(th0 / lsb) * lsb;
        R   = F.report(th, Ma, P);
        g   = min(diff(sort(th)));
        def = max(0, P.minGap - g);
        out(end+1,:) = {Ma, b, rad2deg(lsb), 100*R.THD49, ...
                        100*(R.THD49 - R0.THD49), R.b1err, R.MAXQ, ...
                        double(R.compliant), double(g < P.minGap - 1e-12), ...
                        rad2deg(g), rad2deg(def), def/lsb}; %#ok<AGROW>
    end
end
end

% =========================================================================
% MODULATION-INDEX GRID
% =========================================================================
function [hdr, out] = matol_study(F, P, MAS, TH, C) %#ok<INUSL>
%MATOL_STUDY  Between two stored operating points: the nearest stored entry
%   and the linear interpolation of the angles (sorted, not re-solved), at
%   25%, 50% and 75% of the interval.
hdr = {'Ma_target','Ma_lo','Ma_hi','method','THD49_pct','b1','b1_err', ...
       'b1_err_pct','MAXQ','compliant'};
out = {};
[ms, ix] = sort(MAS(:));  TT = TH(ix,:);
for i = 1:numel(ms)-1
    lo = ms(i);  hi = ms(i+1);
    for f = [0.25 0.5 0.75]
        tgt = lo + f*(hi-lo);
        if f < 0.5, thn = TT(i,:); else, thn = TT(i+1,:); end
        Rn = F.report(thn, tgt, P);
        out(end+1,:) = {tgt, lo, hi, 'nearest', 100*Rn.THD49, Rn.b1, Rn.b1err, ...
                        100*Rn.b1err/max(tgt,eps), Rn.MAXQ, double(Rn.compliant)}; %#ok<AGROW>
        thi = sort((1-f)*TT(i,:) + f*TT(i+1,:));
        Ri  = F.report(thi, tgt, P);
        out(end+1,:) = {tgt, lo, hi, 'interp', 100*Ri.THD49, Ri.b1, Ri.b1err, ...
                        100*Ri.b1err/max(tgt,eps), Ri.MAXQ, double(Ri.compliant)}; %#ok<AGROW>
    end
end
end

% =========================================================================
% MINIMUM GAP (re-optimised)
% =========================================================================
function [hdr, out] = mingap_study(C)
%MINGAP_STUDY  Re-optimise at each minimum gap, loosest first. Every cell
%   is solved from gap_seeds cold multi-starts plus the best solution of the
%   previous (looser) gap projected onto the tighter set; a backward pass
%   then applies the nesting of the feasible sets (a solution feasible at a
%   tighter gap is feasible, with the same objective, at every looser gap).
%   source: 0 cold start, 1 carried, 2 inherited from a tighter gap.
hdr = {'Ma','minGap_deg','THD49_pct','dTHD49_pp','MAXQ','compliant','accepted','NFE','source'};
out = {};
Fg  = plm_core();
gaps = sort(C.gap_deg(:).');          % loosest first
nMa  = numel(C.gap_Ma);

for mi = 1:nMa
    Ma   = C.gap_Ma(mi);
    prev = [];                         % best angles at the previous, looser gap
    nG = numel(gaps);
    cell_thd = inf(1,nG);  cell_R = cell(1,nG);  cell_c = zeros(1,nG);
    cell_g   = gaps;       cell_nfe = zeros(1,nG);
    for gi = 1:numel(gaps)
        g  = gaps(gi);
        Pg = Fg.defaults(C.N, C.nmax, g);
        if Pg.span <= 0
            warning('sensitivity_study:mingap','minGap %g deg leaves no feasible set, skipped', g);
            continue;
        end

        best = Inf;  bR = []; bX = []; viaCarry = 0;

        % 1. the carried solution, projected onto the tighter set
        if ~isempty(prev)
            [THc, Ic] = Fg.solve(Fg.project(prev, Pg), Ma, Pg, struct('mode','thd'));
            if Ic.ok(1) && Ic.thd(1) < best
                best = Ic.thd(1);  bX = THc(1,:);
                bR = Fg.report(bX, Ma, Pg);  viaCarry = 1;
            end
        end

        % 2. independent cold seeds sharing the budget
        per = max(1, floor(C.gap_NFE / max(C.gap_seeds,1)));
        for sd = 1:max(C.gap_seeds,1)
            rng_set(sd*7919 + mi*131 + gi);
            n0 = Fg.nfe();
            while Fg.nfe() - n0 < per
                [TH, I] = Fg.solve(Fg.zsample(rand(24, C.N), Pg), Ma, Pg, struct('mode','thd'));
                ok = find(I.ok);
                for j = ok(:).'
                    if I.thd(j) < best
                        best = I.thd(j);  bX = TH(j,:);
                        bR = Fg.report(bX, Ma, Pg);  viaCarry = 0;
                    end
                end
            end
        end

        if isempty(bR)
            cell_thd(gi) = NaN;  cell_R{gi} = [];  cell_nfe(gi) = Fg.nfe();
            continue;
        end
        prev = bX;
        cell_thd(gi) = best;  cell_R{gi} = bR;  cell_c(gi) = viaCarry; %#ok<AGROW>
        cell_g(gi) = g;  cell_nfe(gi) = Fg.nfe(); %#ok<AGROW>
    end

    % backward pass over the nested feasible sets
    for gi = numel(cell_thd)-1 : -1 : 1
        if isfinite(cell_thd(gi+1)) && cell_thd(gi+1) < cell_thd(gi)
            cell_thd(gi) = cell_thd(gi+1);
            cell_R{gi}   = cell_R{gi+1};
            cell_c(gi)   = 2;
        end
    end

    b0 = cell_thd(1);
    for gi = 1:numel(cell_thd)
        if ~isfinite(cell_thd(gi)) || isempty(cell_R{gi})
            out(end+1,:) = {Ma, cell_g(gi), NaN, NaN, NaN, 0, 0, cell_nfe(gi), 0}; %#ok<AGROW>
            continue;
        end
        R = cell_R{gi};
        out(end+1,:) = {Ma, cell_g(gi), 100*cell_thd(gi), 100*(cell_thd(gi)-b0), ...
                        R.MAXQ, double(R.compliant), 1, cell_nfe(gi), cell_c(gi)}; %#ok<AGROW>
    end
end
end

function rng_set(s)
try, rng(s,'twister'); catch, rand('state',s); randn('state',s); end %#ok<RAND>
end

% =========================================================================
% INPUT / OUTPUT
% =========================================================================
function [MAS, TH] = load_angles(fn, arm, budget, seed, N)
%LOAD_ANGLES  Accepted angle sets of one arm, budget and seed (radians, sorted).
fid = fopen(fn,'r');
if fid < 0, error('sensitivity_study:open','cannot open %s', fn); end
h = strsplit(strtrim(fgetl(fid)), ',');
ci = @(nm) find(strcmp(h, nm), 1);
c_arm = ci('arm'); c_b = ci('nfe_target'); c_s = ci('seed');
c_ma  = ci('Ma');  c_ok = ci('accepted');
c_th  = zeros(1,N);
for k = 1:N, c_th(k) = ci(sprintf('theta%02d_deg', k)); end
MAS = [];  TH = [];
while true
    l = fgetl(fid);
    if ~ischar(l), break; end
    p = strsplit(strtrim(l), ',', 'CollapseDelimiters', false);
    if ~strcmpi(p{c_arm}, arm), continue; end
    if str2double(p{c_b}) ~= budget || str2double(p{c_s}) ~= seed, continue; end
    if str2double(p{c_ok}) ~= 1, continue; end
    v = zeros(1,N);  bad = false;
    for k = 1:N
        v(k) = str2double(p{c_th(k)});
        if ~isfinite(v(k)), bad = true; break; end
    end
    if bad, continue; end
    MAS(end+1,1) = str2double(p{c_ma}); %#ok<AGROW>
    TH(end+1,:)  = deg2rad(sort(v));    %#ok<AGROW>
end
fclose(fid);
end

function write_book(base, sheets)
%WRITE_BOOK  One CSV per sheet always; an .xlsx when the host can write one.
for i = 1:size(sheets,1)
    wcell(sprintf('%s_%s.csv', base, sheets{i,1}), [sheets{i,2}(:).'; sheets{i,3}]);
end
fx = [base '.xlsx'];  ok = false;
if exist('writecell','file') == 2 || exist('writetable','file') == 2
    try
        if exist(fx,'file') == 2, delete(fx); end
        for i = 1:size(sheets,1)
            M = blankify([sheets{i,2}(:).'; sheets{i,3}]);
            if exist('writecell','file') == 2
                writecell(M, fx, 'Sheet', sheets{i,1});
            else
                vn = cellfun(@vname, M(1,:), 'UniformOutput', false);
                writetable(cell2table(M(2:end,:), 'VariableNames', vn), fx, 'Sheet', sheets{i,1});
            end
        end
        ok = true;
    catch e
        fprintf('  (xlsx skipped: %s -- the CSVs were still written)\n', e.message);
    end
end
if ok, fprintf('\nwrote %s   (one sheet per study)\n\n', fx);
else,  fprintf('\nwrote %s_*.csv\n\n', base); end
end

function M = blankify(M)
for i = 1:numel(M)
    v = M{i};
    if isnumeric(v) && (isempty(v) || ~all(isfinite(v(:)))), M{i} = ''; end
end
end

function wcell(fn, M)
fid = fopen(fn,'w');
for r = 1:size(M,1)
    for c = 1:size(M,2)
        if c > 1, fprintf(fid, ','); end
        v = M{r,c};
        if ischar(v)
            if any(v == ','), fprintf(fid, '"%s"', v); else, fprintf(fid, '%s', v); end
        elseif isempty(v) || (isnumeric(v) && ~isfinite(v))
        else
            d = double(v);
            if d == fix(d) && abs(d) < 1e15, fprintf(fid, '%d', d);
            else,                            fprintf(fid, '%.6g', d); end
        end
    end
    fprintf(fid, '\n');
end
fclose(fid);
end

function s = vname(s)
s = regexprep(strtrim(s), '\W', '_');
if isempty(s) || ~isletter(s(1)), s = ['v' s]; end
end
