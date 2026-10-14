function sensitivity_summary(base)
%SENSITIVITY_SUMMARY  Print the tables of a sensitivity_study run.
%
%   sensitivity_summary('results/sensitivity')   % base name, no extension
%
%   Reads <base>_deadtime.csv, _npripple.csv, _mingap.csv, _quant.csv and
%   _matol.csv. Distortion columns are summarised over operating points as
%   median and worst case; minimum-gap violation counts are reported as a
%   range over the seeds (whether a solution sits exactly on the gap bound
%   is a property of the individual angle set).

if nargin < 1, error('sensitivity_summary:arg','usage: sensitivity_summary(''sensitivity'')'); end
W = 92;

% ---------------------------------------------------------------- dead time
[h, r] = rd([base '_deadtime.csv']);
if ~isempty(r)
    banner(W, 'DEAD TIME');
    td = uq(col(h,r,'t_dead_us'));
    fprintf('%9s %9s %11s %11s %12s %12s %14s\n', 't_dead us','delta deg', ...
            'med dTHD pp','max dTHD pp','med b1err %','max b1err %','minGap viol');
    for i = 1:numel(td)
        m = col(h,r,'t_dead_us') == td(i);
        d  = col(h,r,'dTHD49_pp');           d  = d(m);
        b  = abs(col(h,r,'b1_err_pct_of_Ma')); b = b(m);
        dg = col(h,r,'delta_deg');           dg = dg(m);
        [lo, hi, ns] = viol_range(h, r, m);
        fprintf('%9.0f %9.4f %11.3f %11.3f %12.3f %12.3f %8d - %-3d (%d seeds)\n', ...
                td(i), dg(1), md(d), mx(d), md(b), mx(b), lo, hi, ns);
    end
    phi_table(h, r, W);
end

% ---------------------------------------------------------------- NP ripple
[h, r] = rd([base '_npripple.csv']);
if ~isempty(r)
    banner(W, 'NEUTRAL-POINT RIPPLE');
    v = col(h,r,'NPR_rms');
    fprintf('  rms ripple over all operating points and load angles:\n');
    fprintf('    min %.5f   median %.5f   max %.5f   ratio %.1fx\n', ...
            mn(v), md(v), mx(v), mx(v)/max(mn(v),1e-12));
    fprintf('  (the dc component of the neutral-point current is identically zero\n');
    fprintf('   under half-wave symmetry; only the odd triplen ripple remains)\n');
end

% ---------------------------------------------------------------- minimum gap
[h, r] = rd([base '_mingap.csv']);
if ~isempty(r)
    banner(W, 'MINIMUM GAP (re-optimised)');
    g = uq(col(h,r,'minGap_deg'));
    fprintf('%11s %13s %13s %13s\n','minGap deg','med dTHD pp','max dTHD pp','worst Ma');
    for i = 1:numel(g)
        m = col(h,r,'minGap_deg') == g(i);
        d = col(h,r,'dTHD49_pp');  ma = col(h,r,'Ma');
        dm = d(m);  mam = ma(m);
        [~, j] = max(dm);
        fprintf('%11.2f %13.4f %13.4f %13.2f\n', g(i), md(dm), mx(dm), mam(j));
    end
    bad = sum(col(h,r,'dTHD49_pp') < -1e-6);
    if bad > 0
        fprintf('\n  %d rows show a tighter gap below the loosest one (under-converged cells).\n', bad);
    else
        fprintf('\n  Monotone in the minimum gap, as the nesting of the feasible sets requires.\n');
    end
end

% ---------------------------------------------------------------- quantisation
[h, r] = rd([base '_quant.csv']);
if ~isempty(r)
    banner(W, 'TIMER QUANTISATION');
    b = uq(col(h,r,'bits'));
    fprintf('%7s %10s %13s %14s %16s\n','bits','LSB deg','max dTHD pp','minGap viol','max deficit/LSB');
    for i = 1:numel(b)
        m = col(h,r,'bits') == b(i);
        d = col(h,r,'dTHD49_pp');  L = col(h,r,'LSB_deg');
        q = col(h,r,'deficit_over_LSB');
        [lo, hi, ~] = viol_range(h, r, m);
        Lm = L(m);
        fprintf('%7d %10.4f %13.4f %8d - %-3d %16.3f\n', ...
                b(i), Lm(1), mx(d(m)), lo, hi, mx(q(m)));
    end
end

% ---------------------------------------------------------------- Ma grid
[h, r] = rd([base '_matol.csv']);
if ~isempty(r)
    banner(W, 'MODULATION-INDEX GRID (between stored points)');
    meth = unique(txt(h,r,'method'));
    fprintf('%10s %12s %12s %13s %13s %12s\n', ...
            'method','med THD %','max THD %','med b1err %','max b1err %','compliant');
    for i = 1:numel(meth)
        m = strcmp(txt(h,r,'method'), meth{i});
        t = col(h,r,'THD49_pct');  be = abs(col(h,r,'b1_err_pct'));
        c = col(h,r,'compliant');
        fprintf('%10s %12.3f %12.3f %13.3f %13.3f %7d/%-4d\n', meth{i}, ...
                md(t(m)), mx(t(m)), md(be(m)), mx(be(m)), sum(c(m)==1), sum(m));
    end
end
fprintf('\n%s\n\n', repmat('=',1,W));
end

% =========================================================================
function phi_table(h, r, W) %#ok<INUSD>
p = uq(col(h,r,'phi_deg'));  td = uq(col(h,r,'t_dead_us'));
t3 = td(min(numel(td), max(1, find(td >= 3, 1))));   % a representative dead time
fprintf('\n  at t_dead = %g us:\n', t3);
fprintf('%9s', 'phi deg'); fprintf('%12.0f', p); fprintf('\n');
fprintf('%9s', 'max dTHD');
for i = 1:numel(p)
    m = col(h,r,'t_dead_us') == t3 & col(h,r,'phi_deg') == p(i);
    d = col(h,r,'dTHD49_pp');
    fprintf('%12.3f', mx(d(m)));
end
fprintf('\n%9s', 'max b1err');
for i = 1:numel(p)
    m = col(h,r,'t_dead_us') == t3 & col(h,r,'phi_deg') == p(i);
    b = abs(col(h,r,'b1_err_pct_of_Ma'));
    fprintf('%12.3f', mx(b(m)));
end
fprintf('\n');
end

function [lo, hi, ns] = viol_range(h, r, m)
%VIOL_RANGE  Minimum-gap violations counted per seed, reported as a range.
v = col(h,r,'minGap_violated');
if isfield_(h,'angle_seed')
    s = col(h,r,'angle_seed');  us = uq(s);  c = zeros(1,numel(us));
    for i = 1:numel(us), c(i) = sum(v(m & s == us(i))); end
    lo = min(c); hi = max(c); ns = numel(us);
else
    lo = sum(v(m)); hi = lo; ns = 1;
end
end

function banner(W, s)
fprintf('\n%s\n%s\n%s\n', repmat('=',1,W), s, repmat('=',1,W));
end

function [h, r] = rd(fn)
h = struct(); r = {};
fid = fopen(fn,'r');
if fid < 0, return; end
hdr = strsplit(strtrim(fgetl(fid)), ',');
for i = 1:numel(hdr), h.(vn(hdr{i})) = i; end
while true
    l = fgetl(fid);
    if ~ischar(l), break; end
    r{end+1} = strsplit(strtrim(l), ',', 'CollapseDelimiters', false); %#ok<AGROW>
end
fclose(fid);
end

function v = col(h, r, nm)
f = vn(nm);
if ~isfield(h, f), v = nan(numel(r),1); return; end
k = h.(f);  v = nan(numel(r),1);
for i = 1:numel(r)
    if k <= numel(r{i}), v(i) = str2double(r{i}{k}); end
end
end

function v = txt(h, r, nm)
k = h.(vn(nm));  v = cell(numel(r),1);
for i = 1:numel(r), v{i} = r{i}{k}; end
end

function t = isfield_(h, nm), t = isfield(h, vn(nm)); end
function s = vn(s), s = regexprep(strtrim(s),'\W','_'); if ~isletter(s(1)), s=['v' s]; end, end
function u = uq(x), u = unique(x(isfinite(x))).'; end
function y = md(x), x = sort(x(isfinite(x))); if isempty(x), y=NaN; else, y=x(ceil(end/2)); end, end
function y = mx(x), x = x(isfinite(x)); if isempty(x), y=NaN; else, y=max(x); end, end
function y = mn(x), x = x(isfinite(x)); if isempty(x), y=NaN; else, y=min(x); end, end
