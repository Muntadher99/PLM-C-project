function R = evaluate_candidates(TH, Ma, P, F)
%EVALUATE_CANDIDATES  External restoration, acceptance and measurement applied
%   identically to the angles returned by every arm (paper, Section 4).
%
%   R = evaluate_candidates(TH, Ma, P, F)
%
%   TH : K x N switching angles in radians (a NaN row = arm returned nothing)
%   Ma : K x 1 modulation index of each row
%   P  : problem constants, F.defaults(N, nmax)
%   F  : plm_core() handles
%
%   Steps, the same for every arm:
%     1. RESTORE  five minimum-norm Newton steps onto {b1 = Ma}, each
%                 followed by the Euclidean projection onto (4). Not charged
%                 to the arm's NFE (Table 2, "All arms").
%     2. ACCEPT   |b1 - Ma| < P.FEAS_B1 (1e-8) and constraint violation < 1e-9.
%     3. MEASURE  every reported quantity from the restored angles.
%
%   R fields (K x 1 unless noted)
%     theta_deg (K x N), accepted, THD49, THD99, WTHD49, WTHD999, OOB_max,
%     OOB_order, b3_ratio, MAXQ, b1, b1err, viol, nviol, compliant

if nargin < 3 || isempty(P), P = F.defaults(); end
[K, N] = size(TH);
Ma = Ma(:);
bad = any(~isfinite(TH), 2);

R.theta_deg = nan(K, N);
R.accepted  = false(K,1);
f = {'THD49','THD99','WTHD49','WTHD999','OOB_max','OOB_order','b3_ratio', ...
     'MAXQ','b1','b1err','viol','nviol'};
for i = 1:numel(f), R.(f{i}) = nan(K,1); end
R.compliant = false(K,1);

for k = 1:K
    if bad(k), continue; end
    x = restore(TH(k,:), Ma(k), P, F);                              % 1. restore
    Rr = F.report(x, Ma(k), P);
    Sp = F.spectrum(x, P, 999);
    R.accepted(k)  = (Rr.b1err < P.FEAS_B1) && (Rr.viol < 1e-9);   % 2. accept
    R.theta_deg(k,:) = rad2deg(x);                                  % 3. measure
    R.THD49(k)=Rr.THD49; R.MAXQ(k)=Rr.MAXQ;  R.b1(k)=Rr.b1;
    R.b1err(k)=Rr.b1err; R.viol(k)=Rr.viol;  R.nviol(k)=Rr.nviol;
    R.compliant(k)=Rr.compliant;
    R.THD99(k)=Sp.THD99;  R.WTHD49(k)=Sp.WTHD49;  R.WTHD999(k)=Sp.WTHDfull;
    R.OOB_max(k)=Sp.OOB_max; R.OOB_order(k)=Sp.OOB_order; R.b3_ratio(k)=Sp.b3_ratio;
end
end

% -------------------------------------------------------------------------
function x = restore(x, Ma, P, F)
%RESTORE  Newton correction onto {b1 = Ma} with projection onto (4).
x = F.project(x(:).', P);
s = P.s;
for it = 1:5
    g = (4/pi)*(cos(x)*s(:)) - Ma;
    if abs(g) < 1e-15, break; end
    G1 = -(4/pi) * sin(x) .* s;
    x  = F.project(x - (g./max(sum(G1.^2,2),1e-300)).*G1, P);
end
end
