function A = audit_angle_sets(Ma, TH_deg, labels, outcsv)
%AUDIT_ANGLE_SETS  Independent audit of switching-angle sets: every quantity
%   a reader needs to verify a design from the angles alone (paper, Section
%   5.7; Reviewer 1 comment 2, Reviewer 2 comment 2).
%
%   A = audit_angle_sets(Ma, TH_deg)                % Ma: n x 1, TH_deg: n x N
%   A = audit_angle_sets(Ma, TH_deg, labels, 'witnesses_audit.csv')
%
%   Recomputed for every row from Eq. (2) and the Table 1 limits:
%     b1 and |b1 - Ma|, THD49, WTHD49, the sixteen loadings
%     q_n = |b_n| / (L_n |b1|) for n = 5, 7, ..., 49, q_max and the binding
%     orders (q_n within 1e-4 of q_max), the smallest gap theta_{i+1} -
%     theta_i, the margins theta_1 - theta_min and theta_max - theta_N, the
%     worst constraint violation, feasibility (|b1 - Ma| < 1e-8 and
%     violation < 1e-9) and compliance (q_max <= 1 + 1e-4).
%   Nothing is taken from the solver that produced the angles.
%
%   Output CSV (one row per angle set): label, Ma, theta01_deg..thetaNN_deg,
%   b1, b1_err, THD49_pct, WTHD49_pct, q5..q49, q_max, binding_orders,
%   min_gap_deg, first_margin_deg, last_margin_deg, viol_rad, feasible,
%   compliant. The struct A carries the same columns.

if nargin < 3 || isempty(labels), labels = repmat({''}, numel(Ma), 1); end
if ischar(labels), labels = repmat({labels}, numel(Ma), 1); end
if nargin < 4, outcsv = ''; end
Ma = Ma(:);  n = numel(Ma);  N = size(TH_deg, 2);
F = plm_core();  P = F.defaults(N, 49);
TH = sort(deg2rad(TH_deg), 2);

B1 = F.b1(TH, P);  B = F.bn(TH, P);
q  = abs(B) ./ (abs(B1) .* P.LIM);
A.labels  = labels(:);
A.Ma      = Ma;
A.theta   = rad2deg(TH);
A.b1      = B1;
A.b1_err  = abs(B1 - Ma);
A.THD49   = 100 * sqrt(sum(B.^2, 2)) ./ abs(B1);
A.WTHD49  = 100 * sqrt(sum((B ./ P.ORD).^2, 2)) ./ abs(B1);
A.q       = q;
A.orders  = P.ORD;
[A.qmax, jm] = max(q, [], 2); %#ok<ASGLU>
A.binding = cell(n,1);
for k = 1:n
    A.binding{k} = strjoin(arrayfun(@num2str, P.ORD(q(k,:) >= A.qmax(k) - 1e-4), ...
                           'UniformOutput', false), ' ');
end
A.min_gap_deg   = rad2deg(min(diff(TH, 1, 2), [], 2));
A.first_margin  = rad2deg(TH(:,1) - P.lb);
A.last_margin   = rad2deg(P.ub - TH(:,end));
gv = max([zeros(n,1), P.minGap - diff(TH,1,2)], [], 2);
bv = max(max(0, P.lb - TH(:,1)), max(0, TH(:,end) - P.ub));
A.viol      = max(gv, bv);
A.compliant = A.qmax <= 1 + P.COMP_TOL;
A.feasible  = (A.b1_err < P.FEAS_B1) & (A.viol < 1e-9);

if ~isempty(outcsv)
    fid = fopen(outcsv, 'w');
    fprintf(fid, 'label,Ma');
    fprintf(fid, ',theta%02d_deg', 1:N);
    fprintf(fid, ',b1,b1_err,THD49_pct,WTHD49_pct');
    fprintf(fid, ',q%d', P.ORD);
    fprintf(fid, ',q_max,binding_orders,min_gap_deg,first_margin_deg,last_margin_deg,viol_rad,feasible,compliant\n');
    for k = 1:n
        fprintf(fid, '%s,%.4f', A.labels{k}, Ma(k));
        fprintf(fid, ',%.10f', A.theta(k,:));
        fprintf(fid, ',%.15g,%.3e,%.6f,%.6f', A.b1(k), A.b1_err(k), A.THD49(k), A.WTHD49(k));
        fprintf(fid, ',%.6f', q(k,:));
        fprintf(fid, ',%.6f,%s,%.6f,%.6f,%.6f,%.3e,%d,%d\n', A.qmax(k), A.binding{k}, ...
            A.min_gap_deg(k), A.first_margin(k), A.last_margin(k), A.viol(k), ...
            A.feasible(k), A.compliant(k));
    end
    fclose(fid);
    fprintf('audit_angle_sets: wrote %d audited angle sets to %s\n', n, outcsv);
end
end
