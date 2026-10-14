function W = nonideal_switching_model(theta_deg, Ma, NI, cfg, dev)
%NONIDEAL_SWITCHING_MODEL  Switching-function model of the three-phase NPC
%   inverter with converter non-idealities, for one quarter-wave angle set
%   (paper, Section 6.2: "a switching-function model reproduces the simulated
%   THD49 to a median of 0.009 pp"). No Simulink required.
%
%   W = nonideal_switching_model(theta_deg, Ma)                     % ideal
%   W = nonideal_switching_model(theta_deg, Ma, struct('td',2e-6))  % 2 us dead time
%
%   NI fields (any subset; absent = ideal)
%     drops   true: device forward drops (Vce0 + rce|i|, Vf0 + rf|i|)
%     td      dead time (s): every commutation that needs an IGBT turn-on is
%             delayed by td; commutations forced by a turn-off are not
%     eps_dc  static neutral-point offset, fraction of Vdc/2
%     Ctot    total DC-link capacitance C1 + C2 (F): adds the neutral-point
%             ripple v1(t) = Vdc/2 + (1/Ctot) * integral(i_NP)
%     tq      timer resolution (s): commanded instants rounded to multiples of tq
%     iters   fixed-point iterations current <-> voltage [4]
%   cfg: Vdc, f1, R, Lload, M (as npc_loss_model);  dev: as npc_loss_model.
%
%   W fields: THD49 (%), maxq, b1 (recovered fundamental / (Vdc/2)),
%   b1_err_pct, THD50_all (%), even_max (% of V1, largest even order <= 50),
%   np_ripple_pp (V), dthd (pp, against the analytic THD49 of the same
%   angles), I1_pk (A).

if nargin < 3 || isempty(NI),  NI = struct(); end
if nargin < 4 || isempty(cfg), cfg = struct(); end
if nargin < 5 || isempty(dev), dev = struct(); end
cfg = fillc(cfg, 'Vdc', 800);  cfg = fillc(cfg, 'f1', 50);
cfg = fillc(cfg, 'R', 80);     cfg = fillc(cfg, 'Lload', 25e-3);
cfg = fillc(cfg, 'M', 2^17);
NI  = fillc(NI, 'drops', false);  NI = fillc(NI, 'td', 0);
NI  = fillc(NI, 'eps_dc', 0);     NI = fillc(NI, 'Ctot', Inf);
NI  = fillc(NI, 'tq', 0);         NI = fillc(NI, 'iters', 4);
dev = dev_defaults(dev);

M = cfg.M;  f1 = cfg.f1;  dt = 1/(f1*M);  w1 = 2*pi*f1;  E2 = cfg.Vdc/2;
ang = sort(deg2rad(theta_deg(:)).');
if NI.tq > 0
    d = 2*pi*f1*NI.tq;  ang = round(ang/d)*d;          % timer quantisation
end
th = 2*pi*(0:M-1).'/M;
S0 = [state(th, ang), state(mod(th-2*pi/3,2*pi), ang), state(mod(th-4*pi/3,2*pi), ang)];
nd = round(NI.td/dt);                                    % dead time in samples
h  = [0:M/2, -(M/2-1):-1].';
Z  = cfg.R + 1i*h*w1*cfg.Lload;

% ---- fixed-point iteration: current -> delays/drops/ripple -> voltage -> current
S = S0;  v1 = E2*(1+NI.eps_dc)*ones(M,1);  v2 = cfg.Vdc - v1;
I = zeros(M,3);
for p = 1:3, I(:,p) = current(S0(:,p)*E2, S0*E2, Z, p); end
for it = 1:max(1, NI.iters)
    for p = 1:3
        S(:,p) = apply_deadtime(S0(:,p), I(:,p), nd);
    end
    V0 = zeros(M,3);
    for p = 1:3
        s = S(:,p);  i = I(:,p);  ai = abs(i);  sg = sign(i);
        v = (s==1).*v1 - (s==-1).*v2;
        if NI.drops
            vT = dev.Vce0 + dev.rce*ai;  vD = dev.Vf0 + dev.rf*ai;
            pos = i > 0;
            drop = (s==1).*( pos.*2.*vT + ~pos.*2.*vD) + ...
                   (s==0).*( vT + vD ) + ...
                   (s==-1).*(~pos.*2.*vT + pos.*2.*vD);
            v = v - sg.*drop;
        end
        V0(:,p) = v;
    end
    for p = 1:3, I(:,p) = current(V0(:,p), V0, Z, p); end
    if isfinite(NI.Ctot)
        inp = sum(I .* (S == 0), 2);                     % current drawn from NP
        inp = inp - mean(inp);                            % periodic component
        rip = cumsum(inp)*dt / NI.Ctot;  rip = rip - mean(rip);
        v1 = E2*(1+NI.eps_dc) + rip;  v2 = cfg.Vdc - v1;
    end
end

% ---- spectra (definitions of the paper) ------------------------------------
vab = V0(:,1) - V0(:,2);
Vh  = abs(fft(vab)) / (M/2);  V1 = Vh(2);
ORD = 5:2:49;  ORD = ORD(mod(ORD,3) ~= 0);
Lim = lim_of(ORD);
W.THD49 = 100*sqrt(sum(Vh(ORD+1).^2))/V1;
W.maxq  = max((Vh(ORD+1).'/V1)./Lim);
o2 = 2:50; o2 = o2(mod(o2,3) ~= 0);
W.THD50_all = 100*sqrt(sum(Vh(o2+1).^2))/V1;
ev = 2:2:50;  W.even_max = 100*max(Vh(ev+1))/V1;
W.b1 = V1/(sqrt(3)*E2);
W.b1_err_pct = 100*(W.b1 - Ma)/Ma;
W.np_ripple_pp = max(v1) - min(v1);
F = plm_core();  P = F.defaults(numel(ang), 49);
R = F.report(deg2rad(theta_deg(:)).', Ma, P);
W.dthd = W.THD49 - 100*R.THD49;
W.I1_pk = abs(fft(I(:,1)))/(M/2);  W.I1_pk = W.I1_pk(2);
end

% =========================================================================
function s = state(th, ang)
%STATE  Pole state of one phase at electrical angles th.
q = mod(th, 2*pi);  sgn = ones(size(q));  sgn(q >= pi) = -1;
q = mod(q, pi);  q(q > pi/2) = pi - q(q > pi/2);
s = sgn .* mod(sum(q >= ang, 2), 2);
end

function i = current(vp, V0, Z, p) %#ok<INUSL>
%CURRENT  Periodic steady-state phase current of the star RL load.
vN = V0(:,p) - mean(V0, 2);
Ih = fft(vN) ./ Z;  Ih(1) = 0;
i = real(ifft(Ih));
end

function s = apply_deadtime(s0, i, nd)
%APPLY_DEADTIME  Delay by nd samples every commanded transition whose
%   completion needs an IGBT turn-on:
%   0->+1 with i>0, +1->0 with i<0, 0->-1 with i<0, -1->0 with i>0.
%   A pulse shorter than the dead time is swallowed.
s = s0;
if nd <= 0, return; end
M = numel(s0);
prev = s0([end 1:end-1]);
k = find(s0 ~= prev);
for j = k(:).'
    a = prev(j);  b = s0(j);  ij = i(j);
    need = (a==0 && b==1 && ij>0) || (a==1 && b==0 && ij<0) || ...
           (a==0 && b==-1 && ij<0) || (a==-1 && b==0 && ij>0);
    if need
        idx = mod((j:j+nd-1)-1, M) + 1;
        s(idx) = a;
    end
end
end

function L = lim_of(orders)
%LIM_OF  Table 1 limits as fractions of the fundamental.
L = zeros(size(orders));
for k = 1:numel(orders)
    n = orders(k);
    switch n
        case 5, L(k)=0.06; case 7, L(k)=0.05; case 11, L(k)=0.035; case 13, L(k)=0.03;
        case 17, L(k)=0.02; case 19, L(k)=0.015; case 23, L(k)=0.015; case 25, L(k)=0.015;
        otherwise, L(k) = (0.2 + 32.5/n)/100;
    end
end
end

function c = fillc(c, f, v)
if ~isfield(c, f) || isempty(c.(f)), c.(f) = v; end
end

function d = dev_defaults(d)
%DEV_DEFAULTS  Device data from device_parameters; fields passed in override it.
D0 = device_parameters();
f = fieldnames(D0);
for k = 1:numel(f), d = fillc(d, f{k}, D0.(f{k})); end
end
