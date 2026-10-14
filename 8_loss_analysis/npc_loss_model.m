function L = npc_loss_model(sw, Ma, cfg, dev)
%NPC_LOSS_MODEL  Semiconductor loss estimate of a three-phase three-level NPC
%   inverter driven by programmed (SHM) angles or by carrier PWM (paper,
%   Section 6.3; analytical-numerical method of [30]).
%
%   L = npc_loss_model(struct('theta_deg', th), Ma)            % SHM angle set
%   L = npc_loss_model(struct('pwm','pd','mf',63,'thi',true), Ma)  % PD carrier PWM
%   L = npc_loss_model(sw, Ma, cfg, dev)                       % other converter/device
%
%   Model
%   * Pole-voltage state s in {+1, 0, -1} on cfg.M samples per period.
%     SHM: quarter-wave angles with alternating 0->+1 / +1->0 transitions,
%     the waveform of Eq. (2). PWM: phase-disposition carriers [31] with
%     optional min-max injection (sw.thi = true).
%   * Load: balanced star RL (R, Lload per phase); the periodic steady-state
%     current is computed per harmonic, I_h = V_h / (R + j h w1 L), by FFT.
%   * Conduction paths of an NPC leg (i > 0 flows out of the leg):
%       s=+1: i>0 T1,T2   i<0 D1,D2        s=0: i>0 D5,T2   i<0 T3,D6
%       s=-1: i<0 T3,T4   i>0 D3,D4
%     IGBT v = Vce0 + rce|i|, diode v = Vf0 + rf|i|.
%   * Each commutation (blocking Vdc/2) is charged the turn-on, turn-off and
%     recovery energies of the devices involved:
%       0 ->+1: i>0 Eon(T1)+Erec(D5)    i<0 Eoff(T3)
%      +1-> 0 : i>0 Eoff(T1)            i<0 Eon(T3)+Erec(D1)
%       0 ->-1: i<0 Eon(T4)+Erec(D6)    i>0 Eoff(T2)
%      -1-> 0 : i<0 Eoff(T4)            i>0 Eon(T2)+Erec(D4)
%     E(i,v) = Eref * g(|i|) * (v/Vref)^kv, g interpolated from the datasheet
%     curve (device_parameters); diode recovery as a power law in current.
%   * Device parameters at one junction temperature (25 C unless dev.Tvj = 175).
%
%   cfg fields: Vdc [800 V], f1 [50 Hz], R [80 ohm], Lload [25 mH], M [2^17]
%
%   L fields (three-phase converter, powers in W): P_cond_T, P_cond_D,
%   P_sw_T (Eon + Eoff), P_rec_D, P_loss, P_out, eta (%), fsw_dev (turn-ons per
%   second of the most switched device), THD49 (%), THD50_all (%), maxq,
%   I1_pk (A), ITHD (%), b1 (recovered fundamental in units of Vdc/2).

if nargin < 3 || isempty(cfg), cfg = struct(); end
if nargin < 4 || isempty(dev), dev = struct(); end
cfg = fillc(cfg, 'Vdc', 800);          % V, Section 6.2
cfg = fillc(cfg, 'f1', 50);            % Hz
cfg = fillc(cfg, 'R', 80);             % ohm per phase, Section 6.2
cfg = fillc(cfg, 'Lload', 25e-3);      % H per phase, Section 6.2
cfg = fillc(cfg, 'M', 2^17);           % samples per fundamental period
dev = fill_dev(dev);

M  = cfg.M;  w1 = 2*pi*cfg.f1;
th = 2*pi*(0:M-1).'/M;                 % electrical angle of each sample

% ---- switching states of the three phases -------------------------------
if isfield(sw, 'theta_deg') || isfield(sw, 'theta')
    if isfield(sw, 'theta_deg'), ang = sort(deg2rad(sw.theta_deg(:)).'); else, ang = sort(sw.theta(:).'); end
    sA = shm_state(th, ang);
    sB = shm_state(mod(th - 2*pi/3, 2*pi), ang);
    sC = shm_state(mod(th - 4*pi/3, 2*pi), ang);
    kind = 'shm';
else
    if ~isfield(sw, 'mf'), error('npc_loss_model:mf', 'PWM needs sw.mf (carrier/fundamental ratio)'); end
    thi = isfield(sw, 'thi') && sw.thi;
    m = [sin(th), sin(th - 2*pi/3), sin(th - 4*pi/3)] * Ma;
    if thi, m = m - (max(m,[],2) + min(m,[],2))/2; end
    ph = 0; if isfield(sw, 'carrier_phase'), ph = sw.carrier_phase; end
    cu = tri01(sw.mf * th + ph);                    % upper carrier in [0,1]
    cl = cu - 1;                                    % lower carrier in [-1,0]
    S = zeros(M, 3);
    for p = 1:3
        S(:,p) = (m(:,p) >= cu) - (m(:,p) < cl);
    end
    sA = S(:,1); sB = S(:,2); sC = S(:,3);
    kind = 'pwm';
end

% ---- voltages and periodic steady-state current --------------------------
E2 = cfg.Vdc/2;
v0 = [sA sB sC] * E2;                               % pole voltages
vN = v0 - mean(v0, 2);                              % phase voltages of the star load
h  = [0:M/2, -(M/2-1):-1].';                        % harmonic index of each FFT bin
Z  = cfg.R + 1i*h*w1*cfg.Lload;
IA = fft(vN(:,1)) ./ Z;  IA(1) = 0;
iA = real(ifft(IA));

% ---- losses of leg A -----------------------------------------------------
dt = 1/(cfg.f1*M);
pos = iA > 0;  ai = abs(iA);
pT = dev.Vce0*ai + dev.rce*ai.^2;
pD = dev.Vf0*ai  + dev.rf*ai.^2;
nT = 2*((sA==1 & pos) | (sA==-1 & ~pos)) + 1*((sA==0));   % IGBTs conducting
nD = 2*((sA==1 & ~pos) | (sA==-1 & pos)) + 1*((sA==0));   % diodes conducting
Pc_T = sum(nT .* pT) * dt * cfg.f1;
Pc_D = sum(nD .* pD) * dt * cfg.f1;

s_prev = sA([end 1:end-1]);
ev = find(sA ~= s_prev);
Esw = 0;  Erec = 0;  nOnT1 = 0;  nOnT4 = 0;
for k = ev(:).'
    a = s_prev(k);  b = sA(k);
    i = 0.5*(iA(k) + iA(max(k-1,1)));
    Eon  = escale(dev.Eon_ref,  abs(i), E2, dev, 'Eon_g');
    Eoff = escale(dev.Eoff_ref, abs(i), E2, dev, 'Eoff_g');
    Er   = escale(dev.Erec_ref, abs(i), E2, dev, '');
    if     a == 0 && b == 1           % T1 on / T3 off
        nOnT1 = nOnT1 + 1;
        if i > 0, Esw = Esw + Eon;  Erec = Erec + Er; else, Esw = Esw + Eoff; end
    elseif a == 1 && b == 0           % T1 off / T3 on
        if i > 0, Esw = Esw + Eoff; else, Esw = Esw + Eon; Erec = Erec + Er; end
    elseif a == 0 && b == -1          % T4 on / T2 off
        nOnT4 = nOnT4 + 1;
        if i < 0, Esw = Esw + Eon;  Erec = Erec + Er; else, Esw = Esw + Eoff; end
    elseif a == -1 && b == 0          % T4 off / T2 on
        if i < 0, Esw = Esw + Eoff; else, Esw = Esw + Eon; Erec = Erec + Er; end
    else                               % +1 <-> -1 in one sample: two commutations
        Esw = Esw + Eon + Eoff;  Erec = Erec + Er;
    end
end
Psw_T = Esw * cfg.f1;  Prec_D = Erec * cfg.f1;

% ---- three-phase totals ---------------------------------------------------
L.kind     = kind;
L.P_cond_T = 3*Pc_T;   L.P_cond_D = 3*Pc_D;
L.P_sw_T   = 3*Psw_T;  L.P_rec_D  = 3*Prec_D;
L.P_loss   = L.P_cond_T + L.P_cond_D + L.P_sw_T + L.P_rec_D;
L.P_out    = 3 * cfg.R * mean(iA.^2);
L.eta      = 100 * L.P_out / (L.P_out + L.P_loss);
L.fsw_dev  = max(nOnT1, nOnT4) * cfg.f1;

% ---- harmonic content of the line voltage (definitions of the paper) ------
vab = v0(:,1) - v0(:,2);
Vh  = abs(fft(vab)) / (M/2);
V1  = Vh(2);
ORD = 5:2:49;  ORD = ORD(mod(ORD,3) ~= 0);
Lim = grid_lim(ORD);
L.THD49 = 100 * sqrt(sum(Vh(ORD+1).^2)) / V1;
o2 = 2:50;  o2 = o2(mod(o2,3) ~= 0);
L.THD50_all = 100 * sqrt(sum(Vh(o2+1).^2)) / V1;
L.maxq  = max((Vh(ORD+1).'/V1) ./ Lim);
Ih = abs(IA) / (M/2);
L.I1_pk = Ih(2);
L.ITHD  = 100 * sqrt(sum(Ih(3:M/2).^2)) / Ih(2);
L.b1    = V1 / (sqrt(3) * E2);          % recovered fundamental, units of Vdc/2
end

% =========================================================================
function s = shm_state(th, ang)
%SHM_STATE  State of one phase at electrical angles th (0..2pi) for the
%   quarter-wave angles ang.
q = mod(th, 2*pi);
sgn = ones(size(q));  sgn(q >= pi) = -1;
q = mod(q, pi);  q(q > pi/2) = pi - q(q > pi/2);          % fold to [0, pi/2]
cnt = sum(q >= ang, 2);                                    % transitions passed
s = sgn .* mod(cnt, 2);
end

function c = tri01(x)
%TRI01  Unit triangle carrier in [0,1] with period 2*pi, value 0 at x = 0.
c = abs(mod(x/(2*pi), 1)*2 - 1);
c = 1 - c;
end

function E = escale(Eref, i, v, dev, gfield)
%ESCALE  Switching energy at current i and blocked voltage v. The current
%   dependence is the tabulated datasheet curve (dev.E_I, dev.(gfield)),
%   extended linearly below its first point and clipped at zero, or a power
%   law (|i|/Iref)^k (k = dev.ki_rec for the diode recovery energy).
if ~isempty(gfield) && isfield(dev, 'E_I') && isfield(dev, gfield) && ~isempty(dev.(gfield))
    g = max(interp1(dev.E_I, dev.(gfield), i, 'linear', 'extrap'), 0);
else
    k = dev.ki;  if isempty(gfield) && isfield(dev, 'ki_rec'), k = dev.ki_rec; end
    g = (i/dev.Iref)^k;
end
E = Eref * g * (v/dev.Vref)^dev.kv;
end

function L = grid_lim(orders)
%GRID_LIM  Table 1 limits as fractions of the fundamental.
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

function d = fill_dev(d)
%FILL_DEV  Device data from device_parameters; fields passed in override it.
T = 25;  if isfield(d, 'Tvj') && ~isempty(d.Tvj), T = d.Tvj; end
D0 = device_parameters(T);
f = fieldnames(D0);
for k = 1:numel(f), d = fillc(d, f{k}, D0.(f{k})); end
end
