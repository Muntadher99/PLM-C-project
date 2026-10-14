function d = device_parameters(Tvj)
%DEVICE_PARAMETERS  Semiconductor data of the device-level simulation and the
%   loss estimate (paper, Sections 6.2 and 6.3; reference [29]).
%
%   d = device_parameters()      typical values at Tvj = 25 C (used in the paper)
%   d = device_parameters(175)   typical values at Tvj = 175 C (bracketing case)
%
%   Device: Infineon IKW20N65ET7, 650 V / 20 A TRENCHSTOP IGBT7 with
%   co-packed diode, datasheet Rev. 1.20 (2023). Each device of the 800 V NPC
%   leg blocks 400 V, the datasheet test voltage of the switching energies.
%
%   Used by npc_loss_model, nonideal_switching_model, nonideal_simulink_runs
%   and nonideal_model_init.
%
%   How the values were obtained from the datasheet
%   * Conduction: least-squares line v = V0 + r*i through the typical output
%     characteristic (IGBT, VGE = 15 V) and the diode forward characteristic
%     over 0.5-6 A, the load-current range of Section 6.2. At 20 A, 25 C the
%     curves give 1.34 V / 1.64 V against 1.35 V / 1.65 V in the tables.
%   * Switching: E(20 A) from the switching table (VCC = 400 V, RG = 12 ohm,
%     VGE = 0/15 V); current dependence from the "switching energy vs
%     collector current" curve (175 C, the only one plotted), normalised to 1
%     at 20 A and tabulated at the knots E_I; below 4 A it is extended
%     linearly. With this shape the 25 C values give E(10 A) = 0.183 / 0.198 mJ
%     (Eon / Eoff) against 0.16 / 0.21 mJ in the 25 C table.
%   * Diode recovery: Erec = Qrr*VR/4, with Qrr(I) a power law through the two
%     tabulated points (20 A and 10 A).
%   * Tf: current fall time from the switching table; tail time Tt = 2*Tf
%     (not specified in the datasheet; Simulink IGBT block only).

if nargin < 1 || isempty(Tvj), Tvj = 25; end
d.Iref = 20;   d.Vref = 400;   d.kv = 1.0;
% normalised E(IC) shape, 175 C curve, value 1 at 20 A
d.E_I      = [4     5     6     8     10    12    15    20    25    30    40    50    60];
d.Eon_g    = [0.254 0.295 0.336 0.420 0.508 0.598 0.740 1.000 1.298 1.644 2.525 3.587 4.698];
d.Eoff_g   = [0.257 0.307 0.357 0.455 0.551 0.645 0.782 1.000 1.205 1.405 1.827 2.279 2.739];
d.ki       = 1.0;        % fallback exponent, used only if E_I is removed

switch Tvj
    case 25
        d.name     = 'Infineon IKW20N65ET7 (650 V, 20 A, IGBT7 + diode), typical, Tvj = 25 C';
        % conduction (fit 0.5-6 A; max deviation from the curve 0.04 / 0.06 V)
        d.Vce0     = 0.700;   d.rce = 0.0467;    % IGBT  : V, ohm
        d.Vf0      = 0.841;   d.rf  = 0.0655;    % diode : V, ohm
        % switching, IC = 20 A
        d.Eon_ref  = 0.36e-3;  d.Eoff_ref = 0.36e-3;          % J
        d.Erec_ref = 0.25 * 0.44e-6 * 400;                     % J (Qrr = 0.44 uC)
        d.ki_rec   = log(0.44/0.28)/log(2);                    % 0.65 (Qrr 0.28 uC at 10 A)
        d.Tf       = 20e-9;    % s, current fall time at 20 A
    case 175
        d.name     = 'Infineon IKW20N65ET7 (650 V, 20 A, IGBT7 + diode), typical, Tvj = 175 C';
        d.Vce0     = 0.482;   d.rce = 0.0749;
        d.Vf0      = 0.549;   d.rf  = 0.0797;
        d.Eon_ref  = 0.58e-3;  d.Eoff_ref = 0.61e-3;
        d.Erec_ref = 0.25 * 1.18e-6 * 400;                     % Qrr = 1.18 uC
        d.ki_rec   = log(1.18/0.85)/log(2);                    % 0.47
        d.Tf       = 75e-9;
    otherwise
        error('device_parameters:Tvj', 'datasheet values exist for Tvj = 25 or 175 C only');
end
d.Tt = 2*d.Tf;   % s, tail time (assumed, Simulink only)
d.Tvj = Tvj;
end
