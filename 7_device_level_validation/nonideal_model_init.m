%NONIDEAL_MODEL_INIT  Set every variable NPC_SHM.slx reads in the base
%   workspace, so that the model can be opened, updated and run by hand.
%   Values: datasheet devices with a 2 us dead time and a balanced DC link;
%   nonideal_simulink_runs overwrites them case by case.
d0 = device_parameters();
Ron   = d0.rce;     % IGBT on-state slope resistance (ohm)
Vf    = d0.Vce0;    % IGBT forward (knee) voltage (V)
Rond  = d0.rf;      % diode on-state slope resistance (ohm)
Vfd   = d0.Vf0;     % diode forward (knee) voltage (V)
Tf    = d0.Tf;      % IGBT current 10% fall time (s)
Tt    = d0.Tt;      % IGBT current tail time (s)
Rs    = 1e5;        % snubber resistance (ohm)
Cs    = Inf;        % snubber capacitance (Inf = purely resistive snubber)
Td    = 2e-6;       % dead time (s); 0 = ideal gating
C1    = 470e-6;  C2 = 470e-6;      % split DC-link capacitors (F)
VC1_0 = 400;     VC2_0 = 400;      % their initial voltages (V)
theta_deg = [4.323113 6.926852 21.464615 27.313345 32.230475 39.919268 44.218853 46.360690 50.066649 59.402712 61.055868 63.820445 67.679385 81.809510 86.375605];  % PLM-C, Ma = 0.80, seed 1, 20000 NFE
Ma    = 0.80;
fprintf('nonideal_model_init: model variables set (Td = %.1f us, device %s)\n', Td*1e6, d0.name);
