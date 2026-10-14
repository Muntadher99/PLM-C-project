function [Sa1,Sa2,Sa3,Sa4,Sb1,Sb2,Sb3,Sb4,Sc1,Sc2,Sc3,Sc4] = shmpwm3ph_15_deadtime(t, alpha, Td)
%SHMPWM3PH_15_DEADTIME  Gate logic of the MATLAB Function block in
%   NPC_SHM.slx: twelve NPC gate signals from the quarter-wave angles, with
%   dead time Td on every complementary pair.
%
%   This file is a readable copy of the code embedded in the model; the model
%   does not call it from the path.
%
%   Inputs : t      simulation time (s)
%            alpha  switching angles (rad), any order
%            Td     dead time (s); Td = 0 gives ideal gating
%   Outputs: Sx1..Sx4 gates of leg x (x = a, b, c), 1 = on, 0 = off;
%            Sx1 upper outer, Sx2 upper inner, Sx3 lower inner, Sx4 lower
%            outer (S_a1..S_c4 of Fig. 1)
%#codegen
alpha = sort(alpha(:)).';
f = 50;
w = 2*pi*f;
[Sa1,Sa2,Sa3,Sa4] = leg_gates(w, t, 0,        alpha, Td);
[Sb1,Sb2,Sb3,Sb4] = leg_gates(w, t, 2*pi/3,   alpha, Td);
[Sc1,Sc2,Sc3,Sc4] = leg_gates(w, t, 4*pi/3,   alpha, Td);
end

function [S1,S2,S3,S4] = leg_gates(w, t, phi, alpha, Td)
% Commanded level now and Td seconds ago.
L0 = phase_level(w*t - phi,        alpha);
L1 = phase_level(w*(t - Td) - phi, alpha);
[P0, N0] = level_to_PN(L0);
[P1, N1] = level_to_PN(L1);
% Level +1 -> S1,S2 on;  0 -> S2,S3 on;  -1 -> S3,S4 on, i.e.
% S1 = P, S2 = not N, S3 = not P, S4 = N. A switch turns on only after it
% has been commanded on for Td and turns off at once, so the complementary
% pairs S1/S3 and S2/S4 are both off for Td at every commutation.
S1 = P0 * P1;
S3 = (1 - P0) * (1 - P1);
S4 = N0 * N1;
S2 = (1 - N0) * (1 - N1);
end

function L = phase_level(theta, alpha)
theta = mod(theta, 2*pi);
if theta <= pi/2
    L = first_quarter_level(theta, alpha);
elseif theta <= pi
    L = first_quarter_level(pi - theta, alpha);
elseif theta <= 3*pi/2
    L = -first_quarter_level(theta - pi, alpha);
else
    L = -first_quarter_level(2*pi - theta, alpha);
end
end

function L = first_quarter_level(theta, alpha)
% 0 at the start of the quarter, toggling 0/1 at every angle.
L = 0;
for k = 1:length(alpha)
    if theta >= alpha(k)
        L = 1 - L;
    else
        break;
    end
end
end

function [P, N] = level_to_PN(L)
if L > 0.5
    P = 1; N = 0;
elseif L < -0.5
    P = 0; N = 1;
else
    P = 0; N = 0;
end
end
