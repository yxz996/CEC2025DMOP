[file name]: DC5_DF_PF.m
[file content begin]
function PF = dcp1_DF_PF(tau, taut, nt, N)
% INPUT:
%       tau:    current generation counter
%       taut:   frequency of change
%       nt:     severity of change
%       N:      number of points to sample
%
% OUTPUT:
%       PF:     Pareto front solutions

% Calculate time instant (consistent with DC5)
T0 = 50;
tau_tmp = max(tau + taut - (T0 + 1), 0);
t = floor(tau_tmp / taut) / nt;

% Generate reference points (aligned with DC5's GetOptimum)
G = abs(sin(0.5 * pi * t));
x = linspace(0, 1, N)';
P1 = x + G;
P2 = 1 - x + G;

% Apply constraints (mimic DC5's filtering)
Con = (2*sin(5*pi*(sin(-0.15*pi)*P2 + cos(-0.15*pi)*P1))).^6 ...
      - cos(-0.15*pi)*P2 + sin(-0.15*pi)*P1;
valid = Con <= 0;
PF = [P1(valid), P2(valid)];

% Non-dominated sorting (interface with NDSort1)
[FrontNo, ~] = NDSort1(PF, 1);
PF = PF(FrontNo == 1, :);
end