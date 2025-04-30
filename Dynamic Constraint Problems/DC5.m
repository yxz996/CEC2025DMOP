[file name]: DC5_DF.m
[file content begin]
function [f, c] = DC5_DF(x, tau, taut, nt)
% INPUT:
%       x:      variable vector
%       tau:    current generation counter
%       taut:   frequency of change
%       nt:     severity of change
%
% OUTPUT:
%       f:      objective vector
%       c:      constraint violation value

% Calculate time instant (aligned with DC5's t calculation)
T0 = 50;
tau_tmp = max(tau + taut - (T0 + 1), 0);
t = floor(tau_tmp / taut) / nt;

% Problem parameters
G = abs(sin(0.5 * pi * t));

% Objective calculation
g = 1 + sum((x(2:end) - G).^2);
f1 = g * x(1) + G;
f2 = g * (1 - x(1)) + G;
f = [f1, f2];

% Constraint calculation (aligned with DC5's Constraint function)
c = (2 * sin(5 * pi * (sin(-0.15*pi)*f2 + cos(-0.15*pi)*f1))).^6 ...
    - cos(-0.15*pi)*f2 + sin(-0.15*pi)*f1;
end