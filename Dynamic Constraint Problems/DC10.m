[file name]: DC10_DF.m
[file content begin]
function [f, c] = dcp8_DF(x, tau, taut, nt, N)
% INPUT:
%       x:      variable vector (D-dimensional)
%       tau:    current generation counter
%       taut:   change frequency
%       nt:     change severity
%       N:      population size (for t calculation)
%
% OUTPUT:
%       f:      objective vector [f1, f2]
%       c:      constraint violations [c1, c2]

% Time parameter calculation (aligned with DC10's t formula)
t = floor(tau/N/taut)/nt;

% Dynamic parameter calculation
G = sin(0.5 * pi * t);

% Objective calculation
y = x(2:end) - G;
g = 1 + sum((abs(G)*y.^2 - cos(pi*y) + 1).^2);
f1 = g * x(1);
f2 = g * (1 - x(1));
f = [f1, f2];

% Constraint calculation (保持原约束逻辑)
c11 = (f1)^1.5 + (f2)^1.5 - 1.2^1.5;
c12 = sqrt(f1) + sqrt(f2) - (0.95 + 0.5*abs(G));
c21 = 0.8*f1 + f2 - (2.5 + 0.08*sin(2*pi*(f2 - f1)));
c22 = (0.93 + abs(G)/3)*f1 + f2 - (2.7 + abs(G)/2 + 0.08*sin(2*pi*(f2 - f1)));

c = -[c11*c12, c21*c22]; % 原约束取反
end