[file name]: DC8_DF.m
[file content begin]
function [f, c] = dcp5_DF(x, tau, taut, nt, N)
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

% Time parameter calculation (aligned with DC8's formula)
t = floor(tau/N/taut)/nt;

% Dynamic parameter calculation
G = 0.5 * abs(sin(0.5 * pi * t));

% Objective calculation
g = 1 + sum(((x(2:end)-0.5).^2 - cos(pi*(x(2:end)-0.5)) + 1).^2);
f1 = g * x(1);
f2 = g * (1 - x(1));
f = [f1, f2];

% Constraint calculation (maintain original sign conventions)
c1 = -((0.2+G)*f1^2 + f2 - 2) .* (0.7*f1^2 + f2 - 2.5); % 原约束取反
c2 = -(f1^2 + f2^2 - (0.6+G)^2); % 原约束取反
c = [c1, c2];
end