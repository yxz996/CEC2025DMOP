[file name]: DC2_DF.m
[file content begin]
function [f, c] = DC2_DF(x, tau, taut, nt, N)
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

% Time parameter calculation (aligned with DC2's formula)
t = floor(tau/N/taut)/nt;

% Dynamic parameter calculation
G = sin(0.5 * pi * t);

% Objective calculation
y = x(2:end) - G;
g = 1 + sum(y.^2 + sin(0.5*pi*y).^2);
f1 = g * (x(1) + 0.2*G*sin(pi*x(1)));
f2 = g * (1 - x(1) + 0.2*G*sin(pi*x(1)));
f = [f1, f2];

% Constraint calculation (保持原约束逻辑)
c1 = -( (f1 + 2*f2 - 1) .* (f1 + 0.5*f2 - 0.5) ); % 原约束取反
c2 = f1.^2 + f2.^2 - 1.4^2;
c = [c1, c2];
end