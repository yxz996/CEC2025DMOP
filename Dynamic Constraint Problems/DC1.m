[file name]: DC1_DF.m
[file content begin]
function [f, c] = DC1_DF(x, tau, taut, nt, N)
% INPUT:
%       x:      variable vector (D-dimensional)
%       tau:    current generation counter
%       taut:   change frequency
%       nt:     change severity
%       N:      population size (for t calculation)
%
% OUTPUT:
%       f:      objective vector [f1, f2]
%       c:      constraint violation value

% Time parameter calculation (aligned with DC1's t formula)
t = floor(tau/N/taut)/nt;

% Objective calculation
g = 1 + sum(x(2:end).^2);
f1 = g * x(1);
f2 = g * sqrt(1 - x(1)^2);
f = [f1, f2];

% Constraint calculation (保持原约束逻辑)
G = cos(pi*t);
c1 = (3 - G - exp(f1) - 0.3*sin(4*pi*f1) - f2) .* ...
     (4.1 - (1 + 0.3*f1^2 + f1) - 0.3*sin(4*pi*f1) - f2);
c = c1; % 原约束条件直接传递
end