[file name]: DC4_DF.m
[file content begin]
function [f, c] = DC4_DF(x, tau, taut, nt, N)
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

% Time parameter calculation (aligned with DC4's formula)
t = floor(tau/N/taut)/nt;

% Dynamic parameter calculation
G = sin(0.5 * pi * t);
W = 0.5 * pi * t;

% Objective calculation
g = 1 + sum((x(2:end) - G).^2);
f1 = g * x(1);
f2 = g * (1 - x(1));
f = [f1, f2];

% Constraint calculation (保持原约束逻辑)
A = (sin(W)-cos(W))/(sin(W)+cos(W));
B = -2.2*(1-cos(W));

% 第一组约束
c11 = A*f1 - f2 + (B + 1.3)/(sin(W)+cos(W));
c12 = A*f1 - f2 + (B + 1.8)/(sin(W)+cos(W));
c1 = c11 .* c12;

% 第二组约束
c31 = A*f1 - f2 + (B + 2.6)/(sin(W)+cos(W));
c32 = A*f1 - f2 + (B + 3.1)/(sin(W)+cos(W));
c3 = c31 .* c32;

c = -[c1, c3]; % 原约束取反
end