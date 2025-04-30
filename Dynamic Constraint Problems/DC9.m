[file name]: DC9_DF.m
[file content begin]
function [f, c] = DC9_DF(x, tau, taut, nt, N)
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

% Time parameter calculation (aligned with DC9's t formula)
t = floor(tau/N/taut)/nt;

% Dynamic parameter calculation
G = sin(0.5 * pi * t);
w = 10 - abs(floor(10 * G)); % 波动频率参数

% Objective calculation
g = 1 + sum((x(2:end) - G).^2);
A = 0.02 * sin(w * pi * x(1)); % 震荡项
f1 = g * (x(1) + A);
f2 = g * (1 - x(1) + A);
f = [f1, f2];

% Constraint calculation (保持原约束逻辑)
phase_shift = sin(5 * pi * (f1 - f2 + 1)).^2;
c = -(f1 + f2 - phase_shift - G); % 原约束取反
end