[file name]: DC3_DF.m
[file content begin]
function [f, c] = DC3_DF(x, tau, taut, nt, N)
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

% Time parameter calculation (aligned with DC3's t formula)
t = floor(tau/N/taut)/nt;

% Dynamic parameter calculation
G = sin(0.5 * pi * t);

% Objective calculation
g = 1 + sum((x(2:end) - G).^2);
f1 = g * x(1);
f2 = g * sqrt(1 - x(1)^2);
f = [f1, f2];

% Constraint calculation (保持原约束逻辑)
H = 1;
if G < 0
    H = -1;
end
theta = atan((f2./f1).^H);
ll = cos(5*theta.^4).^6;
c = -(1.1 - (f1./(0.9 + (0.1+0.7*abs(G))*ll)).^2 - (f2./(0.9 + (0.8-0.7*abs(G))*ll)).^2);
end