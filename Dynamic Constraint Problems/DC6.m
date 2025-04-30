[file name]: DC6_DF.m
[file content begin]
function [f, c] = dcp2_DF(x, tau, taut, nt, N)
% INPUT:
%       x:      variable vector
%       tau:    current generation counter
%       taut:   frequency of change
%       nt:     severity of change
%       N:      population size (for t calculation)
%
% OUTPUT:
%       f:      objective vector
%       c:      constraint violations

% Time parameter calculation (aligned with DC6's t formula)
t = floor(tau/N/taut)/nt;

% Dynamic parameter
G = sin(0.5 * pi * t);

% Objective calculation
g = 1 + sum((x(2:end)-G).^2 + sin(0.5*pi*(x(2:end)-G)).^2;
f1 = g * (x(1) + 0.25*G*sin(pi*x(1)));
f2 = g * (1 - x(1) + 0.25*G*sin(pi*x(1)));
f = [f1, f2];

% Constraint calculation (reverse sign for <=0 convention)
c1 = -(4*f1 + f2 - 1).*(0.3*f1 + f2 - 0.3);
c2 = (1.85 - f1 - f2 - (0.3*sin(3*pi*(f2-f1))).^2).*(f1 + f2 - 1.3);
c = [c1, c2];
end