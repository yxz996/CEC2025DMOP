[file name]: DC7_DF.m
[file content begin]
function [f, c] = DC7_DF(x, tau, taut, nt, N)
% INPUT:
%       x:      variable vector (D-dimensional)
%       tau:    current generation counter
%       taut:   frequency of change
%       nt:     severity of change
%       N:      population size (for t calculation)
%
% OUTPUT:
%       f:      objective vector [f1, f2]
%       c:      constraint violation value

% Time parameter calculation (aligned with DC7's t formula)
t = floor(tau/N/taut)/nt;

% Dynamic parameter calculation
G = 2 * floor(10*abs(mod(t+1,2)-1 + 1e-4); % 保持原式精确实现

% Objective calculation
g = 1 + sum(((x(2:end)-0.5).^2 - cos(pi*(x(2:end)-0.5) + 1).^2); % 向量化计算
f1 = g * (x(1) + 0.25*sin(pi*x(1)));
f2 = g * (1 - x(1) + 0.25*sin(pi*x(1)));
f = [f1, f2];

% Constraint calculation (保持原约束逻辑)
theta = atan(f2/f1); % 处理f1=0的情况
if f1 == 0
    theta = pi/2; 
end
c11 = (f1)^2 + (f2)^2 - (1.5 + 0.4*sin(4*theta)^16)^2;
c12 = (f1)^2 + (f2)^2 - (1.3 - 0.45*sin(G*theta)^2)^2;
c = -c11*c12; % 原约束条件取反
end