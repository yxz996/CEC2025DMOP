[file name]: DC9_DF_PF.m
[file content begin]
function PF = DC9_DF_PF(tau, taut, nt, N, sample_points)
% INPUT:
%       tau:        current generation counter
%       taut:       change frequency
%       nt:         change severity
%       N:          population size (for t calculation)
%       sample_points: number of PF points to generate
%
% OUTPUT:
%       PF:         Pareto front solutions

% 时间参数计算
t = floor(tau/N/taut)/nt;
G = sin(0.5 * pi * t);
w = 10 - abs(floor(10 * G));

% 生成基础PF点集
x1 = linspace(0, 1, sample_points)';
A = 0.02 * sin(w * pi * x1);
candidates = [x1 + A, 1 - x1 + A];

% 应用约束条件
phase_shift = sin(5 * pi * (candidates(:,1) - candidates(:,2) + 1).^2;
valid = candidates(:,1) + candidates(:,2) - phase_shift >= G; % 原约束取反后的条件
PF = candidates(valid, :);

% 非支配排序 (接口与NDSort1兼容)
[FrontNo, ~] = NDSort1(PF, 1);
PF = PF(FrontNo == 1, :);

% 时间偏移补偿 (原代码逻辑)
PF = PF + t; % 对应原代码中的pf+t
end