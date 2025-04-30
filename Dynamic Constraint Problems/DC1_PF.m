[file name]: DC1_DF_PF.m
[file content begin]
function PF = DC1_DF_PF(tau, taut, nt, N, sample_points)
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
G = cos(pi*t);

% 生成单位圆参考点
theta = linspace(0, pi/2, sample_points)';
V = [cos(theta), sin(theta)];

% 生成约束边界点
x_range = linspace(0, 1, sample_points)';
y_bound = 3 - G - exp(x_range) - 0.3*sin(4*pi*x_range);

% 合并候选解
candidates = [V; [x_range, y_bound]];

% 应用约束条件筛选
valid = (3 - G - candidates(:,1) - exp(candidates(:,1)) - 0.3*sin(4*pi*candidates(:,1)) - candidates(:,2) <= 0 & ...
        (4.1 - (1 + 0.3*candidates(:,1).^2 + candidates(:,1)) - 0.3*sin(4*pi*candidates(:,1)) - candidates(:,2) <= 0 & ...
        (candidates(:,1).^2 + candidates(:,2).^2 <= 1); % 单位圆约束

PF = candidates(valid, :);

% 非支配排序 (接口与NDSort1兼容)
[FrontNo, ~] = NDSort1(PF, 1);
PF = PF(FrontNo == 1, :);

% 边界修正 (移除不满足目标函数形状的解)
if ~isempty(PF)
    pf_norm = PF(:,1).^2 + (PF(:,2).^2)./ (1 - PF(:,1).^2);
    PF(pf_norm > 1.05, :) = [];
end
end