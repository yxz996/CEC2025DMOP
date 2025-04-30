[file name]: DC10_DF_PF.m
[file content begin]
function PF = DC10_DF_PF(tau, taut, nt, N, sample_points)
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

% 生成均匀参考点
theta = linspace(0, pi/2, sample_points)';
X = [cos(theta), sin(theta)];
pf1 = 1.2 * X ./ repmat((sum(X.^1.5, 2)).^(2/3), 1, 2);

% 生成直线段PF
x1 = linspace(0, 1, sample_points)';
pf_linear = [x1, 1 - x1];

% 合并候选解
candidates = [pf1; pf_linear];

% 应用约束条件筛选
valid = (candidates(:,1).^1.5 + candidates(:,2).^1.5 <= 1.2^1.5) & ...
        (sqrt(candidates(:,1)) + sqrt(candidates(:,2)) >= 0.95 + 0.5*abs(G)) & ...
        (0.8*candidates(:,1) + candidates(:,2) <= 2.5 + 0.08*sin(2*pi*(candidates(:,2)-candidates(:,1)))) & ...
        ((0.93 + abs(G)/3)*candidates(:,1) + candidates(:,2) <= 2.7 + abs(G)/2 + 0.08*sin(2*pi*(candidates(:,2)-candidates(:,1))));

PF = candidates(valid, :);

% 非支配排序 (接口与NDSort1兼容)
[FrontNo, ~] = NDSort1(PF, 1);
PF = PF(FrontNo == 1, :);
end