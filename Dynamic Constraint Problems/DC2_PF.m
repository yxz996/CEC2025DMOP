[file name]: DC2_DF_PF.m
[file content begin]
function PF = DC2_DF_PF(tau, taut, nt, N, sample_points)
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

% 生成基础PF点集（分段线性）
x1 = linspace(0, 1, sample_points)';
pf_linear = [x1 + 0.2*G*sin(pi*x1), 1 - x1 + 0.2*G*sin(pi*x1)];

% 生成约束边界点
x2 = linspace(0, 1, sample_points/2)';
pf_constraint = [x2, min((1 - x2)/2, 2*(0.5 - x2))];

% 合并候选解
candidates = [pf_linear; pf_constraint];

% 应用约束条件筛选
valid = (candidates(:,1) + 2*candidates(:,2) >= 1) & ...
        (candidates(:,1) + 0.5*candidates(:,2) >= 0.5) & ...
        (candidates(:,1).^2 + candidates(:,2).^2 <= 1.96);

PF = candidates(valid, :);

% 非支配排序 (接口与NDSort1兼容)
[FrontNo, ~] = NDSort1(PF, 1);
PF = PF(FrontNo == 1, :);

% 边界修正 (原代码逻辑)
PF(PF(:,1)<0,1) = 0;
PF(PF(:,2)<0,2) = 0;
PF = unique(PF, 'rows');
end