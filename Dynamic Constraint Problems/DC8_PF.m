[file name]: DC8_DF_PF.m
[file content begin]
function PF = DC8_DF_PF(tau, taut, nt, N, sample_points)
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
G = 0.5 * abs(sin(0.5 * pi * t));

% 生成基础PF点集（直线部分）
x1 = linspace(0, 1, sample_points)';
pf_linear = [x1, 1 - x1];

% 生成单位圆上的点（球面部分）
theta = linspace(0, pi/2, sample_points)';
pf_sphere = (0.6+G) * [cos(theta), sin(theta)];

% 合并候选解
candidates = [pf_linear; pf_sphere];

% 应用约束条件筛选
valid = ((0.2+G)*candidates(:,1).^2 + candidates(:,2) >= 2) & ...
        (0.7*candidates(:,1).^2 + candidates(:,2) >= 2.5) & ...
        (candidates(:,1).^2 + candidates(:,2).^2 >= (0.6+G)^2);

PF = candidates(valid, :);

% 非支配排序 (接口与NDSort1兼容)
[FrontNo, ~] = NDSort1(PF, 1);
PF = PF(FrontNo == 1, :);

% 边界修正 (原代码逻辑)
if ~isempty(PF)
    [~, idx] = min(PF(:,1)); PF(idx(1),1) = 0;
    [~, idx] = min(PF(:,2)); PF(idx(1),2) = 0;
end
end