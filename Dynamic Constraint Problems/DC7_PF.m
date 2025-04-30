[file name]: DC7_DF_PF.m
[file content begin]
function PF = DC7_DF_PF(tau, taut, nt, N, sample_points)
% INPUT:
%       tau:        current generation counter
%       taut:       change frequency
%       nt:         change severity
%       N:          population size (for t calculation)
%       sample_points: PF采样点数
%
% OUTPUT:
%       PF:         Pareto front solutions

% 时间参数计算
t = floor(tau/N/taut)/nt;
G = 2 * floor(10*abs(mod(t+1,2)-1 + 1e-4);

% 生成参考点
x1 = linspace(0, 1, sample_points)';
P1 = x1 + 0.25*sin(pi*x1);
P2 = 1 - x1 + 0.25*sin(pi*x1);

% 应用约束条件
theta = atan2(P2, P1); % 使用atan2避免除零错误
valid = ((P1.^2 + P2.^2 < (1.5 + 0.4*sin(4*theta).^16).^2) & ...
        ((P1.^2 + P2.^2 > (1.3 - 0.45*sin(G*theta).^2).^2);

PF = [P1(valid), P2(valid)];

% 非支配排序 (接口与NDSort1兼容)
[FrontNo, ~] = NDSort1(PF, 1);
PF = PF(FrontNo == 1, :);

% 边界修正 (原代码逻辑)
if ~isempty(PF)
    [~, idx] = min(PF(:,1)); PF(idx(1),1) = 0;
    [~, idx] = min(PF(:,2)); PF(idx(1),2) = 0;
end
end