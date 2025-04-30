[file name]: DC3_DF_PF.m
[file content begin]
function PF = DC3_DF_PF(tau, taut, nt, N, sample_points)
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
H = 1;
if G < 0
    H = -1;
end

% 生成单位圆参考点
theta = linspace(0, pi/2, sample_points)';
V = [cos(theta), sin(theta)];
V = V ./ vecnorm(V, 2, 2); % 归一化处理

% 应用动态约束
ll = cos(5*(atan((V(:,2)./V(:,1)).^H)).^4).^6;
denom1 = 0.9 + (0.1 + 0.7*abs(G))*ll;
denom2 = 0.9 + (0.8 - 0.7*abs(G))*ll;
valid = 1.1 - (V(:,1)./denom1).^2 - (V(:,2)./denom2).^2 >= 0;

% 构建候选解并筛选
candidates = V(valid, :);
PF = candidates;

% 非支配排序 (接口与NDSort1兼容)
[FrontNo, ~] = NDSort1(PF, 1);
PF = PF(FrontNo == 1, :);

% 后处理：移除异常点
invalid = any(imag(PF) ~= 0, 2) | any(PF < 0, 2);
PF(invalid, :) =