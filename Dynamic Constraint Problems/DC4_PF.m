[file name]: DC4_DF_PF.m
[file content begin]
function PF = DC4_DF_PF(tau, taut, nt, N, sample_points)
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
W = 0.5 * pi * t;
A = (sin(W)-cos(W))/(sin(W)+cos(W));
B = -2.2*(1-cos(W));

% 生成基础PF点集
V = linspace(0, 1, sample_points)';
base_PF = [V, 1 - V];

% 计算约束边界参数
phase = 0.5*t - floor(0.5*t);
if phase >= 0.25 && phase <= 0.75
    % 中间相位生成扩展点
    pf_ext = [zeros(sample_points,1)+1e-6, linspace(1,0,sample_points)';
              linspace(1,0,sample_points)', zeros(sample_points,1)+1e-6];
elseif phase < 0.25
    % 低相位生成倾斜边界
    f12 = -(B - 1.3)/(sin(W)+cos(W));
    f11 = (B - 1.3)/(A);
    pf_ext = [V*f11, V*f12];
else
    % 高相位生成另一倾斜边界
    f12 = -(B - 3.1)/(sin(W)+cos(W));
    f11 = (B - 3.1)/(A);
    pf_ext = [V*f11, V*f12];
end

% 合并候选解
candidates = [base_PF; pf_ext];

% 应用约束条件筛选
valid = (A*candidates(:,1) - candidates(:,2) + (B + 1.3)/(sin(W)+cos(W))) .* ...
       