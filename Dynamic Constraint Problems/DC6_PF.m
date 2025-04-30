[file name]: DC6_DF_PF.m
[file content begin]
function PF = dcp2_DF_PF(tau, taut, nt, N, sample_points)
% INPUT:
%       tau:        current generation counter
%       taut:       change frequency  
%       nt:         change severity
%       N:          population size (for t calculation)
%       sample_points: number of PF points to generate
%
% OUTPUT:
%       PF:         Pareto front solutions

% Time parameter calculation
t = floor(tau/N/taut)/nt;
G = sin(0.5*pi*t);

% Generate reference points
x1 = linspace(0, 1, sample_points)';
P1 = x1 + 0.25*G*sin(pi*x1);
P2 = 1 - x1 + 0.25*G*sin(pi*x1);

% Apply constraints
valid = (4*P1 + P2 -1).*(0.3*P1 + P2 -0.3) >= 0 & ...       % Reverse c1 condition
        (1.85 - P1 - P2 - (0.3*sin(3*pi*(P2-P1))).^2).*... % Original c2 condition
        (P1 + P2 -1.3) <= 0;

% Construct candidate solutions
PF = [P1(valid), P2(valid)];

% Add special case solutions (from DC6's GetOptimum)
x2 = 0:0.001:2;
pf2 = [repmat(x2', 2001,1), kron(0:0.001:2, ones(2001,1))];
cc = (1.85 - pf2(:,2) - pf2(:,1) - (0.3*sin(3*pi*(pf2(:,2)-pf2(:,1)))).^2);
pf2(abs(cc) > 2e-3,:) = [];
PF = [PF; pf2];

% Non-dominated sorting with NDSort1 interface
[FrontNo, ~] = NDSort1(PF, 1);  % 1表示只处理第一个约束
PF = PF(FrontNo == 1, :);

% Boundary adjustment (from original code)
[~,idx] = min(PF(:,1)); PF(idx,1) = 0;
[~,idx] = min(PF(:,2)); PF(idx,2) = 0;
end