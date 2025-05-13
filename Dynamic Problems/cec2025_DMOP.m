% ========================================================|%
% The 10 test functions are for cec2025 competition on    |
% continous dynamic multiobjective optimisation. This document is   |
% free to disseminate for academic use.                                  |
% --------------------------------------------------------|%
% The "time" term in the test suite is defined as:        |
%          t=1/nt*floor(tau/taut)                         |
% where - nt:    severity of change                       |
%       - taut:  frequency of change                      |
%       - tau:   current generation counter               |
% --------------------------------------------------------|%
% Any questions can be directed to                        |
%    Dr. Xiaozhong Yu at xzyu@smail.xtu.edu.cn.            |
%                                                         |
% ========================================================|%

function f=cec2025_DMOP(probID, x, tau, taut, nt)
% INPUT:
%       probID: test problem identifier (i.e. 'DF1')
%       x:      variable vector
%       tau:    current generation counter
%       taut:   frequency of change
%       nt:     severity of change
%
% OUTPUT:
%       f:      objective vector
%


% the first change occurs after T0 generations, that is, the
%generation at which a change occurs is (T0+1), (T0+taut+1), etc.
T0=50;

% calculate time instant
tau_tmp=max(tau+taut-(T0+1),0);
t=1/nt*floor(tau_tmp/taut);

n=length(x); % number of variables
switch (probID)
    %% Multimodality problems
    case 'DP1'
        G=sin(0.5*pi*t);
        a=0.2+2.8*abs(G);
        y=x(2:end)-G;
        g=1+sum((abs(G)*y.^2-10*cos(2*pi*y)+10));
        f(1)=g*(x(1)+0.1*sin(3*pi*x(1))).^a;
        f(2)=g*(1-x(1)+0.1*sin(3*pi*x(1))).^a;
    case 'DP2'
        p=floor(5*abs(sin(pi*t)));
        y=x(2:end)-cos(t);
        g=sum((4*y.^2-cos(2*p*pi*y)+1));
        f(1)=(1+g)*(1-x(1)+0.05.*sin(6*pi*x(1)));
        f(2)=(1+g)*(x(1)+0.05.*sin(6*pi*x(1)));
    case 'DP3'
        k = 10 * cos(2.5 * pi * t);
        a = 0.5 * abs(sin(pi * t));
        g = sum((x(2:end) - 0.5).^2 .* (1 + abs(cos(8 * pi * x(2:end)))));
        f(1) = (1 + g) * cos(0.5 * pi * x(1));
        % 平滑因子（sigmoid）
        s = 1 ./ (1 + exp(-100 * (a - x(1)))); 
        term1 = abs(k * (cos(0.5 * pi * x(1)) - cos(0.5 * pi * a)) + sin(0.5 * pi * a));
        term2 = sin(0.5 * pi * x(1));   
        f(2) = (1 + g) * (s * term1 + (1 - s) * term2);
    case 'DP4'
        G=sin(0.5*pi*t);
        k=2*floor(10*abs(G));
        y=x(2:end)-G;
        g=sum((4*y.^2-cos(k*pi*y)+1));
        f(1)=(1+g)*(x(1)+0.1*sin(3*pi*x(1)));
        f(2)=(1+g)*(1-x(1)+0.1*sin(3*pi*x(1)));
    case 'DP5'    
        G=sin(0.5*pi*t);
        a=0.2+2.8*abs(G);
        y=x(2:end)-G;
        g = sum(y.^2-10*cos(2*pi*y)+10);
        f(1)=(1+g)*(x(1)+0.1*sin(3*pi*x(1)))^a;
        f(2)=(1+g)*(1-x(1)+0.1*sin(3*pi*x(1)))^a;
        
    %% inregular changes
    case 'DP6'
        G=sin(0.5*pi*t);
        w=floor(10*G);
        g=1+sum((x(2:end)-G).^2);
        f(1)=g*(x(1)+0.02*sin(w*pi*x(1)));
        f(2)=g*(1-x(1)+0.02*sin(w*pi*x(1)));
    case 'DP7'
        G=sin(0.5*pi*t);
        a=2.25+2*cos(2*pi*t);
        tmp=G*sin(4*pi*x(1))/(1+abs(G));
        g=1+sum((x(2:end)-tmp).^2);
        f(1)=g*(x(1)+0.1*sin(3*pi*x(1)));
        f(2)=g*(1-x(1)+0.1*sin(3*pi*x(1))).^a;  
    case 'DP8'
        N=1+floor(10*abs(sin(0.5*pi*t)));
        g=1;
        for i=2:n
            tmp=x(i)-cos(4*t+x(1)+x(i-1));
            g=g+tmp.^2;
        end
        f(1)=g*(x(1)+max(0, (0.1+0.5/N)*sin(2*N*pi*x(1))));
        f(2)=g*(1-x(1)+max(0, (0.1+0.5/N)*sin(2*N*pi*x(1))));
    case 'DP9'
        G=sin(0.5*pi*t);
        W=10^(1+abs(G));
        g=sum((x(2:end)-G).^2);
        f(1)=(1+g)*(x(1)+0.05*sin(W*pi*x(1)));
        f(2)=(1+g)*(1-x(1)+0.05*sin(W*pi*x(1)));
    case 'DP10'
        G=abs(sin(0.5*pi*t));
 		p = max(floor(6 * G), 1);  % 确保 p 至少为1，避免退化

		% 防止 cot() 奇异点
		safe_term = atan(1 ./ tan(3 * pi * t^2 + eps));  % 加 eps 防止 tan 为 0

		g = sum((x(2:end) - 1/pi * abs(safe_term)).^2);
		f(1) = (1 + g) * (cos(0.5 * pi * x(1)))^2 + G;
		f(2) = (sin(0.5 * pi * x(1)))^2 + sin(0.5 * pi * x(1)) * (cos(p * pi * x(1)))^2 + G;

    otherwise
        disp('no such test problem.')
end
end
