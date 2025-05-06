% ========================================================|%
% PF calculation for 10 cec2015 test functions on         |
% dynamic multiobjective optimisation. This document is   |
% free to disseminate for academic use.                   |
% --------------------------------------------------------|%
% The "time" term in the test suite is defined as:        |
%          t=1/nt*floor(tau/taut)                         |
% where - nt:    severity of change                       |
%       - taut:  frequency of change                      |
%       - tau:   current generation counter               |
% --------------------------------------------------------|%
% Any questions can be directed to                        |
%     Dr. Xiaozhong Yu at xzyu@smail.xtu.edu.cn.               |
%                                                         |
% ========================================================|%

function h=cec2025_DMOP_PF(probID,tau, taut, nt)
% INPUT:
%       probID: test problem identifier (i.e. 'DF1')
%       tau:    current generation counter
%       taut:   frequency of change
%       nt:     severity of change
%
% OUTPUT:
%       h:      nondominated solutions
%


% the first change occurs after T0 generations, that is, the
%generation at which a change occurs is (T0+1), (T0+taut+1), etc.
T0=50;

% calculate time instant
tau_tmp=max(tau+taut-(T0+1),0);
t=1/nt*floor(tau_tmp/taut);

g=1;
H=50; % number of divisions along each objective.
switch (probID)
    %% Multimodality problems
    case 'DP1'
        x=linspace(0,1,1500);
        G=sin(0.5*pi*t);
        a=0.2+2.8*abs(G);
        f1=g*(x+0.1*sin(3*pi*x)).^a;
        f2=g*(1-x+0.1*sin(3*pi*x)).^a;
        [h]=get_PF({f1,f2}, false);
    case 'DP2'
        x=linspace(0,1,1500);       
        g=0;
        f1=(1+g)*(1-x+0.05.*sin(6*pi*x));
        f2=(1+g)*(x+0.05.*sin(6*pi*x));
        [h]=get_PF({f1,f2}, false);
    case 'DP3'
        x=linspace(0,1,1500);
        k=10*cos(2.5*pi*t);
        a=0.5*abs(sin(pi*t));
        g=0;
        f1=(1+g)*cos(0.5*pi*(x));
        if x<=a
            f2=(1+g)*abs(k*cos(0.5*pi*x)-cos(0.5*pi*a))+sin*(0.5*pi*a);
        else
            f2=(1+g)*sin(0.5*pi*x);
        end
        [h]=get_PF({f1,f2}, false);
    case 'DP4'
        x=linspace(0,1,1500);
        g=0;
        f1=(1+g)*(x+0.1*sin(3*pi*x));
        f2=(1+g)*(1-x+0.1*sin(3*pi*x));
        [h]=get_PF({f1,f2}, false);
    case 'DP5'
        x=linspace(0,1,1500); 
        G=sin(0.5*pi*t);
        a=0.2+2.8*abs(G);
        g=0;
        f1=(1+g)*(x+0.1*sin(3*pi*x))^a;
        f2=(1+g)*(1-x+0.1*sin(3*pi*x))^a;
        [h]=get_PF({f1,f2}, false);
    %% inregular changes        
    case 'DP6'
        x=linspace(0,1,1500);
        G=sin(0.5*pi*t);
        w=floor(10*G);
        f1=g*(x+0.02*sin(w*pi*x));
        f2=g*(1-x+0.02*sin(w*pi*x));
        [h]=get_PF({f1,f2}, false);
    case 'DP7'
        x=linspace(0,1,1500);
        a=2.25+2*cos(2*pi*t);
        f1=g*(x+0.1*sin(3*pi*x));
        f2=g*(1-x+0.1*sin(3*pi*x)).^a;
        [h]=get_PF({f1,f2}, false);
    case 'DP8'
        x=linspace(0,1,1500);
        N=1+floor(10*abs(sin(0.5*pi*t)));
        f1=g*(x+max(0, (0.1+0.5/N)*sin(2*N*pi*x)));
        f2=g*(1-x+max(0, (0.1+0.5/N)*sin(2*N*pi*x)));
        [h]=get_PF({f1,f2}, true); 
    case 'DP9'
        x=linspace(0,1,1500);
        G=sin(0.5*pi*t);
        W=10^(1+abs(G));
        g=0;
        f1=(1+g)*(x+0.05*sin(W*pi*x));
        f2=(1+g)*(1-x+0.05*sin(W*pi*x));
        [h]=get_PF({f1,f2}, true);
    case 'DP10'
        x=linspace(0,1,1500);
        G=abs(sin(0.5*pi*t));
        p=floor(6*G);
        g=0;
        f1=(1+g)*(cos(0.5*pi*x))^2 + G;
        f2=(sin(0.5*pi*x)).^2 +sin(0.5*pi*x).*(cos(p*pi*x)).^2+ G;
        [h]=get_PF({f1,f2}, true);
    otherwise
        disp('no such test problem.')
end
end

% helper function: identify nondominated solutions
% f={f1,f2,f3,...}
% nondominate: nondominated sorting (true), otherwise (false)
function [h]=get_PF(f, nondominate)
    ncell=length(f);
    s=numel(f{1});
    h=[];
    for i=1:ncell
        fi=reshape(f{i},s,1);
        h=[h,fi];
    end
    
    if nondominate
        in=get_skyline(h);
        h=h(in,:);
    end
end


% helper function: find the indices of nondominated solutions
function [Im]=get_skyline(x)

    n=size(x,1);
    
    if n==1
        Im=1;
    else
        Iold=1:n;
        Im=[]; %Im holds the index of the solutions in the first front;
        
        % sieving method
        %[x1_index,x1]=min(sum(x'));
        z=x;
        n1=n-1;
        while (n1>0)
            [V]=get_vector(z);
            I=[];
            J=[];
            %J=find(V==1);
            %I=find(x<0,1);
            for k=1:n1
                if V(k)==1
                    J=[J,k];
                elseif V(k)==-1
                    I=[I,k];
                end
            end
            
            if isempty(I)
                Im=[Im,Iold(1)];
               % plot(x(Iold(1),1),x(Iold(1),2),'ro');
               % hold on;
            end
            In=setdiff(2:n1+1,J+1);
            
            In=Iold(In);
            if length(In)==1
                Im=[Im,In];
                break;
            end
            
            z=x(In,:);
            n1=length(In)-1; %because the size of get_vector(z)is smaller than In.
            Iold=In;
        end
    end
end

function [myvector,ob_count]=get_vector(x)
    ob_count=0;
    [N,D]=size(x);
    
    myvector=[];
    for j=2:N
        dom_less=0;
        dom_equal=0;
        dom_more=0;
        %p=sign( x(1,:)-x(2,:))
        for k = 1 : D
            if (x(1,k) < x(j,k))
                dom_less = dom_less + 1;
            elseif (x(1, k) == x(j,k))
                dom_equal = dom_equal + 1;
            else
                dom_more = dom_more + 1;
            end
            ob_count=ob_count+1;
            if dom_less>1&&dom_more>1
                break;
            end
        end
        if k<D
            myvector(j-1)=0;
            continue;
        end
        if dom_less == 0 && dom_equal ~= D
            myvector(j-1)=-1;
        elseif dom_more == 0 && dom_equal ~= D
            myvector(j-1)=1;
        else
            myvector(j-1)=0;
        end
    end
end
