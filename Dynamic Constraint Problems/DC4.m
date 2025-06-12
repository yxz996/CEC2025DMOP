classdef DC4 < PROBLEM 
% <multi> <real> <large/none> <dynamic> <constrained>
% taut --- 20 --- Number of generations for static optimization
% nt   --- 5 --- Number of distinct steps


    properties
        Optimums={};   
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            [obj.taut,obj.nt] = obj.ParameterSet(20,5);
            obj.M = 2;
            if isempty(obj.D); obj.D = 10; end
            obj.lower    = [0,-ones(1,obj.D-1)];
            obj.upper    = [1,ones(1,obj.D-1)];
            obj.encoding = ones(1,obj.D);
        end
        %% Evaluate solutions
        function Population = Evaluation(obj,varargin)
            PopDec     = obj.CalDec(varargin{1});
            PopObj     = obj.CalObj(PopDec);
            PopCon     = Constraint(obj,PopObj);
            % Attach the current number of function evaluations to solutions
            Population = SOLUTION(PopDec,PopObj,PopCon,zeros(size(PopDec,1),1)+obj.FE);
            obj.FE     = obj.FE + length(Population);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            t = floor(obj.FE/obj.N/obj.taut)/obj.nt; 
            G = sin(0.5*pi*t); 
            g = 1 + sum((PopDec(:,2:end)-G).^2,2);
            PopObj(:,1) =g.*PopDec(:,1);
            PopObj(:,2) =g.*(1-PopDec(:,1));
        end

        %% Generate points on the Pareto front  
        function R = GetOptimum(obj,N) 
            tt = floor((0:obj.maxFE)/obj.N/obj.taut)/obj.nt;
            H = sin(0.01.*pi.*tt);
            H = unique(round(H*1e6)/1e6);

            V(:,1) = (0:1/N:1)';
            V(:,2) = 1 - V(:,1) ;

            for i = 1 : length(H)
                pf=[];  pf2=[];  
                t  = (i-1)/obj.nt;
                W = 0.5*pi*t;
                pf = V;
                c11 = (sin(W)-cos(W))/(sin(W)+cos(W)).* pf(:,1) - pf(:,2) + (-2.2*(1-cos(W))+1.3)./(cos(W)+sin(W));
                c12 = (sin(W)-cos(W))/(sin(W)+cos(W)).* pf(:,1) - pf(:,2) + (-2.2*(1-cos(W))+1.8)./(cos(W)+sin(W));
                c1a = c11 .* c12;
                c31 = (sin(W)-cos(W))/(sin(W)+cos(W)).* pf(:,1) - pf(:,2) + (-2.2*(1-cos(W))+2.6)./(cos(W)+sin(W));
                c32 = (sin(W)-cos(W))/(sin(W)+cos(W)).* pf(:,1) - pf(:,2) + (-2.2*(1-cos(W))+3.1)./(cos(W)+sin(W));
                c3a = c31 .* c32;
                pf(c1a<0 | c3a<0,:) = [];
                if 0.5*t-floor(0.5*t) >= 0.25 && 0.5*t-floor(0.5*t) <= 0.75
                    pf2(:,1) = [zeros(N+1,1); 1+V(:,1)];
                    pf2(:,2) = [1+V(:,1); zeros(N+1,1)];
                    c11 = (sin(W)-cos(W))/(sin(W)+cos(W)).* pf2(:,1) - pf2(:,2) + (-2.2*(1-cos(W))+1.3)./(cos(W)+sin(W));
                    c12 = (sin(W)-cos(W))/(sin(W)+cos(W)).* pf2(:,1) - pf2(:,2) + (-2.2*(1-cos(W))+1.8)./(cos(W)+sin(W));
                    c1a = c11 .* c12;
                    c31 = (sin(W)-cos(W))/(sin(W)+cos(W)).* pf2(:,1) - pf2(:,2) + (-2.2*(1-cos(W))+2.6)./(cos(W)+sin(W));
                    c32 = (sin(W)-cos(W))/(sin(W)+cos(W)).* pf2(:,1) - pf2(:,2) + (-2.2*(1-cos(W))+3.1)./(cos(W)+sin(W));
                    c3a = c31 .* c32;
                    pf2(c1a<0 | c3a<0,:) = [];
                elseif 0.5*t-floor(0.5*t) < 0.25
                    f12 = -(2.2*(1-cos(W))-1.3)/(cos(W)+sin(W));
                    f11 = (2.2*(1-cos(W))-1.3)/(sin(W)-cos(W));
                    f22 = -(2.2*(1-cos(W))-1.8)/(cos(W)+sin(W));
                    f21 = (2.2*(1-cos(W))-1.8)/(sin(W)-cos(W));
                    pf2 = [V.*repmat([f11 f12],size(V,1),1);V.*repmat([f21 f22],size(V,1),1)];
                    c11 = pf2(:,1) + pf2(:,2) -1;
                    pf2(c11<0,:) = [];
                elseif 0.5*t-floor(0.5*t) > 0.75   
                    f12 = -(2.2*(1-cos(W))-3.1)/(cos(W)+sin(W));
                    f11 = (2.2*(1-cos(W))-3.1)/(sin(W)-cos(W));
                    f22 = -(2.2*(1-cos(W))-2.6)/(cos(W)+sin(W));
                    f21 = (2.2*(1-cos(W))-2.6)/(sin(W)-cos(W));
                    pf2 = [V.*repmat([f11 f12],size(V,1),1);V.*repmat([f21 f22],size(V,1),1)];
                    c11 = pf2(:,1) + pf2(:,2) -1;
                    pf2(c11<0,:) = [];
                end
                pf = [pf; pf2];
                pf(NDSort1(pf,1)~=1,:) = [];

                obj.Optimums(i,:) = {H(i),pf};
            end
            R=cat(1,obj.Optimums{:,2});
        end

        %% Generate the image of Pareto front  
        function R = GetPF(obj)

            R = obj.GetOptimum(100);
        end



        %% Display a population in the objective space 
        function DrawObj(obj,Population) 
                t      = floor(Population.adds/obj.N/obj.taut)/obj.nt;
                H      = sin(0.01.*pi.*t);
                H      = round(H*1e6)/1e6;
                change = [0;find(H(1:end-1)~=H(2:end));length(H)];
                allH   = cell2mat(obj.Optimums(:,1));
                tempStream = RandStream('mlfg6331_64','Seed',2);
                for i = 1 : length(change)-1
                    color = rand(tempStream,1,3);
                    showdata=Population(change(i)+1:change(i+1)).objs+(i-1)*0.5;
                    Draw(showdata,'o','MarkerSize',5,'Marker','o','Markerfacecolor',sqrt(color),'Markeredgecolor',color,{'\it f\rm_1','\it f\rm_2',[]});
                    ax=Draw([],'-','LineWidth',1,'Color',color);
                    drawPF=obj.Optimums{find(H(change(i)+1)==allH,1),2}+(i-1)*0.5;
                    for i2=1:length(drawPF(:,1))-1
                        fx2=drawPF(i2+1,1);
                        fx1=drawPF(i2,1);
                        dis=fx2-fx1;
                        if dis>0&&dis<0.09
                            plot(ax,[drawPF(i2,1),drawPF(i2+1,1)] , [drawPF(i2,2),drawPF(i2+1,2)],'-','Color',[.0 .0 .0],'LineWidth',1);
                        end
                    end
                end
        end

        %% Calculate the metric value  
        function score = CalMetric(obj,metName,Population)
                t      = floor(Population.adds/obj.N/obj.taut)/obj.nt;
                H      = sin(0.01.*pi.*t);
                H      = round(H*1e6)/1e6;
                change = [0;find(H(1:end-1)~=H(2:end));length(H)];
                Scores = zeros(1,length(change)-1);
                allH   = cell2mat(obj.Optimums(:,1));
                for i = 1 : length(change)-1
                    subPop    = Population(change(i)+1:change(i+1));
                    Scores(i) = feval(metName,subPop,obj.Optimums{find(H(change(i)+1)==allH,1),2});
                end
                score = mean(Scores);
        end

    end
end

%% Calculate constraint violations values 
function PopCon = Constraint(obj,PopObj)
    t = floor(obj.FE/obj.N/obj.taut)/obj.nt;
    G=sin(0.5*pi*t);
     W = 0.5*pi*t;                
     c11 = (sin(W)-cos(W))/(sin(W)+cos(W)).* PopObj(:,1) - PopObj(:,2) + (-2.2*(1-cos(W))+1.3)./(cos(W)+sin(W));
     c12 = (sin(W)-cos(W))/(sin(W)+cos(W)).* PopObj(:,1) - PopObj(:,2) + (-2.2*(1-cos(W))+1.8)./(cos(W)+sin(W));
     c1 = c11 .* c12;
     c31 = (sin(W)-cos(W))/(sin(W)+cos(W)).* PopObj(:,1) - PopObj(:,2) + (-2.2*(1-cos(W))+2.6)./(cos(W)+sin(W));
     c32 = (sin(W)-cos(W))/(sin(W)+cos(W)).* PopObj(:,1) - PopObj(:,2) + (-2.2*(1-cos(W))+3.1)./(cos(W)+sin(W));
     c3 = c31 .* c32;

    PopCon(:,1)=-c1;
    PopCon(:,2)=-c3;
end