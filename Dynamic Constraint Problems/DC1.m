classdef DC1 < PROBLEM 
% <multi> <real> <large/none> <dynamic> <constrained>
% taut --- 20 --- Number of generations for static optimization
% nt   --- 5 --- Number of distinct steps

    properties
        Optimums={};   % Point sets on all Pareto fronts
        output1=NaN(1,200);
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
           
            Population = SOLUTION(PopDec,PopObj,PopCon,zeros(size(PopDec,1),1)+obj.FE);
            obj.FE     = obj.FE + length(Population);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            t = floor(obj.FE/obj.N/obj.taut)/obj.nt; 
            g = 1+sum(PopDec(:,2:end).^2,2);
            PopObj(:,1) =g.*PopDec(:,1);
            PopObj(:,2) =g.*sqrt(1-PopDec(:,1).^2);
        end

        %% Generate points on the Pareto front  
        function R = GetOptimum(obj,N) 
            tt = floor((0:obj.maxFE)/obj.N/obj.taut)/obj.nt;
            H = sin(0.01.*pi.*tt);
            H = unique(round(H*1e6)/1e6);

            x = linspace(0,1,N)';
            V(:,1) = x;
            V(:,2) = 1-x;
            V=V./repmat(sqrt(sum(V.^2,2)),1,2);

            for i = 1 : length(H)
                pf=[]; pf2=[]; 
                pf = V;
                t  = (i-1)/obj.nt;
                G = cos(pi*t);
                c11 = 3 - G - exp(pf(:,1)) - 0.3*sin(4*pi*pf(:,1)) - pf(:,2);
                c12 = 4.1 - (1 + 0.3*pf(:,1).^2 + pf(:,1)) - 0.3*sin(4*pi*pf(:,1)) - pf(:,2);
                Con = c11.*c12;
                pf(Con>0,:) = [];
                syms x y
                s=solve(3 - G - exp(x) - 0.3*sin(4*pi*x) - y==0,y==0,x,y);
                X=min(double(s.x));
                pf2(:,1) = linspace(0,X,N)';
                pf2(:,2) = 3 - G - exp(pf2(:,1)) - 0.3*sin(4*pi*pf2(:,1));
                OC = pf2(:,1).^2 + pf2(:,2).^2 -1;%目标函数
                pf2(OC<0,:) = [];
                pf = [pf;pf2];
                pf(NDSort1(pf,1)~=1,:) = [];

                obj.Optimums(i,:) = {H(i),pf};
            end
            R=cat(1,obj.Optimums{:,2});
        end

        %% Generate the image of Pareto front  
        function R = GetPF(obj)
            R = obj.GetOptimum(100);
        end


        %% Display a population in the objective space 显示种群和truePF
        function DrawObj(obj,Population) 
                t      = floor(Population.adds/obj.N/obj.taut)/obj.nt;
                H      = sin(0.01.*pi.*t);
                H      = round(H*1e6)/1e6;
                change = [0;find(H(1:end-1)~=H(2:end));length(H)];
                allH   = cell2mat(obj.Optimums(:,1));
                tempStream = RandStream('mlfg6331_64','Seed',2);
%                 flag=1;
                for i = 1 : length(change)-1
                    color = rand(tempStream,1,3);

                    Popp=Population(change(i)+1:change(i+1));
                    infeasible=any(Popp.cons>0,2);
                    showdata1=Popp(~infeasible).objs+(i-1)*0.5;
                    showdata2=Popp(infeasible).objs+(i-1)*0.5;
                    Draw(showdata1,'o','MarkerSize',5,'Marker','o','Markerfacecolor',sqrt(color),'Markeredgecolor',color,{'\it f\rm_1','\it f\rm_2',[]});
                    color=[.0 .0 .0];
                    Draw(showdata2,'o','MarkerSize',5,'Marker','o','Markerfacecolor',sqrt(color),'Markeredgecolor',color,{'\it f\rm_1','\it f\rm_2',[]});

                    ax=Draw([],'-','LineWidth',1,'Color',color);
                    drawPF=obj.Optimums{find(H(change(i)+1)==allH,1),2}+(i-1)*0.5;
                    for i2=1:length(drawPF(:,1))-1
                        fx2=drawPF(i2+1,1);
                        fx1=drawPF(i2,1);
                        dis=fx2-fx1;
                        if dis>0&&dis<0.03
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
                tt = int32(unique(Population.adds/100));
                tt=tt(end);
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
    G=cos(pi*t);
    Con1=3 - G - exp(PopObj(:,1)) - 0.3*sin(4*pi*PopObj(:,1)) - PopObj(:,2);
    Con2=4.1 - (1 + 0.3*PopObj(:,1).^2 + PopObj(:,1)) - 0.3*sin(4*pi*PopObj(:,1)) - PopObj(:,2);
    PopCon(:,1)=Con1.*Con2;
end