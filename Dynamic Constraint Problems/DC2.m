classdef DC2 < PROBLEM 
% <multi> <real> <large/none> <dynamic> <constrained>
% taut --- 20 --- Number of generations for static optimization
% nt   --- 5 --- Number of distinct steps

    properties
        Optimums={};   
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
            G = sin(0.5*pi*t);
            g = 1 + sum((PopDec(:,2:end)- G).^2 + sin((PopDec(:,2:end)- G)*0.5*pi).^2,2);
            PopObj(:,1) =g.*(PopDec(:,1) + 0.2*G*sin(pi*PopDec(:,1)));
            PopObj(:,2) =g.*(1-PopDec(:,1) + 0.2*G*sin(pi*PopDec(:,1)));
            PopObj(PopObj < 1e-18) = 0;
        end

        %% Generate points on the Pareto front 
        function R = GetOptimum(obj,N) 
            tt = floor((0:obj.maxFE)/obj.N/obj.taut)/obj.nt;
            H = sin(0.01.*pi.*tt);
            H = unique(round(H*1e6)/1e6);

            x1 = linspace(0,1,N)';
            pf2=[]; pf2(:,1)=x1;
            pf2(:,2) = [1-2*x1(x1<1/3); 0.5-0.5*x1(x1>=1/3)];

            for i = 1 : length(H)
                pf=[];  
                t  = (i-1)/obj.nt;
                G = sin(0.5*pi*t);
                pf(:,1) = x1 + 0.2*G*sin(pi*x1);
                pf(:,2) = 1-x1 + 0.2*G*sin(pi*x1);
                pf(pf==min(pf)) = 0;  pf(pf==max(pf)) = 1;
                c11 = pf(:,1) + 2*pf(:,2) - 1;
                c12 = pf(:,1) + 0.5*pf(:,2) - 0.5;
                c1 = c11 .* c12;
                pf(c1<0,:) = [];
                LengthPF = size(pf,1);
                pf = [pf; pf2];
                Select = NDSort(-pf,1) == inf;
                Select(1:LengthPF) = false;
                pf(Select,:) = [];
                pf(NDSort1(pf,1)~=1,:) = [];

                pf=unique(pf,'rows');

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
                    Draw(drawPF,'-','LineWidth',1,'Color',[0 0 0]);
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
    Con1=PopObj(:,1) + 2*PopObj(:,2) - 1;
    Con2=PopObj(:,1) + 0.5*PopObj(:,2) - 0.5;
    PopCon(:,1)=-Con1.*Con2;
    PopCon(:,2)=PopObj(:,1).^2 + PopObj(:,2).^2 - 1.4.^2;
end