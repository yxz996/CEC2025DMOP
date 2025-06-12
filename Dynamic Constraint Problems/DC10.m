classdef DC10< PROBLEM %
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
            G = sin(0.5*pi*t);
            y = (PopDec(:,2:end)-G);
            g = 1 + sum((abs(G)*y.^2 - cos(pi*y) +1).^2,2);
            PopObj(:,1) =g.*PopDec(:,1);
            PopObj(:,2) =g.*(1-PopDec(:,1));
        end

        %% Generate points on the Pareto front  
        function R = GetOptimum(obj,N) 
            tt = floor((0:obj.maxFE)/obj.N/obj.taut)/obj.nt;
            H = sin(0.01.*pi.*tt);
            H = unique(round(H*1e6)/1e6);

            x1 = (0:1/(501-1):1)';
            X = UniformPoint(500,2);
            pf1 = 1.2* X./repmat((sum(X.^1.5,2)).^(2/3),1,2);

            for i = 1 : length(H)
                pf=[];
                t = (i-1) / obj.nt;
                G = sin(0.5*pi*t);
                pf(:,1) = x1 ;
                pf(:,2) = 1-x1;
                c1 = pf(:,1).^(0.5) + pf(:,2).^(0.5) - (0.95+0.5*abs(G) );
                pf(c1>0,:) = [];
                pf = [pf;pf1];
                pf(NDSort(pf,1)~=1,:) = [];

    
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
                    Draw(drawPF,'.','LineWidth',1,'Color',[0 0 0]);
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
    G = sin(0.5*pi*t);
    c11 = (PopObj(:,1)).^1.5+ (PopObj(:,2)).^1.5 - 1.2^1.5;
    c12 = PopObj(:,1).^(0.5) + PopObj(:,2).^(0.5) - (0.95+0.5*abs(G) );
    c21 =  0.8*PopObj(:,1) + PopObj(:,2) - (2.5 + 0.08*sin(2*pi*(PopObj(:,2)-PopObj(:,1))));
    c22 = (0.93 + abs(G)/3)*PopObj(:,1) + PopObj(:,2) - (2.7+abs(G)/2 + 0.08*sin(2*pi*(PopObj(:,2)-PopObj(:,1))));
    PopCon(:,1)  = -c11.*c12;
    PopCon(:,2)  = -c21.*c22;
end