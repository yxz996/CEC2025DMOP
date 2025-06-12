classdef DC7 < PROBLEM %
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
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
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
            G = 2 * floor(10*abs(mod(t+1,2)-1)+1e-4);
            g = 1+sum( ((PopDec(:,2:end)-0.5).^2 - cos(pi*(PopDec(:,2:end)-0.5)) + 1).^2 ,2);
            PopObj(:,1) =g.*(PopDec(:,1)+0.25*sin(pi.*PopDec(:,1)));
            PopObj(:,2) =g.*(1-PopDec(:,1)+0.25*sin(pi*PopDec(:,1)));
%             PopObj(PopObj < 1e-18) = 0;
        end

        %% Generate points on the Pareto front  
        function R = GetOptimum(obj,N) 
            tt = floor((0:obj.maxFE)/obj.N/obj.taut)/obj.nt;
            H = sin(0.01.*pi.*tt);
            H = unique(round(H*1e6)/1e6);

            x1 = (0:1/(500-1):1)';
            pf1(:,1) = x1 +  0.25*sin(pi*x1);
            pf1(:,2) = 1- x1 +  0.25*sin(pi*x1);

            for i = 1 : length(H)
                pf=pf1;
                t  = (i-1) / obj.nt;
                G = 2 * floor(10*abs(mod(t+1,2)-1)+1e-4);
                c11 = (pf(:,1)).^2 + (pf(:,2)).^2 - (1.5 + 0.4*sin(4*atan(pf(:,2)./pf(:,1))).^16).^2;
                c12 = (pf(:,1)).^2 + (pf(:,2)).^2 - (1.3 - 0.45*sin(G*atan(pf(:,2)./pf(:,1))).^2).^2;
                c1 = c11.*c12<0;
                pf(c1,:) = [];

    
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
                    Draw(drawPF,'*','LineWidth',1,'Color',[0 0 0]);
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
    G = 2 * floor(10*abs(mod(t+1,2)-1)+1e-4);
    c11 = (PopObj(:,1)).^2 + (PopObj(:,2)).^2 - (1.5 + 0.4*sin(4*atan(PopObj(:,2)./PopObj(:,1))).^16).^2;
    c12 = (PopObj(:,1)).^2 + (PopObj(:,2)).^2 - (1.3 - 0.45*sin(G*atan(PopObj(:,2)./PopObj(:,1))).^2).^2;
    PopCon = -c11.*c12;
end