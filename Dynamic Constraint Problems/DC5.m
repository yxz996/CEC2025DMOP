classdef DC5<PROBLEM
% <multi> <real> <large/none> <constrained> <dynamic>
% taut --- 10 --- Number of generations for static optimization
% nt   --- 10 --- Number of distinct steps
%preEvolution   --- 30 --- Number of pre Evolution

    properties
        taut = 10; % Number of function evaluations for each change
        nt = 10;   % Time scale of dynamic change
        preEvolution = 30;  % The maximum generation of first environment
        Optimums={};            
        FeasibleRegions = {}; 
        t;
        max_time;
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            % Initialize dynamic parameters and problem settings
            [obj.taut, obj.nt, obj.preEvolution] = obj.ParameterSet(10, 10, 30);
            if isempty(obj.M); obj.M = 2; end
            if isempty(obj.D); obj.D = 10; end
            obj.lower = zeros(1,obj.D);
            obj.upper = ones(1,obj.D);
            obj.encoding = ones(1, obj.D);
            obj.taut = obj.taut;
            obj.nt = obj.nt;
            obj.preEvolution = obj.preEvolution*obj.N;
            obj.max_time=floor((obj.maxFE-obj.N)/obj.N/obj.taut)/obj.nt;
            obj.maxFE=obj.maxFE+obj.preEvolution;
            obj.t = 0;
            obj.PreGenerateOptimumsAndFeasibleRegions();
        end

        %% Pre-generate POFs and feasible regions for all time points
        function PreGenerateOptimumsAndFeasibleRegions(obj)
            % Generate POF (considering constraints) and feasible regions for all time points
            t=min(floor((0:obj.maxFE-obj.preEvolution)/obj.N/obj.taut)/obj.nt,obj.max_time);
            t = unique(t);
            N = 10000; % Number of points for POF
            x1 = (0:1/N:1)';
            addtion1 = (0.001:0.001:0.4)';
            for i = 1:length(t)
                ti = t(i);
                pf=[];  
                G = abs(sin(0.5*pi*ti));
                P1=x1+G;
                P2=1-x1+G;
                addtion2=repmat(min(P1),length(addtion1),1);
                addtion3=addtion1+max(P1);
                P1 = [P1;addtion2;addtion3];
                P2 = [P2;addtion3;addtion2];
                c1 = cos(-0.15*pi)*P2 - sin(-0.15*pi)*P1 - (2*sin(5*pi*(sin(-0.15*pi)*P2 + cos(-0.15*pi)*P1))).^6;
                Con =-c1;
                pf(:,1) = P1;
                pf(:,2) = P2;
                pf(Con>0,:) = [];
                pf(NDSort(pf,1)~=1,:) = [];
                pf=sortrows(pf);
                obj.Optimums(i,:)={ti,pf};

                % Generate feasible region
                [x, y] = meshgrid(linspace(0, 2.5, 1000), linspace(0, 2.5, 1000));
                z = nan(size(x));

                c1 = cos(-0.15*pi)*y - sin(-0.15*pi)*x - (2*sin(5*pi*(sin(-0.15*pi)*y + cos(-0.15*pi)*x))).^6;
                Con =-c1;

                feasible = Con<=0 & x+y>=1+2*G; 
                z(feasible) = 0;
                feasible_region = {x, y, z};
                obj.FeasibleRegions(i, :) = {ti, feasible_region};
            end
        end

        %% Calculate objective values
        function PopObj = CalObj(obj, PopDec)
            % Compute dynamic objectives based on time t
            t = obj.t;
            G = abs(sin(0.5*pi*t));
            g=1+sum( (PopDec(:,2:end)-G).^2 ,2);
            PopObj(:,1) = g.*PopDec(:,1) + G;
            PopObj(:,2) = g.*(1-PopDec(:,1))+G;
            PopObj(PopObj<1e-18)=0;
        end

        %% Calculate constraint violations
        function PopCon = CalCon(obj, PopDec, PopObj)
            % Compute dynamic constraints
            t = obj.t;
            G = abs(sin(0.5*pi*t));
            c1 = cos(-0.15*pi)*PopObj(:,2) - sin(-0.15*pi)*PopObj(:,1) - (2*sin(5*pi*(sin(-0.15*pi)*PopObj(:,2) + cos(-0.15*pi)*PopObj(:,1)))).^6;
            PopCon =-c1;
        end

        %% Evaluate population
        function Population = Evaluation(obj, PopDec)
            % Evaluate solutions and update time parameter
            PopDec = obj.CalDec(PopDec);
            PopObj = obj.CalObj(PopDec);
            PopCon = obj.CalCon(PopDec, PopObj);
            Population = SOLUTION(PopDec, PopObj, PopCon, repmat(obj.FE, size(PopDec, 1), 1));
            obj.FE = obj.FE + length(Population);
            
            obj.t=min(floor(max(obj.FE-obj.preEvolution,0) / obj.N / obj.taut)/obj.nt,obj.max_time);
        end

        %% Generate points on the Pareto front
        function R = GetOptimum(obj, N)
            R=obj.Optimums{:,2};
            R=cat(1,R);
        end

        %% Generate the feasible region
        function R = GetPF(obj)
            R=obj.FeasibleRegions;
        end
        %% Calculate the score of a Population or M_score of Populations
        function score = CalMetric(obj,metName,Population)
            look =Population.adds;
            if length(look(1,:))==2
                Pa=Population.adds;
                t= min(floor(max(Pa(:,2)-obj.preEvolution,0)/obj.N/obj.taut)/obj.nt,obj.max_time);
            else
                t= min(floor(max(Population.adds-obj.preEvolution,0)/obj.N/obj.taut)/obj.nt,obj.max_time);
            end
            
            change = [0;find(t(1:end-1)~=t(2:end));length(t)];
            Scores = zeros(1,length(change)-1);
            allt   = cell2mat(obj.Optimums(:,1));

            for i = 1 : length(change)-1
                
                subPop    = Population(change(i)+1:change(i+1));
                Scores(i) = feval(metName,subPop,obj.Optimums{find(t(change(i)+1)==allt,1),2});
            end
            if obj.caonima && length(Scores)>1
                disp(class(obj));
                obj.shabi=[obj.shabi;Scores];
                if size(obj.shabi,1)==5
                    disp(mean(obj.shabi,1));
                    obj.shabi=[];
                end
            end
            score = mean(Scores);
        end
        %% Display a Population or all Populations with feasibile region
        function DrawObj(obj,Population,test_POF)
            look =Population.adds;
            if length(look(1,:))==2
                Pa=Population.adds;
                t= min(floor(max(Pa(:,2)-obj.preEvolution,0)/obj.N/obj.taut)/obj.nt,obj.max_time);
            else
                t= min(floor(max(Population.adds-obj.preEvolution,0)/obj.N/obj.taut)/obj.nt,obj.max_time);
            end
            change = [0;find(t(1:end-1)~=t(2:end));length(t)];
            allt   = cell2mat(obj.Optimums(:,1));
            tempStream = RandStream('mlfg6331_64','Seed',2);
            if nargin==2
                test_POF=0;
            end
            deviation_index=0;
            if obj.M==2
                for i = 1 : length(change)-1
                    color = rand(tempStream,1,3);
                    if test_POF && length(change)>2
                        
                    else
                        if length(change)==2
                            ax=Draw(Population(change(i)+1:change(i+1)).objs,'o','MarkerSize',5,'Marker','o','Markerfacecolor',sqrt(color),'Markeredgecolor',color,{'\it f\rm_1','\it f\rm_2',[]});
                        else
                            ax=Draw(Population(change(i)+1:change(i+1)).objs+deviation_index*t(change(i)+1),'o','MarkerSize',5,'Marker','o','Markerfacecolor',sqrt(color),'Markeredgecolor',color,{'\it f\rm_1','\it f\rm_2',[]});
                        end

                    end
                    if length(change)==2
                        Draw(obj.Optimums{find(t(change(i)+1)==allt,1),2},'o','MarkerSize',1,'Markerfacecolor',sqrt(color),'Markeredgecolor',color);
                    else
                        Draw(obj.Optimums{find(t(change(i)+1)==allt,1),2}+deviation_index*t(change(i)+1),'o','MarkerSize',1,'Markerfacecolor',sqrt(color),'Markeredgecolor',color);
                    end

                end
                
            end
            if length(change)==2
                if ~isempty(obj.PF)
                    PF=obj.PF{find(t(1)==allt,1),2};
                    if ~iscell(obj.PF)
                        if obj.M == 2
                            plot(ax,PF(:,1),PF(:,2),'-k','LineWidth',1);
                        elseif obj.M == 3
                            plot3(ax,PF(:,1),PF(:,2),PF(:,3),'-k','LineWidth',1);
                        end
                    else
                        if obj.M == 2
                            surf(ax,PF{1},PF{2},PF{3},'EdgeColor','none','FaceColor',[.85 .85 .85]);
                        elseif obj.M == 3
                            surf(ax,PF{1},PF{2},PF{3},'EdgeColor',[.8 .8 .8],'FaceColor','none');
                        end
                        set(ax,'Children',ax.Children(flip(1:end)));
                    end
                end
            end
        end
    end

end