% Main Algorithm for NSGA3
Offspring  = repmat(empty_individual,nPop,1);

[pop,FrontNo,CrowdDis,~]        = EnvironmentalSelectionNSGA2([pop],nPop);

%PopCon          = reshape(extractfield(pop,'Constraint'),[],size(pop,1))';
MatingPool      = TournamentSelection(2,nPop,FrontNo,-CrowdDis);
Offspr          = GA(pop(MatingPool),NSGA2);

%% Constraint handling needs to go in MatingPool selection - see Aaron's version of NSGA2
% is added in EnvironmentalSelection (in initialisation) and in the
% calculation of the mating pool here as well

% calculate the function values here

parfor i =1:nPop
    switch runcase
        case('Tom')
            [Cost(i,:),Constraint(i,:)]     = CostFunction(Offspr(i,:),i,SETTINGS)
        case('Suzie')
            [Cost(i,:),Constraint(i,:)]     = CostFunction(Offspr(i,:),i,OptSettings);
        case('Lachlan')
            [Cost(i,:),Constraint(i,:)]         =CostFunction(Offspr(i,:),i,SETTINGS);
    end
end

for cntr=1:nPop
    Offspring(cntr).Position    = Offspr(cntr,:);
    Offspring(cntr).Cost        = Cost(cntr,:);
    Offspring(cntr).Constraint  = Constraint(cntr,:);
end

[pop,FrontNo,CrowdDis,PF]        = EnvironmentalSelectionNSGA2([pop;Offspring],nPop);

PopObj          = reshape(extractfield(pop,'Cost'),[],size(pop,1))';
PopCon          = reshape(extractfield(pop,'Constraint'),[],size(pop,1))';
ConVPop     = sum(max(0, PopCon), 2);
index = find(ConVPop == 0);
feasibleFraction = length(index)/length(ConVPop);

% Estimate simulation time 
if it == 1
    first = datetime('now');
elseif it == 2
    second = datetime('now');
    iterationTime = diff([first, second]);
    simTime = iterationTime*MaxIt; 
    disp(['Expected Simulation Duration: ', num2str(hours(simTime)),' hours'])
end

% Show Iteration Information
disp(['NSGA2  Iteration ' num2str(it) ': Nr of Pareto Members = ' num2str(PF) ' Feasible fraction ' num2str(round(feasibleFraction,2)) ' Finished at ' datestr(now)]);
% Damping Inertia Weight
if Optimisation.plot == 1
    CostPlot = reshape(extractfield(pop,'Cost'),[],size(pop,1))';
    scatter(CostPlot(:,1),CostPlot(:,2))
    pause(0.05)
    drawnow
end