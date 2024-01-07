function [Population,FrontNo,CrowdDis,PF] = EnvironmentalSelectionNSGA2(Population,N)
% The environmental selection of NSGA-II

PopObj    = reshape(extractfield(Population,'Cost'),[],size(Population,1))';
PopCon    = reshape(extractfield(Population,'Constraint'),[],size(Population,1))';

%% Non-dominated sorting
[FrontNo,MaxFNo]    = NDSort(PopObj,PopCon,N);
Next                = FrontNo < MaxFNo;

%% Calculate the crowding distance of each solution
CrowdDis            = CrowdingDistance(PopObj,FrontNo);

%% Select the solutions in the last front based on their crowding distances
Last                = find(FrontNo==MaxFNo);
[~,Rank]            = sort(CrowdDis(Last),'descend');
Next(Last(Rank(1:N-sum(Next)))) = true;

%% Population for next generation
Population  = Population(Next);
FrontNo     = FrontNo(Next);
CrowdDis    = CrowdDis(Next);
PF          = sum(FrontNo==1);

end