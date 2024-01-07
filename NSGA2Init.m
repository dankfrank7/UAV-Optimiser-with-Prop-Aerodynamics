% Initialisation for NSGA3
% Algorithm settings
NSGA2.nPop                = nPop;                 % Population Size
NSGA2.nRep                = nPop;                 % Repository Size
NSGA2.nObj                = nObj;                 % Number of Objectives of the Problem

[NSGA2.proC,NSGA2.disC,NSGA2.proM,NSGA2.disM] = deal(1,20,1,20);

NSGA2.lower               	= zeros(1,nVar);                % normalised bounds lower bound is 0
NSGA2.upper                 = ones(1,nVar);                 % normalised bounds upper bound is 1

% Initialisation of the individual
empty_individual.Position   = [];
empty_individual.Cost       = [];
empty_individual.Constraint = [];

