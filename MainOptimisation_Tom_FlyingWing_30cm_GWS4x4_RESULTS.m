%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   BWB AUV Optimization                                                  %
%   Author: Dries Verstraete                                              %
%   Date: August 2019                                                     %
%   Organisation: Usyd                                                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;
tic
format compact
warning('off','all')

%addpath('.\VAP3point5')
curr_dir = pwd
cd './VAP3point5'
addpath(curr_dir)


disp(['Calculation Started at ' datestr(now)]);

% Delete old files:
delete *.inp *.ngh *.pan *.par *.czy *.out *.sdsa *.txt *.ms2 *.err *.f *.dat *.prs *.f06 *.op2 *IFPDAT *.bdf 

%% This is swhere your results will be save
saveFolder   = strcat(curr_dir,'/OUTPUT/');
saveFilename = 'FlyingWing_NoProp_30span_DRIES_RESULTS'; % Change name of results file HERE!!!!
% 'TestVariablesDate.mat'; 

%
% NOTES: 30 cm class based off PMAV + Propeller
%

%% Flight Condition

% SETTINGS 
toSave = 1; % Saves plots to file 
VAP_IN = []; 
% VAP_IN.PRINT = 1;
% VAP_IN.PLOT = 1;
% VAP_IN.GIF = 1;
% Optimiser Settings 
SETTINGS.AnalyseStability   = true;
SETTINGS.stabSelection      = [1 1 1]; % Cma Clb Cnb (Analyse = 1, Don't analyse = 0)
SETTINGS.Wingspan_constraint = 0.30;
SETTINGS.Length_constraint   = 0.30;
SETTINGS.Cma_lower          = -1.5;
SETTINGS.Cma_upper          = -0.05;
SETTINGS.Clb_lower          = -0.3;
SETTINGS.Clb_upper          = -0.05;
SETTINGS.CNb_lower          = 0.05;
SETTINGS.CNb_upper          = 0.6;
SETTINGS.Cd_extra           = 0; % CRUD drag 
SETTINGS.rpmMax             = 12e4; % rpm
% Set
SETTINGS.maxtime        = 40;
SETTINGS.maxTimeStab    = 40;
SETTINGS.dAngle         = 18; % Not used in non-prop cases
SETTINGS.relax          = 'TRUE';
SETTINGS.steady         = 'TRUE';
SETTINGS.systemWeight   = 0.08;
SETTINGS.elecWeight     = 0.7;
SETTINGS.airfoilCG      = 0.4;
SETTINGS.cogExplicitx   = nan;%.0566; %nan to calculate based off geometry.
SETTINGS.cogExplicitz   = nan;
SETTINGS.weightExplicit = nan; % nan to calculate based off geometry.
% Conditions s
SETTINGS.speed          = 18;   %'nan' to turn on fixed Lift analysis
SETTINGS.density        = 1.225;
SETTINGS.beeta          = 0;    % Sideslipe angle
% Wing
SETTINGS.planformType   = 'zimmerman'; 
SETTINGS.wingAirfoil    = 'Phoenix';
SETTINGS.wingSpanSpacing = 'halfsine'; % or linear
SETTINGS.noWingSpanEle  = 8;
SETTINGS.noWingDVE      = 2*(SETTINGS.noWingSpanEle-1);
SETTINGS.symmetry       = 'FALSE';
SETTINGS.incidence      = 0;
SETTINGS.chordwise_elements = 1;
% Prop
SETTINGS.collective     = 0;
SETTINGS.propellerName  = 'GWS4x4'; % noProp to turn off prop
SETTINGS.blades         = '2';
SETTINGS.noRotorSpanEle = 6;
SETTINGS.pitchShift     = 8;
SETTINGS.propDiam       = 0.1016;
SETTINGS.propAirfoil    = 'naca4409';
SETTINGS.ref_diam       = SETTINGS.propDiam;
SETTINGS.rotation_direction = 'CCW';
SETTINGS.veh_x_hub      = -0.05;
SETTINGS.veh_z_hub      = -0.03;
%  Vertical Stabiliser     
SETTINGS.wingtips       = false;
SETTINGS.vt_doublePanel = true;
SETTINGS.verticalTailType = 'rectangular'; % 'noVerticalTail' removes tail
SETTINGS.noVTSpanEle    = 5;
SETTINGS.z_vtStart      = 0;
SETTINGS.chordwise_elements_vt = 1;
SETTINGS.VerticalTail_Airfoil = 'ht08';
if SETTINGS.vt_doublePanel
    SETTINGS.noVTDVE        = 2*SETTINGS.noVTSpanEle - 2;
else 
    SETTINGS.noVTDVE        = SETTINGS.noVTSpanEle - 1;
end
% Horizontal Stabiliser     
SETTINGS.sweep_ht       = nan; % nan changes sweep to keep trailing edge straight
SETTINGS.horizontalTailType = 'noHorizontalTail'; % 'noHorizontalTail' or 'rectangular'
SETTINGS.taper_ht       = 0.5;
SETTINGS.twist_ht       = 0;
SETTINGS.dihedral_ht    = 0;
SETTINGS.zStart_ht      = 0;            % Absolute
SETTINGS.incidence_ht   = 0;
SETTINGS.HorizontalTail_Airfoil = 'ht08';
SETTINGS.chordwise_elements_ht = 1;
SETTINGS.noHTSpanEle    = 6;
SETTINGS.noHTDVE        = 2*(SETTINGS.noHTSpanEle -1);
SETTINGS.HTSpanSpacing  = 'halfsine';
% Fuselage 
SETTINGS.fuselageType   = 'rectangular'; % NoFuselage or rectangular
SETTINGS.incidence_fu   = 0;
SETTINGS.chordwise_elements_fu = 1;
SETTINGS.Fuselage_airfoil = 'NACA0012';
SETTINGS.FuSpanSpacing  = 'linear';
SETTINGS.noFuSpanEle    = 2;
SETTINGS.fuseAR         = 0.2; % This Overrites fusediam 
SETTINGS.taper_fu       = 0.8;
SETTINGS.sweep_fu       = 50;
SETTINGS.xStart_fu      = -0.03;
SETTINGS.zStart_fu      = nan; % Nan to define based on croot and fuseAR
% Prop Angle DT 
SETTINGS.rpm        = 60*SETTINGS.speed./0.7/SETTINGS.propDiam;
omega = SETTINGS.rpm*2*pi/60;
dt    = deg2rad(SETTINGS.dAngle)/omega; 
SETTINGS.delta_time = dt;
% Variable waiting room 
SETTINGS.taper         = 0.4;
SETTINGS.volumeHtail   = 1;
SETTINGS.ar_ht         = 1;
SETTINGS.l_vt          = 1.05;% x-position of leading edge * c_root
SETTINGS.sweep_vt      = nan; % nan changes sweep to keep trailing edge straight
SETTINGS.sweep_vt2     = nan;
SETTINGS.twist         = 2;
% VARIABLES (NOTE: actual value here does not matter)
VARIABLES.taper_vt      = 0.5;
VARIABLES.taper_vt2     = 0.7;
VARIABLES.volumeTail    = 0.15;% Volume Coefficient
VARIABLES.volumeTail2   = 0.15;% Lower tail 
VARIABLES.ar_vt         = 0.5;    
VARIABLES.ar_vt2        = 0.5;

VARIABLES.dihedral      = 0;
VARIABLES.alfa          = 4; 
VARIABLES.AR            = 0.8;
VARIABLES.c_root        = 0.2;
VARIABLES.sweep         = 15;


[SETTINGS.varMin, SETTINGS.varMax, SETTINGS.variableNames, SETTINGS.nVar] = fcnVARIABLE_SETUP(VARIABLES);

%% Some Optimisation Settings
Optimisation.MaxNrOfIterations      = 100; % INCREASE to 100
Optimisation.NrOfSwitches           = 10;
Optimisation.PopulationSize         = 8;       % this gets multiplied by nVar!!!
Optimisation.initialdesign          = 'latin';  % for now only latin hypercube and sobol implement

% for now only use NSGA2 / NSGA3 / C2ODE
% only run one algorithm
Optimisation.Algorithm              = 'NSGA2'; % NSGA3 CMOPSO MOEADDE GDE3 NSGA3

Optimisation.run                    = 'Tom';
Optimisation.ConstraintHandling     = 'iE';
Optimisation.switchiteration        = 0.7; % fraction of total iterations where switch to second algorithm is made
Optimisation.MaxItLimit             = 1.01; % (fraction) can go up to this many iterations if pareto front does not have enough members
Optimisation.FeasibleCutOff         = 0.8;

% Add your settings here
OptSettings.plot                    = 0;
Optimisation.plot                   = OptSettings.plot;

%% Add all your settings stuff in here - paste lines 10 to 64 of thesis main script in here

%% Problem Definition
runcase = Optimisation.run;
switch runcase
    
    case 'Tom'
        CostFunction=@(x,Id,SETTINGS) fcnMAINSOLVER_TOM(x,Id,SETTINGS,VAP_IN);
        
        nVar                    = SETTINGS.nVar;                   % Number of Decision Variables
        x                       = unifrnd(0,1,1,nVar);
        Id                      = 1;
        
        OptSettings.nObj       = numel(CostFunction(x,Id,SETTINGS));
        
    case('Lachlan')
        
        CostFunction=@(x,Id,OptSettings) fcnMAINSOLVER(x,Id,SETTINGS);  % Cost Function
        % Optimisation Setup - Define all the variables in here
        

        nVar                    = length(VARIABLES);                   % Number of Decision Variables
        x                       = unifrnd(0,1,1,nVar);
        Id                      = 1;
        
        OptSettings.nObj       = numel(CostFunction(x,Id,SETTINGS));
        % OptSettings.nObj       = numel(CostFunction(x,AoA,density,velocity,fields,VarMin,VarMax,Id,OptSettings));
        
    case('Suzie')
        
        CostFunction=@(x,Id,SETTINGS) fcnMAINSOLVER(x,Id,SETTINGS);  % Cost Function
        % Optimisation Setup - Define all the variables in here
        OptimisationSettings;
        
        VarMin     = OptSettings.VarMin;
        VarMax     = OptSettings.VarMax;
        nVar                    = length(fieldnames(VarMin));                   % Number of Decision Variables
        Optimisation.Variables  = fieldnames(VarMin);
        for ijk = 1:nVar
            OptSettings.LowerBound(ijk) = VarMin.(Optimisation.Variables{ijk});
            OptSettings.UpperBound(ijk) = VarMax.(Optimisation.Variables{ijk});
        end
        fields = fieldnames(VarMin);
        OptSettings.fields     = fields;
        x                       = unifrnd(0,1,1,nVar);
        Id                      = 1;
        
        OptSettings.nObj       = numel(CostFunction(x,Id,OptSettings));
        % OptSettings.nObj       = numel(CostFunction(x,AoA,density,velocity,fields,VarMin,VarMax,Id,OptSettings));
  
end

OptSettings.nVar       = nVar;
VarSize                 = [1 OptSettings.nVar];   % Size of Decision Variables Matrix

% Number of Objective Functions
Id = 0;
nVar                    = OptSettings.nVar;
nObj                    = OptSettings.nObj;



%% Optimisation Parameters
MaxIt                   = Optimisation.MaxNrOfIterations                    % Maximum Number of Iterations (number of generations)

%% Deal with multiple algorithms
if isfield(Optimisation,'Algorithm2')
    if isfield(Optimisation,'Algorithm3')
        alglist = {Optimisation.Algorithm,Optimisation.Algorithm2,Optimisation.Algorithm3};
    else
        alglist = {Optimisation.Algorithm,Optimisation.Algorithm2};
    end
else
    alglist={Optimisation.Algorithm};
end

if ~isempty(find(strcmp(alglist,'C2ODE')))
    nPop                = 50
else
    nPop               	= Optimisation.PopulationSize*nVar                	% Population Size (10 times number of design variables)

end

%% Initialization
% Set-up parallelisation:
% Determine number of cores available:
numcores = 0; % Set to 2 (will have 2 available to work)
if numcores == 0
    numcores = feature('numcores');
end

% checking to see if my pool is already open
if numcores > 1
    if isempty(gcp('nocreate')) == 1
        % open parallel pool:
        parpool(numcores);
        pctRunOnAll warning('off','all')
    end
end

conv_histCost               = [];
conv_histPos                = [];
conv_histConstraint         = [];

it = 0;

% This is where all the algorithms get initialised. If you want to add a
% new algorithm set up its initialisation here
InitialisationSequence

if rem(nPop,2)~=0
nPop = nPop+1;
Optimisation.nPop = nPop;
end

% Initialise the population as an empty population
pop                         = repmat(empty_individual,nPop,1);

feasibleFraction            = 0;
%% Start the initialisation
disp(['Initialisation Started at ' datestr(now)]);

% This is where the initial design sequences are set up
% For now it has a latin hypercube and Sobol sequence.
% If you want to add others this is the place to do it
InitialDesigns

for i=1:nPop
    Position(i,:)                               = InitialGeneration(i,:);
end

% And calculate the value of the cost function

parfor i=1:nPop
    switch runcase
        case('Tom')
            [Cost(i,:),Constraint(i,:)]         =CostFunction(Position(i,:),i,SETTINGS)
        case('Suzie')
            [Cost(i,:),Constraint(i,:)]         =CostFunction(Position(i,:),i,OptSettings);
        case('Lachlan')
            [Cost(i,:),Constraint(i,:)]         =CostFunction(Position(i,:),i,SETTINGS);
    end
end


for i=1:nPop
    pop(i).Position     = Position(i,:);
    pop(i).Cost         = Cost(i,:);
    pop(i).Constraint   = Constraint(i,:);
end

% Save convergence history:
conv_histCost                       = [conv_histCost;...
    reshape(extractfield(pop,'Cost'),[],size(pop,1))'  it*ones(size(pop,1),1)];
conv_histPos               = [conv_histPos;...
    reshape(extractfield(pop,'Position'),[],size(pop,1))' it*ones(size(pop,1),1)];
conv_histConstraint         = [conv_histConstraint;...
    reshape(extractfield(pop,'Constraint'),[],size(pop,1))' it*ones(size(pop,1),1)];


if OptSettings.plot == 1
    figure(1)
    hold on
    grid on
    figure(2)
    hold on
    grid on
    %axis([0 1 0 2])
end

AlgSwitch           = 0;
RestartSwitch       = 0;
MaxItLim            = Optimisation.MaxItLimit*MaxIt;

% Some algorithms need an initial calculation of a value to start things
% of. Rather than implementing this in the main file they are set up in
% this separate file for clarity
InitialCalculations

SwitchFrequency = ceil(Optimisation.MaxNrOfIterations./Optimisation.NrOfSwitches);


%% Main Loop
while it<MaxIt+1
        if rem(it,50) == 49
             delete(gcp('nocreate'));
             parpool(numcores);
             pctRunOnAll warning('off','all');
         end
    % Delete old files:
    delete *.inp *.ngh *.pan *.par *.czy *.out *.sdsa *.txt *.ms2 *.err *.f *.dat *.prs *.f06 *.op2 *IFPDAT *.bdf
    if rem(it,SwitchFrequency) == 1
        if it < 0.9*MaxIt
            if length(alglist) > 1
                l = rand;
                if length(alglist) ==2
                    if l < 0.5
                        Optimisation.Algorithm = alglist{1};
                    else
                        Optimisation.Algorithm = alglist{2};
                    end
                else
                    if l <= 1/3
                        Optimisation.Algorithm = alglist{1};
                    elseif l<=2/3
                        Optimisation.Algorithm = alglist{2};
                    else
                        Optimisation.Algorithm = alglist{3};
                    end
                end
            else
                Optimisation.Algorithm = alglist{1};
            end
        elseif ~isempty(find(strcmp(alglist,'NSGA3')))
            Optimisation.Algorithm = 'NSGA3';
        else
            if length(alglist) > 1
                l = rand;
                if length(alglist) ==2
                    if l < 0.5
                        Optimisation.Algorithm = alglist{1};
                    else
                        Optimisation.Algorithm = alglist{2};
                    end
                else
                    if l <= 1/3
                        Optimisation.Algorithm = alglist{1};
                    elseif l<=2/3
                        Optimisation.Algorithm = alglist{2};
                    else
                        Optimisation.Algorithm = alglist{3};
                    end
                end
            else
                Optimisation.Algorithm = alglist{1};
            end
            
        end
    end
    
    % This is the main switch to call the various algorithms
    % Add the call to the main part of the algorithm in here
    MainCalculations
    
    it = it+1;
    switch Optimisation.ConstraintHandling
        case('iE2')
            iE2.it          = it;
    end
    % Save convergence history:
    conv_histCost       = [conv_histCost; reshape(extractfield(pop,'Cost'),[],size(pop,1))'  it*ones(size(pop,1),1)];
    conv_histPos        = [conv_histPos; reshape(extractfield(pop,'Position'),[],size(pop,1))' it*ones(size(pop,1),1)];
    conv_histConstraint = [conv_histConstraint; reshape(extractfield(pop,'Constraint'),[],size(pop,1))' it*ones(size(pop,1),1)];
    
    % Save workspace every 10 iterations 
    if rem(it,10) == 0
        LowerBounds = SETTINGS.varMin;
        UpperBounds = SETTINGS.varMax;
        save(strcat(saveFolder,saveFilename))
    end 
end

switch runcase
    case('Tom')
        LowerBounds = SETTINGS.varMin;
        UpperBounds = SETTINGS.varMax;
        
        % Do not overrite any output file 
%         if isfile(fullfile(saveFolder, strcat(saveFilename,'.mat')))
%             while isfile(fullfile(saveFolder, strcat(saveFilename,'.mat')))
%                 saveFilename = strcat(saveFilename,'_copy');
%             end
%             save(strcat(saveFolder,saveFilename))
%             fprintf('File duplicate detected, workspace saved as %s\n',saveFilename)
%         else
%             save(strcat(saveFolder,saveFilename))
%         end 
        save(strcat(saveFolder,saveFilename))
        
        % Plot output
        if OptSettings.plot 
            hFig = fcnPLOT_OUTPUT(conv_histCost, conv_histPos, conv_histConstraint, CostPlot, SETTINGS, MaxIt, saveFilename, toSave);
        end
        
    case('Lachlan')
        LowerBounds = SETTINGS.lower_bound;
        UpperBounds = SETTINGS.lower_bound;
        save(savefilename)
    case('Suzie')
        LowerBounds = OptSettings.LowerBound;
        UpperBounds = OptSettings.UpperBound;
        save(savefilename)
    case('test')
        figure
        plot(PFt(:,1),PFt(:,2));
        hold on
        grid on
        scatter(reshape(conv_histCost(:,1),[],1),reshape(conv_histCost(:,2),[],1),40,conv_histCost(:,3),'filled')
        colorbar
        scatter(CostPlot(:,1),CostPlot(:,2),'rx')
        axis([0 1 0 2])
    case('testc')
        figure
        plot(PFt(:,1),PFt(:,2));
        hold on
        grid on
        scatter(reshape(conv_histCost(:,1),[],1),reshape(conv_histCost(:,2),[],1),40,conv_histCost(:,3),'filled')
        colorbar
        scatter(CostPlot(:,1),CostPlot(:,2),'rx')
        axis([0 1.5 0 2])
end
toc
% Save results:
% Close parallel pool:
delete(gcp('nocreate'));

cd('X:\Users\Tom\Documents\Uni Work\Y5S1\Thesis\aMainFiles')
rmpath('./VAP3point5/')
return

