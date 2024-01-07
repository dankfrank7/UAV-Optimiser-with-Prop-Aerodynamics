%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          - TESTING SANDBOX -                            % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
close all
cd('X:\Users\Tom\Documents\Uni Work\Y5S1\Thesis\aMainFiles\VAP3point5')
addpath('X:\Users\Tom\Documents\Uni Work\Y5S1\Thesis\aMainFiles')
InitialisePlotting();
set(gcf,'color','w');

% SETTINGS 
toSave = 1; % Saves plots to file 
VAP_IN = []; 
VAP_IN.PRINT = 1;
VAP_IN.PLOT = 1;
VAP_IN.GIF = 1;
% Optimiser Settings 
SETTINGS.AnalyseStability   = true;
SETTINGS.stabSelection      = [1 1 1]; % Cma Clb Cnb (Analyse = 1, Don't analyse = 0)
SETTINGS.Wingspan_constraint = 0.45;
SETTINGS.Length_constraint   = 0.45;
SETTINGS.Cma_lower          = -1.5;
SETTINGS.Cma_upper          = -0.05;
SETTINGS.Clb_lower          = -0.3;
SETTINGS.Clb_upper          = -0.05;
SETTINGS.CNb_lower          = 0.05;
SETTINGS.CNb_upper          = 0.6;
SETTINGS.Cd_extra           = 0; % CRUD drag 
SETTINGS.rpmMax             = 12e4; % rpm
% Set
SETTINGS.relax          = 'TRUE';
SETTINGS.steady         = 'TRUE';
SETTINGS.elecWeight     = 0.7;
SETTINGS.airfoilCG      = 0.35;
SETTINGS.systemWeight   = 0.22;
SETTINGS.cogExplicitx   = nan;%.0566; %nan to calculate based off geometry.
SETTINGS.cogExplicitz   = nan;
SETTINGS.weightExplicit = nan; % nan to calculate based off geometry.
% Conditions s
SETTINGS.speed          = 25;   %'nan' to turn on fixed Lift analysis
SETTINGS.density        = 1.225;
SETTINGS.beeta          = 0;    % Sideslipe angle
% Wing
SETTINGS.planformType   = 'zimmerman'; 
SETTINGS.wingAirfoil    = 'Phoenix';
SETTINGS.wingSpanSpacing = 'halfsine'; % or linear
SETTINGS.noWingSpanEle  = 10;
SETTINGS.noWingDVE      = 2*(SETTINGS.noWingSpanEle-1);
SETTINGS.symmetry       = 'FALSE';
SETTINGS.incidence      = 0;
SETTINGS.chordwise_elements = 1;
% Prop
SETTINGS.collective     = 0;
SETTINGS.propellerName  = 'GWS4x4'; % noProp to turn off prop
SETTINGS.blades         = '2';
SETTINGS.noRotorSpanEle = 6;
SETTINGS.pitchShift     = 2;
SETTINGS.propDiam       = 0.23;
SETTINGS.propAirfoil    = 'naca4409';
SETTINGS.ref_diam       = SETTINGS.propDiam;
SETTINGS.rotation_direction = 'CCW';
SETTINGS.veh_x_hub      = -0.05;
SETTINGS.veh_z_hub      = -0.03;
%  Vertical Stabiliser     
SETTINGS.wingtips = 0;
SETTINGS.vt_doublePanel = 0;
SETTINGS.taper_vt       = 0.6; % Taper from trailing edge
SETTINGS.verticalTailType = 'rectangular'; % 'noVerticalTail' removes tail
SETTINGS.noVTSpanEle    = 6;
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
SETTINGS.horizontalTailType = 'rectangular'; % 'noHorizontalTail' or 'rectangular'
SETTINGS.taper_ht       = 0.536;
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
SETTINGS.noFuSpanEle    = 4;
SETTINGS.fuseAR         = 0.5; % This Overrites fusediam 
SETTINGS.taper_fu       = 0.8;
SETTINGS.sweep_fu       = 50;
SETTINGS.xStart_fu      = -0.03;
SETTINGS.zStart_fu      = -0.03; % Nan to define based on croot and fuseAR
% Variable waiting room 
SETTINGS.rpm        = 60*SETTINGS.speed./0.7/SETTINGS.propDiam; % Change this to advance ratio?
SETTINGS.taper         = 0.4;
SETTINGS.volumeHtail   = 0.3;
SETTINGS.ar_ht         = 1;
SETTINGS.l_vt          = 1.6;% x-position of leading edge (distance from ac_w to ac_vt)
% VARIABLES (NOTE: actual value here does not matter)
VARIABLES.taper_vt      = 0.5;
VARIABLES.taper_vt2     = 0.7;
VARIABLES.volumeTail    = 0.10;% Volume Coefficient
VARIABLES.volumeTail2   = 0.15;% Lower tail 
VARIABLES.ar_vt         = 1.2;    
VARIABLES.ar_vt2        = 0.5;
VARIABLES.sweep_vt       = nan; % nan changes sweep to keep trailing edge straight
VARIABLES.sweep_vt2     = nan;
VARIABLES.dihedral      = 7;
VARIABLES.alfa          = 0; 
VARIABLES.AR            = 4.8;
VARIABLES.c_root        = 0.125;
VARIABLES.sweep         = 0;
VARIABLES.twist         = 0;


% Index 
wing_index = 1:SETTINGS.noWingDVE;
start_ht = SETTINGS.noWingDVE + 1;
end_ht   = SETTINGS.noWingDVE + SETTINGS.noHTDVE;
ht_index = start_ht:end_ht;

% [SETTINGS.varMin, SETTINGS.varMax, SETTINGS.variableNames, nVar, VARIABLES_norm] = fcnVARIABLE_SETUP(VARIABLES);
% 
% hFig = fcnPLOT_GEOM(VARIABLES_norm, SETTINGS);
[SETTINGS.varMin, SETTINGS.varMax, SETTINGS.variableNames, nVar, VARIABLES_norm] = fcnVARIABLE_SETUP(VARIABLES);

ID = 0;
free = 0;
%% ITERATION STUDY
% it_vec = [10, 15, 20, 25, 30, 40];
% for i = 1:length(it_vec)
%     SETTINGS.maxtime        = it_vec(i);
%     SETTINGS.maxTimeStab    = it_vec(i);
%     
%     SETTINGS.delta_time     = 0.001;
% 
%     [inputFilename, GEOMETRY, WEIGHT] = fcnGENERATE_INPUT(VARIABLES_norm, ID, SETTINGS, free);
%     
%     [OUTP, COND, ~, ~, ~, SURF, ~, ~, ~] = fcnVAP_MAIN(inputFilename, VAP_IN);  
%     %STABILITY = fcnSTABILITY(inputFilename, GEOMETRY, SETTINGS, VAP_IN, OUTP, SURF, COND, SETTINGS.maxTimeStab, ID);
%     
%   
%     Cma_vec(i) = STABILITY.Cma;
%     Cnb_vec(i) = STABILITY.CNb;
%     Clb_vec(i) = STABILITY.Clb;
% 
% end

%figure
% plot(it_vec,Cma_vec,'-o')
% ylabel('Cma')
% xlabel('Max Iterations')
% 
% figure
% plot(it_vec,Clb_vec,'-o')
% ylabel('Clb')
% xlabel('Max Iterations')
% 
% figure
% plot(it_vec,Cnb_vec,'-o')
% ylabel('Cnb')
% xlabel('Max Iterations')

%% NORMAL FORCE
it_vec = [5, 10, 15, 20, 30, 40, 50, 60, 70, 80];
it_vec = 40;
it_vec_wing = repelem(it_vec,length(wing_index));
it_vec_tail = repelem(it_vec,length(ht_index));
thetay_w      = linspace(0,pi,SETTINGS.noWingSpanEle);
y_space_w     = sin(thetay_w/2);
thetay_t      = linspace(0,pi,SETTINGS.noHTSpanEle);
y_space_t     = sin(thetay_t/2);

dAngle_vec = [18];
for j = 1:length(dAngle_vec)
    for i = 1:length(it_vec)
        SETTINGS.maxtime        = it_vec(i);
        SETTINGS.maxTimeStab    = it_vec(i);

         omega = SETTINGS.rpm*2*pi/60;
         dt    = deg2rad(dAngle_vec(j))/omega; 
         SETTINGS.delta_time     = dt;

        [inputFilename, GEOMETRY, WEIGHT] = fcnGENERATE_INPUT(VARIABLES_norm, ID, SETTINGS, free);

        [OUTP, COND, ~, ~, ~, SURF, ~, ~, ~] = fcnVAP_MAIN(inputFilename, VAP_IN);  
        STABILITY = fcnSTABILITY(inputFilename, GEOMETRY, SETTINGS, VAP_IN, OUTP, SURF, COND, SETTINGS.maxTimeStab, ID);

        wingCn(i,:) = OUTP.vecCNDIST(wing_index);
        tailCn(i,:) = OUTP.vecCNDIST(ht_index);
        
        Cma_vec(i) = STABILITY.Cma;
        Cnb_vec(i) = STABILITY.CNb;
        Clb_vec(i) = STABILITY.Clb;

        figure(j*5)
        plot3(ones(1,length(wing_index))*it_vec(i),[-flip(y_space_w(1:end-1)), y_space_w(1:end-1)],wingCn(i,:),'-o')
        hold on 
        grid on 
        box on
        xlabel('Max Iteration'); ylabel('2y/b'); zlabel('Wing C_N');
        drawnow;
        figure(j*6)
        plot3(ones(1,length(ht_index))*it_vec(i),[-flip(y_space_t(1:end-1)), y_space_t(1:end-1)],tailCn(i,:),'-o')
        hold on 
        grid on 
        xlabel('Max Iteration'); ylabel('2y/b'); zlabel('Tail C_N');
        drawnow;
    end
    figure(j*7)
    plot(it_vec,Cma_vec,'-o')
    ylabel('Cma')
    xlabel('Max Iterations')

    figure(j*8)
    plot(it_vec,Clb_vec,'-o')
    ylabel('Clb')
    xlabel('Max Iterations')

    figure(j*9)
    plot(it_vec,Cnb_vec,'-o')
    ylabel('Cnb')
    xlabel('Max Iterations')
end





%% PROP ANGLE 
% dAngle_vec = [5, 15, 25, 35];
% 
% for i = 1:length(dAngle_vec)
%     SETTINGS.maxtime        = 25;
%     SETTINGS.maxTimeStab    = 25;
%     
%     omega = SETTINGS.rpm*2*pi/60;
%     dt    = deg2rad(dAngle_vec(i))/omega; 
%     SETTINGS.delta_time     = dt;
% 
%     [inputFilename, GEOMETRY, WEIGHT] = fcnGENERATE_INPUT(VARIABLES_norm, ID, SETTINGS, free);
%     
%     [OUTP, COND, ~, ~, ~, SURF, ~, ~, ~] = fcnVAP_MAIN(inputFilename, VAP_IN);  
%     STABILITY = fcnSTABILITY(inputFilename, GEOMETRY, SETTINGS, VAP_IN, OUTP, SURF, COND, SETTINGS.maxTimeStab, ID);
%     
%     
%     Cma_vec(i) = STABILITY.Cma;
%     Cnb_vec(i) = STABILITY.CNb;
%     Clb_vec(i) = STABILITY.Clb;
% end
% 
% figure
% plot(dAngle_vec,Cma_vec,'-o')
% ylabel('Cma')
% xlabel('Prop Rotation Angle')
% 
% figure
% plot(dAngle_vec,Clb_vec,'-o')
% ylabel('Clb')
% xlabel('Prop Rotation Angle')
% 
% figure
% plot(dAngle_vec,Cnb_vec,'-o')
% ylabel('Cnb')
% xlabel('Prop Rotation Angle')
