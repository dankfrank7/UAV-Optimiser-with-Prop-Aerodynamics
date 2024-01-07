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

VAP_IN = [];
VAP_IN.PRINT = 1;
VAP_IN.PLOT = 1;
VAP_IN.GIF = 1;
% Optimiser Settings 
SETTINGS.AnalyseStability   = false;
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
SETTINGS.rpmMax             = 15e4; % rpm
% Set
SETTINGS.maxtime        = 10;
SETTINGS.maxTimeStab    = 10;
SETTINGS.dAngle         = 18;
SETTINGS.relax          = 'TRUE';
SETTINGS.steady         = 'TRUE';
SETTINGS.systemWeight   = 0.22;
SETTINGS.elecWeight     = 0.7;
SETTINGS.airfoilCG      = 0.35;
SETTINGS.cogExplicitx   = nan;%.0566; %nan to calculate based off geometry.
SETTINGS.cogExplicitz   = nan;
SETTINGS.weightExplicit = nan; % nan to calculate based off geometry.
% Conditions s
SETTINGS.speed          = 18;   %'nan' to turn on fixed Lift analysis
SETTINGS.density        = 1.225;
SETTINGS.beeta          = 0;    % Sideslipe angle
% Wing
SETTINGS.ref_span       = 0.35;
SETTINGS.planformType   = 'rectangular'; 
SETTINGS.wingAirfoil    = 'Phoenix';
SETTINGS.wingSpanSpacing = 'halfsine'; % or linear
SETTINGS.noWingSpanEle  = 8;
SETTINGS.noWingDVE      = 2*(SETTINGS.noWingSpanEle-1);
SETTINGS.symmetry       = 'FALSE';
SETTINGS.incidence      = 0;
SETTINGS.chordwise_elements = 1;
% Prop
SETTINGS.collective     = 0;
SETTINGS.propellerName  = 'APC7x6'; % noProp to turn off prop
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
SETTINGS.wingtips       = true; % must have double panel for this 
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
SETTINGS.taper_ht       = 0.4;
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
SETTINGS.fuselageType   = 'rectangular'; % noFuselage or rectangular
SETTINGS.incidence_fu   = 0;
SETTINGS.chordwise_elements_fu = 1;
SETTINGS.Fuselage_airfoil = 'NACA0012';
SETTINGS.FuSpanSpacing  = 'linear';
SETTINGS.noFuSpanEle    = 2;
SETTINGS.fuseAR         = 0.2; % This Overrites fusediam 
SETTINGS.taper_fu       = 0.8;
SETTINGS.sweep_fu       = 50;
SETTINGS.xStart_fu      = -0.02;
SETTINGS.zStart_fu      = nan; % Nan to define based on croot and fuseAR
% Prop Iter Testing
SETTINGS.rpm        = 60*SETTINGS.speed./0.7/SETTINGS.propDiam;
omega = SETTINGS.rpm*2*pi/60;
dt    = deg2rad(SETTINGS.dAngle)/omega; 
SETTINGS.delta_time     = dt;

% Variable waiting room 
SETTINGS.taper         = 0.5;
SETTINGS.l_vt          = 1.05;% x-position of leading edge * c_root
SETTINGS.sweep_vt      = nan; % nan changes sweep to keep trailing edge straight
SETTINGS.sweep_vt2     = nan;


SETTINGS.taper_vt2     = 0.5;
SETTINGS.twist         = 0;

% VARIABLES (NOTE: actual value here does not matter)
VARIABLES.taper_vt      = 0.35;
VARIABLES.volumeTail    = 0.08;% Volume Coefficient
SETTINGS.volumeTail2    = 0.06;% Lower tail 
VARIABLES.ar_vt         = 1.4;    
SETTINGS.ar_vt2         = 1;

VARIABLES.volumeHtail   = 0.45;
VARIABLES.ar_ht         = 1.8;

VARIABLES.dihedral      = 10;
VARIABLES.alfa          = 2; 
VARIABLES.AR            = 3.8;
VARIABLES.c_root        = 0.15;
VARIABLES.sweep         = 25;

ID = 0;

% Optimisation variables set up 
[SETTINGS.varMin, SETTINGS.varMax, SETTINGS.variableNames, nVar, VARIABLES_norm] = fcnVARIABLE_SETUP(VARIABLES);

% Random variable generation
%VARIABLES_norm                       = unifrnd(0,1,1,nVar);
% VARIABLES_norm = zeros(1,nVar);

% Plot model 
hFig = fcnPLOT_GEOM(VARIABLES_norm, SETTINGS);

[OBJECTIVES, CONSTRAINTS] = fcnMAINSOLVER_TOM(VARIABLES_norm,0,SETTINGS,VAP_IN);

% Input generation
fixed = 1; free = 0;
[inputFilename, GEOMETRY, WEIGHT] = fcnGENERATE_INPUT(VARIABLES_norm, ID, SETTINGS, free);

% Testing stabilty function
if SETTINGS.AnalyseStability  
    dalpha = -10:4:10;
    STABILITY = fcnSTABILITY_ANALYSIS(inputFilename, GEOMETRY, SETTINGS, VAP_IN, dalpha); 
    STABILITY
    
    data = fcnSTAB_VALIDATION();
    
    % Stability plots 
    if SETTINGS.stabSelection(1)
        figure % Cma best cog around 0.041
            plot(dalpha,STABILITY.Cmcg_vec)
            grid on
            xlabel('\alpha [deg]')
            ylabel('C_M')    
            p = polyfit(dalpha,STABILITY.Cmcg_vec,1);
            hold on 
            plot(dalpha,polyval(p,dalpha))
            plot(data.Cma(:,1),data.Cma(:,2),'x')
            legend('VAP 3.5','VAP Fitted','Test')
    end
    if SETTINGS.stabSelection(2)
        figure % Clb best cog around 0.034
            plot(dalpha,STABILITY.Clcg_vec)
            grid on
            xlabel('\beta [deg]')
            ylabel('C_l')    
            p = polyfit(dalpha,STABILITY.Clcg_vec,1);
            hold on 
            plot(dalpha,polyval(p,dalpha))
            plot(data.Clb.beta_vec,data.Clb.alpha0,'x')
            legend('VAP 3.5','VAP Fitted','Test')    
    end
    if SETTINGS.stabSelection(3)
        figure % CNb - does not work if wing incidence is 4 deg
            plot(dalpha, STABILITY.CNcg_vec)
            grid on
            xlabel('\beta [deg]')
            ylabel('C_N')
            p = polyfit(dalpha,STABILITY.CNcg_vec,1);
            hold on 
            plot(dalpha, polyval(p,dalpha))
            plot(data.CNb.beta_vec, data.CNb.alpha0,'x')
            legend('VAP 3.5', 'VAP Fitted', 'Test')
    end
    
else 
    [OUTP, COND, INPU, FLAG, MISC, SURF, VEHI, VISC, WAKE] = fcnVAP_MAIN(inputFilename, VAP_IN);
end


% figure(3)
% hold on
% cog_x = GEOMETRY.cog_x
% plot3(-0.4+cog_x,0,0.005,'x','linewidth',4)













