function [OBJECTIVES, CONSTRAINTS] = fcnMAINSOLVER_TOM(VARIABLES, ID, SETTINGS, VAP_IN)
%% This function calls VAP 3.5 or g a given geometry and solver settings and
% extracts the OBJECTIVES and CONSTRAINTS for use in MainOptimisation Loop
% 
% Author: Tom Ryan
% 6/09/2021

    %% Generate inpute .VAP file based on geometry and settings   
    free  = 0;
    fixed = 1;
    [inputFilename, GEOMETRY, WEIGHT] = fcnGENERATE_INPUT(VARIABLES, ID, SETTINGS, free);
    
    %% Call VAP 3.5 
    [OUTP, COND, ~, ~, ~, SURF, ~, ~, ~] = fcnVAP_MAIN(inputFilename, VAP_IN);    
    
    % Numerically Calculate stability derivatives 
    if SETTINGS.AnalyseStability   
        STABILITY = fcnSTABILITY(inputFilename, GEOMETRY, SETTINGS, VAP_IN, OUTP, SURF, COND, SETTINGS.maxTimeStab, ID); 
    end       
    
    if ~strcmp(SETTINGS.propellerName, 'noProp')
        Ct = OUTP.vecCT_AVG;
        Thrust = SETTINGS.density*(SETTINGS.rpm/60)^2*SETTINGS.propDiam^4*Ct;
    end
    
    %% Extract Objectives 
    Cl = OUTP.vecCLv_AVG;
    Cd = OUTP.vecCD_AVG;
    Lift = 0.5*SETTINGS.density*SETTINGS.speed.^2*GEOMETRY.ref_area*Cl;
    Drag = 0.5*SETTINGS.density*SETTINGS.speed.^2*GEOMETRY.ref_area*(Cd + SETTINGS.Cd_extra);
    
    OBJECTIVES = [-Cl./Cd, WEIGHT];
    
    % CONSTRAINTS 
    % Geom constraints
    spanConstraint = GEOMETRY.ref_span - SETTINGS.Wingspan_constraint;
    % Length Vectors 
    xMax_vt = 0; 
    xMax_ht = 0; 
    xMax_wing = max(GEOMETRY.wing_x + GEOMETRY.chord_vec);
    if ~strcmp(SETTINGS.verticalTailType,'noVerticalTail')
        xMax_vt = max(GEOMETRY.vert_x + GEOMETRY.vertChord_vec);
    end
    if ~strcmp(SETTINGS.horizontalTailType,'noHorizontalTail')
        xMax_ht = max(GEOMETRY.hori_x + GEOMETRY.horChord_vec);
    end
    lengthConstraint = max([xMax_ht, xMax_vt, xMax_wing]) - SETTINGS.Length_constraint;
    
    
    % Trim Constraint
    if strcmp(SETTINGS.propellerName,'noProp')
        trimConstraint = WEIGHT*9.81 - Lift; % NEED TO INCLUDE TOLERANCE
        
        rpmConstraint = [];
    else 
        trimConstraint = [WEIGHT*9.81 - Lift, Drag - Thrust]; % NEED TO INCLUDE TOLERANCE
        
        rpmConstraint = SETTINGS.rpm - SETTINGS.rpmMax;
    end
    
    % Stabiltiy Constraints
    if SETTINGS.AnalyseStability
        CmaConstraint = [STABILITY.Cma - SETTINGS.Cma_upper, SETTINGS.Cma_lower - STABILITY.Cma];
        ClbConstraint = [STABILITY.Clb - SETTINGS.Clb_upper, SETTINGS.Clb_lower - STABILITY.Clb];
        CNbConstraint = [STABILITY.CNb - SETTINGS.CNb_upper, SETTINGS.CNb_lower - STABILITY.CNb];
        CNb_Clb = abs(STABILITY.CNb / STABILITY.Clb);
        DutchRollConstraint = [CNb_Clb - 2/3, 1/3 - CNb_Clb]; 
        stabConstraint = [CmaConstraint, ClbConstraint, CNbConstraint, DutchRollConstraint];

    else 
        stabConstraint = [];
    end
    
    CONSTRAINTS = [spanConstraint, lengthConstraint, trimConstraint, rpmConstraint, stabConstraint];

    
end
    