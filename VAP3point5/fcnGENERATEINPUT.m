function [inputFilename, GEOMETRY] = fcnGENERATE_INPUT(VARIABLES, ID, SETTINGS)
% This function takes a geometry defiend by VARIABLES and SETTINGS and
% generates the .VAP (xml) input file to be called in fcnMAINSOLVER_TOM()
% 
% INPUTS: 
%   VARIABLES -  variables that are varied by the optimisation algorithm  
%   ID - not sure where this fits in yet 
%   SETTINGS - constaint geometry/VAP settings that the opt alg does not
%   touch 
%
% OUTPUT: 
%   *.vap* saves inputFilename.vap to workspace
%   inputFilename - string INCLUDING .vap 
%   GEOMETRY - structure defining geometry of current iteration
%
% THIS IS OLD VERSION DO NOT USE 22/09/2021
    
    % Extract variables
    fn_var = fieldnames(VARIABLES);
    for i = 1:numel(fn_var)
        assignin('caller',fn_var{i},VARIABLES.(fn_var{i}));
    end
    fn_set = fieldnames(SETTINGS);
    for i = 1:numel(fn_set)
        assignin('caller', fn_set{i}, SETINGS.(fn_set{i}));
    end
    
    % Load root file 
    rootFilename = 'Tom_mainROOT';          
    input = strcat(rootFilename, '.VAP'); 
    
    % Convert input XML file to structure
    inputS = fcnXML2STRUCT(input);

    %% SETTINGS 
    tempSettings = inputS.VAP.settings;
    tempSettings.maxtime    = n2c(maxtime);
    tempSettings.delta_time = n2c(delta_time);
    tempSettings.relax      = n2c(relax);
    tempSettings.steady     = n2c(steady);
    inputS.settings = tempSettings;
    
    %% WING 
    tempWing = inputS.vehicle.wing;
    tempWing.symmetry   = n2c(symmetry);
    tempWing.incidence  = n2c(incidence);
    tempWing.chordwise_elements = n2c(chordwise_elements);
    
    switch wingSpanSpacing
        case 'linear'
            y_space     = linspace(0,1,noWingSpanEle); 

        case 'halfsine'
            thetay      = linspace(0,pi,noWingSpanEle);
            y_space     = sin(thetay/2);

        otherwise 
            error('Spanwise spacing type not recognised - TR')
    end
    
    switch planformType 
        case 'rectangular'
            wing_y       = ref_span/2*y_space;
            wing_x       = wing_y.*tand(sweep);
            wing_z       = wing_y.*tand(dihedral);
            twist_vec    = twist*y_space;
            chord_vec    = c_root*(1-taper*y_space);
            ref_cmac = c_root*2/3*((1+taper+taper^2)/(1+taper));
            inputS.VAP.vehicle.ref_cmac = ref_cmac;
            
        otherwise
            error('Planform shape not recognised - TR')
    end 
   
    for i = 1:noWingSpanEle 
        tempWing.panel.section{i}.wing_x = n2c(wing_x(i));
        tempWing.panel.section{i}.wing_y = n2c(wing_y(i));
        tempWing.panel.section{i}.wing_z = n2c(wing_z(i));
        tempWing.panel.section{i}.twist  = n2c(twist_vec(i));
        tempWing.panel.section{i}.chord  = n2c(chord_vec(i));
    end
    inputS.vehicle.wing = tempWing;
    
    ref_area = fcnsWINGAREA(wing_y,chord_vec);
    inputS.VAP.vehicle.ref_area = ref_area;
    inputS.VAP.vehicle.ref_span = ref_span;
    
    %% WEIGHT 
    % Weigt calculated as per paper from Dries, alpha used to set L = W 
    inputS.VAP.vehicle.alpha = n2c(alpha);
   

    %% PROPELLER 
    if strcat(propellerName,'noProp')
        inputS.vehicle.propeller = {};
        
    else
        propRad  = propdDiam/2;
        propGeom = dlmread(strcat('propellers/',propellerName,'.dat'));
        radialRaw = propGeom(:,1);
        chordRaw  = propGeom(:,2);
        pitchRaw  = propGeom(:,3);
        leadingEdgeRaw = fcnPROPLE(propellerName);
        
        radial_vec      = linspace(radialRaw(1) ,0.98 ,noRotorSpanEle);
        chordr_vec      = interp1(radialRaw,chordRaw,radial_vec,'spline')*propRad;
        pitch_vec       = interp1(radialRaw,pitchRaw,radial_vec,'spline') + pitch_shift;
        leadingedge_vec = interp1(leadingEdgeRaw(:,1),leadingEdgeRaw(:,2),radial_vec,'spline')*propRad;
        
        tempRotor = inputS.VAP.vehicle.rotor;
        tempRotor.ref_diam           = n2c(propDiam);
        tempRotor.rotation_direction = n2c(rotation_direction);
        tempRotor.blades             = n2c(blades);
        
        tempPanel = cell(1,noPanels);
        tempPanel.strip_airfoil = n2c(propAirfoil);
        for i = 1:noRotorSpanEle 
            tempPanel.section{i}.rotor_x = n2c(-leadingedge_vec(i)*cosd(pitch_vec(i))); % Global -Z 
            tempPanel.section{i}.rotor_y = n2c(radial_vec(i)*propRad);                  % Global +Y
            tempPanel.section{i}.rotor_z = n2c(sind(pitch_vec(i))*chordr_vec(i)/2);      % Global -X 
            tempPanel.section{i}.chord   = n2c(chordr_vec(i));
            tempPanel.section{i}.twist   = n2c(pitch_vec(i));
        end
        
        tempRotor.panel = tempPanel;
        inputS.VAP.vehicle.rotor = tempRotor; 
    end
        
    %% OUTPUT
    name = 'tempVAPINPUT';
    inputFilename = strcat(name,'.vap');
    fcnSTRUCT2XML(inputS, inputFilename)
    
    GEOMETRY.span      = ref_span;
    GEOMETRY.wing_area = ref_area;
    GEOMETRY.cmac      = ref_cmac; 
end

%% Subroutines
function ref_area = fcnsWINGAREA(wing_y,chord_vec)
% Calculates the wing area based on the sum of trapezedoil sections 

    panelArea = zeroes(1,length(wing_y)-1);

    for i = 1:length(wing_y)-1
        dy = wing_y(i+1)-wing(i);
        panelArea(i) = 0.5*dy*(chord_vec(i+1) + chord_vec(i));
    end
      
    ref_area = sum(panelArea);
end