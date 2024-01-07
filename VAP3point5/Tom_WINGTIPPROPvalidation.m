% This function is used to create the Velduis PROWING validation study 
% 
% Author: Tom Ryan
% 1/09/2021

%% Start with a clear workspace
clear; clc;
close all
InitialisePlotting();

% Load in propeller geometry 
% r/R    c/R     beta (deg)
propellerName = 'wingtipprop';
propGeom = dlmread(strcat('propellers/',propellerName,'.dat'));
radialRaw = propGeom(:,1);
chordRaw = propGeom(:,2);
pitchRaw = propGeom(:,3);   
leadingEdgeRaw = fcnPROPLE(propellerName);
airfoilRaw = fcnPROPAIRFOIL(propellerName); % Pass propeller name for specific distribution OR specific airfoil for homogenous distribution

% Refine geom into usable variables 
noSectionsPerPanel = 6;
noPanels = length(airfoilRaw.af);
noSpanDivisions = noSectionsPerPanel*noPanels;  % Total number of divisions
pitchShift = 2; % +2 Deg looks good FOR APC PROPS, 4 deg for NACA5868-9

propDiam           = 0.236; %%%%
propRad            = propDiam./2;
radial_vec = linspace(radialRaw(1),0.98,noSpanDivisions); % Maybe change this to cosine spaced?
chord_vec  = interp1(radialRaw,chordRaw,radial_vec,'spline')*propRad; % Not normalised by radius 
pitch_vec  = interp1(radialRaw,pitchRaw,radial_vec,'spline') + pitchShift; 
leadingedge_vec = interp1(leadingEdgeRaw(:,1),leadingEdgeRaw(:,2),radial_vec,'spline')*propRad;

% Plotting 
plotProp = false; %%%%
if plotProp
    figure(1)
    plot(radial_vec, chord_vec,'-o')
    xlabel('r/R')
    ylabel('Chord (m)')
    grid minor
    title(propellerName)
    figure(2)
    plot(radial_vec, pitch_vec,'-o')
    xlabel('r/R')
    ylabel('Pitch (deg)')
    grid minor
    title(propellerName)
    figure(6)   
    axis equal 
    plot(radial_vec,leadingedge_vec./propRad,':')
    hold on 
    plot(radial_vec,(leadingedge_vec-chord_vec)./propRad,':k')
    drawnow;
end

%% Root XML file 
filename = 'Tom_WINGTIProot';           % Root file name
input = strcat(filename, '.VAP');       % Full input file with extension

% Convert input XML file to structure
inputS = fcnXML2STRUCT(input);

%% Input Settings
maxIter = 30;
noWingSpanEle = 6;
inputS.VAP.settings.maxtime = n2c(maxIter);
inputS.VAP.settings.delta_time = '0.0002';
inputS.VAP.settings.relax = 'TRUE';
inputS.VAP.vehicle.wing.chordwise_elements = '2';
inputS.VAP.vehicle.wing.panel.spanwise_elements = '2'; % One element per section bc we are doing half cosines spacing
inputS.VAP.vehicle.rotor.blades = '4';
inputS.VAP.vehicle.rotor.collective = '1';

%% Set Up Wing Geometry 
%Overrite root file with sections distributed to the end 
b_wing = 1.46;  % Full span
thetay = linspace(0,pi,noWingSpanEle);
yCos = b_wing/2*sin(thetay/2);
tempWingPanel = inputS.VAP.vehicle.wing.panel;
for i = 1:noWingSpanEle
    tempWingPanel.section{i}.wing_x = '0';
    tempWingPanel.section{i}.wing_y = n2c(yCos(i));
    tempWingPanel.section{i}.wing_z = '0';
    tempWingPanel.section{i}.chord = '0.24';
    tempWingPanel.section{i}.twist = '0';
end
inputS.VAP.vehicle.wing.panel = tempWingPanel;    

%% Set Up Propeller geometry 
% Loop through each propeller section 
r = 1;                          % r goes from 1:noSpanDivisions
tempPanel = cell(1,noPanels);
% Panel
for i = 1:noPanels 
    
    
    tempPanel{i}.strip_airfoil = airfoilRaw.af{i};
    tempPanel{i}.spanwise_elements = '1';
    
    % Section
    for j = 1:noSectionsPerPanel
        tempPanel{i}.section{j}.rotor_x = n2c(cosd(pitch_vec(r))*-leadingedge_vec(r));  %n2c(-chord_vec(i)/2);  % Global -Z
        tempPanel{i}.section{j}.rotor_y = n2c(radial_vec(r)*propDiam/2);                % Global Y 
        tempPanel{i}.section{j}.rotor_z = n2c(sind(pitch_vec(i))*chord_vec(i)/2);       % Global -X 
        tempPanel{i}.section{j}.chord   = n2c(chord_vec(r));
        tempPanel{i}.section{j}.twist   = n2c(pitch_vec(r));
        
        r = r + 1;  % Next division 
    end
end
inputS.VAP.vehicle.rotor.panel = tempPanel;
holdRotor = inputS.VAP.vehicle.rotor;

%% Run Simulation 
alpha_vec = -5:5:15;
J_vec = [0, 0.9, 0.8];
Vel = 40;               % From paper
RPM_range = 60*Vel./J_vec./propDiam;      % RPM (rev/min)

Cl_vec = nan(length(J_vec),length(alpha_vec));
Cd_vec = nan(length(J_vec),length(alpha_vec));
CT_vec = nan(length(J_vec),length(alpha_vec));
CP_vec = nan(length(J_vec),length(alpha_vec));

% Loop through modes (2 advance ratio + 1 no prop) 
for i = 1:length(J_vec)
    
    if J_vec(i) == 0        
        fprintf('No Rotor\n')
        
        inputS.VAP.vehicle.rotor = {};  % Remove rotor if advance ratio is 0
    else
        fprintf('RPM %.0f\n',RPM_range(i))
        fprintf('Advance Ratio %.4f\n',J_vec(i))
        
        inputS.VAP.vehicle.rotor = holdRotor; % Revert rotor 
        
        % Update advance ratio 
        inputS.VAP.vehicle.rotor.rpm    = RPM_range(i); % In RPM
    end

    for j = 1:length(alpha_vec)
        
        fprintf('Alpha %.4f\n',alpha_vec(j))
        
        % Update angle of attack
        inputS.VAP.vehicle.wing.incidence = n2c(alpha_vec(j));
   
        % Convert structure back to XML file
        newfilename = 'tempVAPINPUT';              
        newinput    = strcat(newfilename, '.VAP'); 
        fcnSTRUCT2XML(inputS,newinput);   

        VAP_IN = [];
        [OUTP, ~, ~, ~, ~, ~, ~, ~, ~]  = fcnVAP_MAIN(newinput, VAP_IN);
        
        if J_vec(i) == 0 
            Cl_vec(i,j) = OUTP.vecCLv_AVG;
            Cd_vec(i,j) = OUTP.vecCD_AVG;
            CT_vec(i,j) = 0;
            CP_vec(i,j) = 0;
        else
            Cl_vec(i,j) = OUTP.vecCL_AVG;
            Cd_vec(i,j) = OUTP.vecCD_AVG;
            CT_vec(i,j) = mean(OUTP.vecCT(end-10:end));
            CP_vec(i,j) = mean(OUTP.vecCP(end-10:end));
        end
        clc
    end
end
Cx_vec = Cd_vec - CT_vec;

%% Load Validation Data 
% Select dataset
valData = 'PropWingVal2';
val = fcnINITPROPVAL(valData);
fn = fieldnames(val);

for k = 2:length(fn)
    FullValJ_vec(k) = str2num(fn{k}(end-1:end))/10;
end

for kk = 1:length(J_vec) 
    if isempty(find(J_vec == FullValJ_vec(kk)))
        val = rmfield(val,fn{kk});
    end
end
fnNEW = fieldnames(val);    % Update fieldnames after removing non-relavant ones 

figure
hold on
for iii = 1:length(J_vec)
    plot(Cx_vec(iii,:),Cl_vec(iii,:),'-o')
    tempVal = val.(fnNEW{iii});
    plot(tempVal(:,1), tempVal(:,2), '-x')
end
grid minor 
xlabel('C_D')
ylabel('C_L')
legend('VAP 3.5: Prop OFF','Test: Prop OFF','VAP 3.5: J = 0.9','Test: J = 0.9','VAP 3.5: J = 0.8','Test: J = 0.8')

%% Subroutines 
% Numeric to character vector
function output = n2c(double)
    output = char(string(double));
end



