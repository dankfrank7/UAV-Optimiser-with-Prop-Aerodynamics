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
propellerName = 'naca5868_9_pitch25';
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
pitchShift = 4; % +2 Deg looks good FOR APC PROPS, 4 deg for NACA5868-9

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
filename = 'Tom_PROWIMroot';           % Root file name
input = strcat(filename, '.VAP');       % Full input file with extension

% Convert input XML file to structure
inputS = fcnXML2STRUCT(input);

%% Input Settings
maxIter = 140;
noWingSpanEle = 16;
inputS.VAP.settings.maxtime = n2c(maxIter);
inputS.VAP.settings.delta_time = '0.0002';
inputS.VAP.settings.relax = 'TRUE';
inputS.VAP.vehicle.alpha = '4';
inputS.VAP.vehicle.wing.chordwise_elements = '2';
inptuS.VAP.vehicle.wing.spanwise_elements = n2c(noWingSpanEle);
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

%% Run Simulation 
J_vec = 0.8:0.1:1.2;       % Taken from experimental data range
J = 0.9;                 % For testing
Vel = 49.25;               % From paper
RPM = 60*Vel/J/propDiam;      % RPM (rev/min)
inputS.VAP.vehicle.rotor.rpm    = RPM;  % In RPM

Cn_vec = nan(noWingSpanEle,3);
% 1) CCW == InUp
% 2) CW  == OutUp
% 3) PropOff
alpha_vec = [0, 4, 8];
mode = {'InUp','OutUp','PropOff'};
rotation = {'CCW','CW'};
valData = {};
noRuns = 3;
for i = 1:noRuns % Change this to 3 for the 4deg AoA
    fprintf('Mode: %s\n ', mode{i})
    fprintf('RPM %.0f\n',RPM)
    fprintf('Advance Ratio %.4f\n',J)

%     if i < 3
%         tempRotor = inputS.VAP.vehicle.rotor;
%         tempRotor.rotation_direction = rotation{i};
%     else
%         tempRotor = {}; % Delete propeller if in Prop Off mode 
%     end
    tempRotor = inputS.VAP.vehicle.rotor;
    tempRotor.rotation_direction = rotation{1};
    
    inputS.VAP.vehicle.alpha = n2c(alpha_vec(i));
    
    % Stpre it back in inpout 
    inputS.VAP.vehicle.rotor = tempRotor;
    
    % Convert structure back to XML file
    newfilename = 'tempVAPINPUT';              
    newinput    = strcat(newfilename, '.VAP'); 
    fcnSTRUCT2XML(inputS,newinput);   

    VAP_IN.PRINT = 1;
    VAP_IN.PLOT = 1;
    VAP_IN.GIF = 1;
%     VAP_IN.CIRCPLOT = 1;
    [OUTP, ~, ~, ~, ~, ~, ~, ~, ~]  = fcnVAP_MAIN(newinput, VAP_IN);

    Cn_vec(:,i) = OUTP.vecCNDIST;
    valData{i} = strcat('PROWIMvalidation',inputS.VAP.vehicle.alpha);
end

%% Load Validation Data 
% Select dataset


%fn = fieldnames(val.Test);

yB_vec = linspace(0,1,noWingSpanEle);
figure
hold on 
% VAP 3.5 Data
for i = 1:noRuns 
    plot(yB_vec,Cn_vec(:,i),'-')
end
xlabel('2y/b')
ylabel('C_N')
grid minor 
% Test and VLM Data 
% for i = 1:length(fn)
%     temp_test = val.Test.(fn{i});
%     plot(temp_test(:,1),temp_test(:,2),'x')
%     %temp_vlm = val.VLM.(fn{i});
%     %plot(temp_vlm(:,1),temp_vlm(:,2),'-')
% end
myred           = [216 30 49]/255;
myblue          = [27 99 157]/255;
myblack         = [0 0 0]/255;
colorVec = {myblack, myblue, myred};
for i = 1:length(alpha_vec)
    val = fcnINITPROPVAL(valData{i});
    temp_test = val.Test.InUp;
    plot(temp_test(:,1),temp_test(:,2),'x','color',colorVec{i})
end

lgd = legend('VAP 3.5: 0 deg \alpha','VAP 3.5: 4 deg \alpha','VAP 3.5: 8 deg \alpha', 'Test: 0 deg \alpha', 'Test: 4 deg \alpha', 'Test: 8 deg \alpha');
%lgd = legend('VAP 3.5: InUp','VAP 3.5: OutUp','VAP 3.5: No Prop','Test: No Prop','VLM: No Prop', 'Test: InUp','VLM: InUP', 'Test: OutUp','VLM: OutUP');
%lgd = legend('VAP 3.5: InUp','Test: InUp', 'VLM: InUp');
%lgd.Title.String = strcat('Alpha: ',inputS.VAP.vehicle.alpha,' deg');

%% Subroutines 
% Numeric to character vector
function output = n2c(double)
    output = char(string(double));
end



