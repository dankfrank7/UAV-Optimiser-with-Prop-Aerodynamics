% VAP 3.5 Testing 
% Author: Thomas Ryan
% 3/08/2021
%
% This script generates plots for the Trips validation data

% Start with a clear workspace
clear; clc
close all
InitialisePlotting()

% Root XML file 
filename = 'Tom_TripsROOT';                % Root file name
input = strcat(filename, '.VAP');       % Full input file with extension

%% Varying Input

% Set up vectors
Cd_correction = 0.005; % This has now been shifted into the fcnVISCOUS_WING.m file
alpha_vec   = 0:1:10;
Cl_vec      = zeros(size(alpha_vec));
Cd_vec      = zeros(size(alpha_vec));

% Geometry variables for the Trips validation
root_chord  = 0.17;
taper_vec   = [0.3,0.5,1];
tip_chord   = root_chord*taper_vec;
ref_area    = [0.0366, 0.0487, 0.0867];
ref_span    = [0.332, 0.383, 0.51];
ref_cmac    = [0.121, 0.132, 0.17];
tip_y       = ref_span/2;

% Convert input XML file to structure
inputS = fcnXML2STRUCT(input);

% Change Input File Settings 
inputS.VAP.vehicle.wing.panel.spanwise_elements = '8';
inputS.VAP.settings.maxtime                     = '20';
inputS.VAP.vechicle.wing.chordwise_elements     = '1';
inputS.VAP.settings.relax                       = 'TRUE';
        
% MODE 1 Trips 0.3 Taper
% MODE 2 Trips 0.5 Taper
% MODE 3 Trips 1.0 Taper

for mode = 1:3
    
    % Change geometry 
    inputS.VAP.vehicle.ref_area                     = fcnN2C(ref_area(mode));
    inputS.VAP.vehicle.ref_span                     = fcnN2C(ref_span(mode));
    inputS.VAP.vehicle.ref_cmac                     = fcnN2C(ref_cmac(mode));
    inputS.VAP.vehicle.wing.panel.section{2}.chord  = fcnN2C(tip_chord(mode));
    inputS.VAP.vehicle.wing.panel.section{2}.wing_y.Text = fcnN2C(tip_y(mode));  

    % VAP_IN = [];
    % [OUTP, ~, ~, ~, ~, ~, ~, ~, ~] = fcnVAP_MAIN('TomTesting_2.VAP', VAP_IN);

    for i = 1:length(alpha_vec)

        % Change variables in structure
        fprintf('Mode %.2f\n',mode)
        fprintf('Angle of attack %.2f deg\n',alpha_vec(i))
        inputS.VAP.vehicle.alpha        = fcnN2C(alpha_vec(i)); % Inputs need to be character array 

        % Convert structure back to XML file
        newfilename = 'tempVAPINPUT';              % .VAP will be added to end if not already in the filename
        newinput    = strcat(newfilename, '.VAP'); % MUST BE IN CAPITALS
        fcnSTRUCT2XML(inputS,newinput);            % Will save in the local folder, overwritten each time 

        % Call the main function
        VAP_IN.PRINT = 1;
        VAP_IN.PLOT = 1;
        VAP_IN.GIF = 1;
        [OUTP, ~, ~, ~, ~, ~, ~, ~, ~]  = fcnVAP_MAIN(newinput, VAP_IN);
        Cl_vec(i,mode)  = OUTP.vecCLv_AVG;
        Cd_vec(i,mode)  = OUTP.vecCD_AVG + Cd_correction;
        clc 
    end

    % Lift to drag 
    LtoD_vec = Cl_vec./Cd_vec;

    % Trips Validation Data
    data = fcnINITTRIPS();
    
    % Plotting 
    leg_vec = {'Taper 0.3','Taper 0.5','Taper 1.0'};
    
    figure % Drag polar
    plot(Cd_vec(:,mode),Cl_vec(:,mode),'-o')
    grid minor
    xlabel('C_D')
    ylabel('C_L')
    hold on 
    plot(data.Trips{mode}.CD,data.Trips{mode}.CL,'-o')
    lgd1 = legend('VAP 3.5','CFD','location','best');
    lgd1.Title.String = leg_vec{mode};

    figure % Lift to drag vs Cl
    plot(Cl_vec(:,mode),LtoD_vec(:,mode),'-o')
    grid minor 
    xlabel('C_L')
    ylabel('L/D')
    LD_ref = data.Trips{mode}.CL./data.Trips{mode}.CD;
    hold on 
    plot(data.Trips{mode}.CL, LD_ref,'-o')
    lgd2 = legend('VAP 3.5','CFD','location','best');
    lgd2.Title.String = leg_vec{mode};
end
%% QUESTIONS 4/08/2021

% HOW TO SET FLAG, WHAT IS THE PURPOSE OF 'VAP_IN' IF IT IS IMMEDIATELY
% OVERWRITEN

% WHAT IS THE CIRC PLOT TOGGLE PLOTTING? IS IT LIFT DISTRIBUTION

% WHAT IS THE INFLUENCE OF MULTIPLE DVE PANNELS? DOES IT MATTER? IS IT AS
% IMPORTANT AS IN A VLM/PANEL CODE?

% IS THERE A DIFFERENCE BETWEEN MANTUALLY CREATING SECTIONS AND DEFINING
% ROOT AND TIP + SPANWISE SECTIONS??

% HOW DOES THE THRUST COEFFICIENT INFLUENCE THE DRAG COEFFICIENT, IGNORING
% NET X FORCE?

% IS THE VALIDATION ERROR COMING FROM THE VISCOUS PROPELLER DATA? IS IT AT
% A DIFFERENT REYNOLDS NUMBER

%% NEXT STEPS

% ADD PROPELLER TO INPUT FILE 

% ADD AN ANGLE OF ATTACK SWEEP - Complete 4/08/2021

% ADD GEOMETRY SWEEP - Complete 5/08/2021

%% Subroutines 
% Numeric to character vector
function output = fcnN2C(double)
    output = char(string(double));
end