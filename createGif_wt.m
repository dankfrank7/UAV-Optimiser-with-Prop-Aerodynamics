% This script generates .gifs of the 30 cm models with GWS4x4 propeller 
%
% Author: Tom Ryan

clear; clc
close all
InitialisePlotting()

% Workspace 
curr_dir = pwd;
cd './VAP3point5'
addpath(curr_dir)

loadFilename = 'Wingtips_GWS4x4_30span_DRIES_RESULTS';
load(strcat(curr_dir,'/OUTPUT/',loadFilename,'.mat'))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SETTINGS.maxtime = 80;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Normalised variables 
flyingWing_vars = [0.7013, 0.1354, 0.0607, 0.1260, 0.2173, 0.6369, 0.9783, 0.3, 0.2786, 0.3642, 0.8057];
wingTips_vars = [0.6096, 0.3484, 0.2835, 0.2881, 0.5933, 0.2903, 0.7288, 0.4266, 0.6065, 0.3268, 0.7, 0.1523];

VARIABLES_norm = wingTips_vars;            

% Plot model 
hFig = fcnPLOT_GEOM(VARIABLES_norm, SETTINGS);

% Input generation
fixed = 1; free = 0; ID = 0;
[inputFilename, GEOMETRY, WEIGHT] = fcnGENERATE_INPUT(VARIABLES_norm, ID, SETTINGS, free);

% Run Vap 
VAP_IN.PRINT = 0;
VAP_IN.PLOT  = 1;
VAP_IN.GIF   = 1;

[OUTP, COND, INPU, FLAG, MISC, SURF, VEHI, VISC, WAKE] = fcnVAP_MAIN(inputFilename, VAP_IN);