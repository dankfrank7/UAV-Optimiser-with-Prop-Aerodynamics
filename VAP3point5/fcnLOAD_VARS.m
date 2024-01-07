function fcnLOAD_VARS(VARIABLES, variableNames, varMin, varMax)
% Extracts the optimsation variables from the 1xnVar vector, denormailises 
% and saves them individually to the caller workspace using their names.  
%
% Inputs: VARIABLES, 1xnVar length vector of optimisation variables 
%
% Author: Tom Ryan
% 20/09/2021
    
    % Loop through each elemeng of the input variables 
    for i = 1:length(VARIABLES)
        assignin('caller',variableNames{i}, varMin(i) + VARIABLES(i)*(varMax(i)-varMin(i)))
    end