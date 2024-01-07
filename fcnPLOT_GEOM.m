function hFig = fcnPLOT_GEOM(VARIABLES,SETTINGS)
% This function plots the aircraft design based on VARIABLES hopefully
% without going into the main VAP 3.5 file.
%
% Author: Tom Ryan

    % Generate .vap input file
    [inputFilename, ~, ~] = fcnGENERATE_INPUT(VARIABLES, 1, SETTINGS, 0);

    % Read in the Geometry
    VAP_IN = [];
    [~, COND, VISC, INPU, VEHI] = fcnXMLREAD(inputFilename, VAP_IN);
    
    % Initializing parameters to null/zero/nan
    [WAKE, OUTP, INPU, SURF] = fcnINITIALIZE(COND, INPU);
    
    % Discretizing geometry into DVE
    [INPU, ~, ~, ~, ~, ~, SURF, ~] = fcnGEOM2DVE(INPU, COND, VISC, VEHI, WAKE, OUTP, SURF);
    
    if ~isempty(INPU.vecSYM) && any(INPU.vecSYM)
        sym = true;
    else
        sym = false;
    end

    % Body plot
    verbose = 0; % pLots panel indexing 
    hFig = fcnPLOTBODY_2(verbose, SURF.valNELE, SURF.matDVE, SURF.matVLST, SURF.matCENTER, sym);
end