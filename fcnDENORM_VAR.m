function [VAR_denorm, label] = fcnDENORM_VAR(VAR_norm, variableNames, varMin, varMax, i)
% This function denormalises the variable vector VAR_vec 
% Variable is defined by index i corresponding to name in variableNames and
% bounds defined in varMin, varMax 
%
% Author: Tom Ryan

    % Library for the color bar label
    labelLib.ref_span = 'Span [m]';
    labelLib.c_root   = 'Root Chord [m]';
    labelLib.sweep    = 'Sweep [deg]';
    labelLib.taper    = 'Taper Ratio';
    labelLib.dihedral = 'Dihedral [deg]';
    labelLib.alfa     = 'Alpha [deg]';
    labelLib.rpm      = 'RPM';
    labelLib.AR       = 'Aspect Ratio';
    labelLib.twist    = 'Twist [deg]';
    labelLib.volumeTail = 'V_{vt}';
    labelLib.l_vt      = 'VT Location [c_{root}]';
    labelLib.ar_vt    = 'VT Aspect Ratio';
    labelLib.ar_ht    = 'HR Aspect Ratio';
    labelLib.volumeHtail = 'V_{ht}';
    labelLib.sweep_vt = 'Sweet VT [deg]';

    try 
        label = labelLib.(variableNames{i});
    catch 
        label = variableNames{i};
    end
    
    % Denormalise variables 
    VAR_denorm = varMin(i) + VAR_norm.*(varMax(i)-varMin(i));
end
    
    