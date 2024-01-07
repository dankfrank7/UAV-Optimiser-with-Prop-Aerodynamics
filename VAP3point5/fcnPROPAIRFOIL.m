function airfoilData = fcnPROPAIRFOIL(propellerName)
% This function outputs the NACA 4-series airfoil distrbution across the 
% propeller span as described by MacNeill, 2017. 
% Output is in structure: S.rR = [0.2, 0.3 .... 1], S.af =
% {'nacaXXXX',....}
%
%
% Author: Tom Ryan
% 24/08/2021
    switch propellerName 
        case 'wingtipprop'
            airfoilData.rR = 1;
            airfoilData.af = {'clarky'};
        case 'naca5868_9_pitch15'
            airfoilData.rR = 1;
            airfoilData.af = {'clarky'};
        case 'naca5868_9_pitch25'
            airfoilData.rR = 1;
            airfoilData.af = {'clarky'};
        case 'APC19x12'
            airfoilData.rR = 0.2:0.1:1;
            airfoilData.af = {'naca5524','naca4418','naca4515','naca4513','naca4512','naca4412','naca4410','naca4309','naca4309'};
        case 'APC17x12'
            airfoilData.rR = 0.2:0.1:1;
            airfoilData.af = {'naca4521','naca5515','naca5514','naca5513','naca4412','naca4411','naca4410','naca4409','naca4409'};
        otherwise
            airfoilData.rR = 1;
            airfoilData.af = {propellerName}; 
    end
return