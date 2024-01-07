function [inputFilename, GEOMETRY, WEIGHT] = fcnGENERATE_INPUT(VARIABLES, ID, SETTINGS, fixedLift)
% This function takes a geometry defiend by VARIABLES and SETTINGS and
% generates the .VAP (xml) input file to be called in fcnMAINSOLVER_TOM()
% 
% INPUTS: 
%   VARIABLES -  variables that are varied by the optimisation algorithm (1xnVar vector) 
%                these come in normalised form and are immediately
%                denormailised
%   ID        - not sure where this fits in yet 
%   SETTINGS  - constaint geometry/VAP settings that the opt alg does not
%               touch (struct)
%   fixedLift - Boolean, TRUE if fixedLift analysis, FALSE if velocity
%               defined
%
% OUTPUT: 
%   *.vap* saves inputFilename.vap to workspace
%   inputFilename - string INCLUDING .vap 
%   GEOMETRY - structure defining geometry of current iteration
% 
% Author: Tom Ryan 
% 7/09/2021
%
% Status: THIS IS THE NEW VERSION 25/10/2021
    
    % Extract variables
    fcnSTRUCT2VARS(SETTINGS)
    
    % Load denormailised variables and save them as their name 
    fcnLOAD_VARS(VARIABLES, variableNames, varMin, varMax)
    
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
    inputS.VAP.settings = tempSettings;    
    
    if fixedLift 
        inputS.VAP.vehicle.speed = 'nan';   % Triggers fixed lift analysis
    else
        inputS.VAP.vehicle.speed = n2c(speed);
    end
    inputS.VAP.vehicle.alpha = n2c(alfa);
    inputS.VAP.vehicle.beta  = n2c(beeta);
    
    %% WING 
    tempWing = inputS.VAP.vehicle.wing;
    tempWing.symmetry            = n2c(symmetry);
    tempWing.incidence           = n2c(incidence);
    tempWing.chordwise_elements  = n2c(chordwise_elements);
    tempWing.panel.strip_airfoil = n2c(wingAirfoil);
    
    switch wingSpanSpacing
        case 'linear'
            y_space     = linspace(0,1,noWingSpanEle); 

        case 'halfsine'
            thetay      = linspace(0,pi,noWingSpanEle);
            y_space     = sin(thetay/2);

        otherwise 
            error('Wing Spanwise spacing type not recognised - TR')
    end
    
    switch planformType 
        case 'rectangular'
            wing_y       = ref_span/2*y_space;
            wing_x       = wing_y.*tand(sweep);
            wing_z       = wing_y.*tand(dihedral);
            twist_vec    = twist*y_space;
            chord_vec    = c_root + c_root*y_space*(taper-1);
            ref_cmac     = c_root*2/3*((1+taper+taper^2)/(1+taper));
            y_cmac         = ref_span/6*(1+2*taper)/(1+taper);
            x_ac         = y_cmac*tand(sweep)+0.25*ref_cmac;
            inputS.VAP.vehicle.ref_cmac = ref_cmac;
            
        case 'zimmerman'
            b_min = c_root/4;
            b_maj = 3*b_min;
            a = c_root*AR/2;
            
            wing_y      = y_space*a;
            for m = 1:noWingSpanEle
                x_min(m) = sqrt(b_min^2*(1-wing_y(m)^2/a^2));
                x_maj(m) = -sqrt(b_maj^2*(1-wing_y(m)^2/a^2));
            end
            chord_vec = x_min-x_maj;
            wing_x      = b_min - chord_vec/4 + wing_y.*tand(sweep);
            wing_z      = wing_y.*tand(dihedral);
            twist_vec   = twist*y_space;
            
            ref_span = a*2;
            
            % Mean Aerodynamic Chord
            a_def = a;
            b_def = 0;
            syms x_int
            fun1 = sqrt(b_min^2*(1-x_int^2/a^2));
            fun2 = -sqrt(b_maj^2*(1-x_int^2/a^2));
            c = fun1-fun2;
            c2 = c^2;
            ref_cmac = real(double((2*int(c2,b_def,a_def))/(2*int(c,b_def,a_def)))); 
            
            % Interpolate to find aerodynamic centre
            ind_cmac = interp1(chord_vec, 1:noWingSpanEle, ref_cmac,'linear');
            x_ac     = interp1(1:noWingSpanEle,wing_x,ind_cmac,'linear')+0.25*ref_cmac;
            y_cmac   = interp1(1:noWingSpanEle,wing_y,ind_cmac,'linear');
            
        case 'invZimmerman'
            b_min = c_root/4;
            b_maj = 3*b_min;
            a = c_root*AR/2;
            
            wing_y      = y_space*a;
            for m = 1:noWingSpanEle
                x_min(m) = -sqrt(b_min^2*(1-wing_y(m)^2/a^2));
                x_maj(m) = sqrt(b_maj^2*(1-wing_y(m)^2/a^2));
            end
            chord_vec = x_maj-x_min;
            wing_x      = b_min - chord_vec/4 + wing_y.*tand(sweep);
            wing_z      = wing_y.*tand(dihedral);
            twist_vec   = twist*y_space;
            
            ref_span = a*2;
            
            % Mean Aerodynamic Chord
            a_def = a;
            b_def = 0;
            syms x_int
            fun1 = sqrt(b_min^2*(1-x_int^2/a^2));
            fun2 = -sqrt(b_maj^2*(1-x_int^2/a^2));
            c = fun1-fun2;
            c2 = c^2;
            ref_cmac = real(double((2*int(c2,b_def,a_def))/(2*int(c,b_def,a_def)))); 
            
            % Interpolate to find aerodynamic centre
            ind_cmac = interp1(chord_vec, 1:noWingSpanEle, ref_cmac,'linear');
            x_ac     = interp1(1:noWingSpanEle,wing_x,ind_cmac,'linear')+0.25*ref_cmac;
            y_cmac   = interp1(1:noWingSpanEle,wing_y,ind_cmac,'linear');
        
        case 'elliptical'
            b = c_root/2;
            a = c_root*AR/2;
            
            wing_y      = y_space*a;
            for m = 1:noWingSpanEle
                x_min(m) = sqrt(b^2*(1-wing_y(m)^2/a^2));
                x_maj(m) = -sqrt(b^2*(1-wing_y(m)^2/a^2));
            end
            chord_vec = x_min-x_maj;
            wing_x      = b - chord_vec/2 + wing_y.*tand(sweep);
            wing_z      = wing_y.*tand(dihedral);
            twist_vec   = twist*y_space;
            
            ref_span = a*2;
            
            % Mean Aerodynamic Chord
            a_def = a;
            b_def = 0;
            syms x_int
            fun1 = sqrt(b^2*(1-x_int^2/a^2));
            fun2 = -sqrt(b^2*(1-x_int^2/a^2));
            c = fun1-fun2;
            c2 = c^2;
            ref_cmac = real(double((2*int(c2,b_def,a_def))/(2*int(c,b_def,a_def)))); 
            
            % Interpolate to find aerodynamic centre
            ind_cmac = interp1(chord_vec, 1:noWingSpanEle, ref_cmac,'linear');
            x_ac     = interp1(1:noWingSpanEle,wing_x,ind_cmac,'linear')+0.25*ref_cmac;
            y_cmac   = interp1(1:noWingSpanEle,wing_y,ind_cmac,'linear');
            
            
        otherwise
            error('Planform shape not recognised - TR')
    end 
   
    for side = 0:1
        RHS = 1;
        if side
            RHS = -1;
        end
        for i = 1:noWingSpanEle 
            ii = side*noWingSpanEle + i;
            tempWing.panel.section{ii}.wing_x = n2c(wing_x(i));
            tempWing.panel.section{ii}.wing_y = n2c(RHS*wing_y(i));
            tempWing.panel.section{ii}.wing_z = n2c(wing_z(i));
            tempWing.panel.section{ii}.twist  = n2c(twist_vec(i));
            tempWing.panel.section{ii}.chord  = n2c(chord_vec(i));
        end
    end
    
    % Sort wing panels from left tip to right tip
    tempWing.panel.section = tempWing.panel.section([2*noWingSpanEle:-1:(noWingSpanEle+1),2:noWingSpanEle]);
    inputS.VAP.vehicle.wing = tempWing;
    
    [ref_area, panelArea] = fcnsWINGAREA(wing_y,chord_vec);
    inputS.VAP.vehicle.ref_area = ref_area;
    inputS.VAP.vehicle.ref_span = ref_span;   
    
    %% VERTICAL STABILISER
    switch verticalTailType 
        case 'noVerticalTail'
            tempVert = {};
            
        case 'rectangular'
            tempVert            = tempWing;     % Make a main wing copy 
            tempVert.panel.section    = {};           % Delete sections 
            tempVert.symmetry            = n2c(symmetry);
            tempVert.incidence           = 0;
            tempVert.chordwise_elements  = n2c(chordwise_elements_vt);
            tempVert.panel.strip_airfoil = n2c(VerticalTail_Airfoil);
            tempVert.panel.spanwise_elements = '1';

            z_space_vt = linspace(0,1,noVTSpanEle);
            
            L_vt        = l_vt*c_root; % Now position from the root leading edge 
            S_vt        = volumeTail*ref_span*ref_area./L_vt;
            tail_span   =  sqrt(ar_vt/2*S_vt);
            c_root_vt   = S_vt/(1+taper_vt)/tail_span;
            if isnan(sweep_vt)
                sweep_vt    = atand(-c_root_vt*(taper_vt-1)/tail_span);
            end
            ref_cmac_vt = 2*c_root_vt*2/3*((1+taper_vt+taper_vt^2)/(1+taper_vt));
            z_cmac_vt   = tail_span/6*(1+2*taper_vt)/(1+taper_vt);
            x_ac_vt     = z_cmac_vt*tand(sweep_vt)+0.25*ref_cmac_vt;
            if wingtips
                x_vtStart = wing_x(end);
                vert_y    = wing_y(end)*ones(1,noVTSpanEle);
                z_vtStart = wing_z(end);
            else
                x_vtStart   = L_vt;
                vert_y = zeros(1,noVTSpanEle);
            end
            
                      
            vert_z          = z_vtStart + tail_span*z_space_vt;
            vert_x          = x_vtStart + (vert_z-z_vtStart)*tand(sweep_vt);
            vertChord_vec   = c_root_vt + c_root_vt*z_space_vt*(taper_vt-1);
            
            for i = 1:noVTSpanEle
                tempVert.panel.section{i}.wing_x = n2c(vert_x(i));
                tempVert.panel.section{i}.wing_y = n2c(vert_y(i));
                tempVert.panel.section{i}.wing_z = n2c(vert_z(i));
                tempVert.panel.section{i}.twist  = 0;
                tempVert.panel.section{i}.chord  = vertChord_vec(i);
            end
            
            if vt_doublePanel % For the PMAV double tail design 
                
                S_vt2 = volumeTail2*ref_span*ref_area./L_vt;
                tail_span2 = sqrt(ar_vt2/2*S_vt2);
                c_root_vt2 = S_vt2/(1+taper_vt2)/tail_span2;
                if isnan(sweep_vt2)
                    sweep_vt2    = atand(-c_root_vt2*(taper_vt2-1)/tail_span2);
                end
                ref_cmac_vt2 = 2*c_root_vt2*2/3*((1+taper_vt2+taper_vt2^2)/(1+taper_vt2));
                z_cmac_vt2   = tail_span2/6*(1+2*taper_vt2)/(1+taper_vt2);
                x_ac_vt2     = z_cmac_vt2*tand(sweep_vt2)+0.25*ref_cmac_vt2;
                
                vert2_z         = flip(z_vtStart - tail_span2*z_space_vt);
                vert2_x         = x_vtStart + abs(vert2_z-z_vtStart)*tand(sweep_vt2);
                vertChord2_vec  = flip(c_root_vt2 + c_root_vt2*z_space_vt*(taper_vt2-1));

                for i = 1:noVTSpanEle
                    % Lower
                    tempVert.panel.section{i}.wing_x = n2c(vert2_x(i));
                    tempVert.panel.section{i}.wing_y = n2c(vert_y(i));
                    tempVert.panel.section{i}.wing_z = n2c(vert2_z(i));
                    tempVert.panel.section{i}.twist  = 0;
                    tempVert.panel.section{i}.chord  = vertChord2_vec(i);
                    % Upper 
                    tempVert.panel.section{i+noVTSpanEle}.wing_x = n2c(vert_x(i));
                    tempVert.panel.section{i+noVTSpanEle}.wing_y = n2c(vert_y(i));
                    tempVert.panel.section{i+noVTSpanEle}.wing_z = n2c(vert_z(i));
                    tempVert.panel.section{i+noVTSpanEle}.twist  = 0;
                    tempVert.panel.section{i+noVTSpanEle}.chord  = vertChord_vec(i);  
                end
                
                if wingtips 
                    tempVert2 = tempVert;
                    for i = 1:noVTSpanEle
                        tempVert2.panel.section{i}.wing_y = n2c(-vert_y(i));
                        tempVert2.panel.section{i+noVTSpanEle}.wing_y = n2c(-vert_y(i));
                    end
                    
                    tempVert.panel.section(noVTSpanEle) = [];
                    tempVert2.panel.section(noVTSpanEle) = [];
                    inputS.VAP.vehicle.wing = {tempWing, tempVert, tempVert2};
                    
                else
                    % Delete double up
                    tempVert.panel.section(noVTSpanEle) = [];
                    inputS.VAP.vehicle.wing = {tempWing, tempVert};
                end
                

                
                [VT_area, VTpanelArea] = fcnsWINGAREA([vert2_z(1:end-1), vert_z], [vertChord2_vec(1:end-1), vertChord_vec]);
            
                GEOMETRY.VT_area     = VT_area;
                GEOMETRY.VTpanelArea = VTpanelArea;
                GEOMETRY.span_vt     = tail_span2 + tail_span;
                GEOMETRY.x_ac_vt     = x_ac_vt;
                GEOMETRY.z_ac_vt     = z_cmac_vt;
                GEOMETRY.cmac_vt     = ref_cmac_vt; 
                GEOMETRY.x_ac_vt2    = x_ac_vt2;
                GEOMETRY.z_ac_vt2    = z_cmac_vt2;
                GEOMETRY.cmac_vt2    = ref_cmac_vt2; 

                GEOMETRY.vert_x         = [vert2_x(1:end-1), vert_x];
                GEOMETRY.vert_z         = [vert2_z(1:end-1), vert_z]; 
                GEOMETRY.vertChord_vec  = [vertChord2_vec(1:end-1), vertChord_vec];
                
            else
                [VT_area, VTpanelArea] = fcnsWINGAREA(vert_z, vertChord_vec);
            
                GEOMETRY.VT_area     = VT_area;
                GEOMETRY.VTpanelArea = VTpanelArea;
                GEOMETRY.span_vt     = tail_span;
                GEOMETRY.x_ac_vt     = x_ac_vt;
                GEOMETRY.z_ac_vt     = z_cmac_vt;
                GEOMETRY.cmac_vt     = ref_cmac_vt; 

                GEOMETRY.vert_x         = vert_x;
                GEOMETRY.vert_z         = vert_z; 
                GEOMETRY.vertChord_vec  = vertChord_vec;            
                
                inputS.VAP.vehicle.wing = {tempWing, tempVert};
            end
 
        otherwise 
            error('Vertical Tail type not recognised - TR')
    end

    %% HORIZONTAL STABILISER 
    switch horizontalTailType
        case 'noHorizontalTail'
            tempHori = {};

        case 'rectangular'
            tempHori = tempWing;     % Make a main wing copy 
            tempHori.panel.section       = {};           % Delete sections
            tempHori.symmetry            = n2c(symmetry);
            tempHori.incidence           = n2c(incidence_ht);
            tempHori.chordwise_elements  = n2c(chordwise_elements_ht);
            tempHori.panel.strip_airfoil = n2c(HorizontalTail_Airfoil);
            tempHori.panel.spanwise_elements = '1';
            
            switch HTSpanSpacing
                case 'linear'
                    y_space_ht     = linspace(0,1,noHTSpanEle); 

                case 'halfsine'
                    thetay_ht      = linspace(0,pi,noHTSpanEle);
                    y_space_ht     = sin(thetay_ht/2);

                otherwise 
                    error('HT Spanwise spacing type not recognised - TR')
            end
            
            L_ht            = l_vt*c_root; % BASED OFF VERTICAL TAIL LOCATION
            S_ht            = volumeHtail/2*ref_span*ref_area./L_ht;
            span_ht         = sqrt(ar_ht*S_ht);
            c_root_ht       = S_ht/(1+taper_ht)/span_ht;
            if isnan(sweep_ht)
                sweep_ht        = atand(-2*c_root_ht*(taper_ht-1)/span_ht);
            end
            
            cmac_ht         = c_root_ht*2/3*((1+taper_ht+taper_ht^2)/(1+taper_ht));
            y_cmac_ht       = span_ht/6*(1+2*taper_ht)/(1+taper_ht);
            x_ac_ht         = y_cmac_ht*tand(sweep_ht)+0.25*cmac_ht;
            xStart_ht       = L_ht;
            
            hori_y          = span_ht/2*y_space_ht;
            hori_x          = xStart_ht + hori_y*tand(sweep_ht);
            hori_z          = zStart_ht + hori_y.*tand(dihedral_ht);
            horiTwist_vec   = twist_ht*y_space_ht;
            horChord_vec    = c_root_ht + c_root_ht*y_space_ht*(taper_ht-1);            
            
            for side = 0:1
                RHS = 1;
                if side
                    RHS = -1;
                end
                for i = 1:noHTSpanEle
                    ii = side*noHTSpanEle + i;
                    tempHori.panel.section{ii}.wing_x = n2c(hori_x(i));
                    tempHori.panel.section{ii}.wing_y = n2c(RHS*hori_y(i));
                    tempHori.panel.section{ii}.wing_z = n2c(hori_z(i));
                    tempHori.panel.section{ii}.twist  = n2c(horiTwist_vec(i));
                    tempHori.panel.section{ii}.chord  = n2c(horChord_vec(i));
                end
            end    
            
            % Sort HT panels from right tip to left tip 
            tempHori.panel.section = tempHori.panel.section([2*noHTSpanEle:-1:(noHTSpanEle+1),2:noHTSpanEle]);
            
            if strcmp(verticalTailType, 'noVerticalTail')
                inputS.VAP.vehicle.wing = {tempWing, tempHori};
            else
                inputS.VAP.vehicle.wing = {tempWing, tempHori, tempVert};
            end
            
            [HT_area, HTpanelArea] = fcnsWINGAREA(hori_y,horChord_vec);
            
            GEOMETRY.HT_area     = HT_area;
            GEOMETRY.HTpanelArea = HTpanelArea;
            GEOMETRY.span_ht     = span_ht;
            GEOMETRY.x_ac_ht     = x_ac_ht;
            GEOMETRY.cmac_ht     = cmac_ht; 
            GEOMETRY.aspect_ratio_ht = span_ht^2./HT_area;
            
            GEOMETRY.hori_y     = hori_y;
            GEOMETRY.hori_x     = hori_x;
            GEOMETRY.hori_z     = hori_z; 
            GEOMETRY.horChord_vec  = horChord_vec;
     
        otherwise 
            error('Horizontal Tail type not recognised - TR')
    end

    %% FUSELAGE
    switch fuselageType 
        case 'noFuselage'
            tempFuse = {};
            
        case 'rectangular'
            tempFuse1 = tempWing;     % Make a main wing copy 
            tempFuse1.panel.section       = {};           % Delete sections
            tempFuse1.symmetry            = n2c(symmetry);
            tempFuse1.incidence           = n2c(incidence_fu);
            tempFuse1.chordwise_elements  = n2c(chordwise_elements_fu);
            tempFuse1.panel.strip_airfoil = n2c(Fuselage_airfoil);
            tempFuse1.panel.spanwise_elements = '1';
            tempFuse2 = tempFuse1;
            
            switch FuSpanSpacing
                case 'linear'
                    y_space_fu     = linspace(0,1,noFuSpanEle); 

                case 'halfsine'
                    thetay_fu      = linspace(0,pi,noFuSpanEle);
                    y_space_fu     = sin(thetay_fu/2);

                otherwise 
                    error('Fuselage spanwise spacing type not recognised - TR')
            end
             
            fuseLength = c_root; % Overrite for optimisation
            fuseDiam = fuseAR*fuseLength;
            if isnan(zStart_fu)
                zStart_fu = -fuseDiam/2-0.001; % For the opitmisation process 
            end
            fuse_w = fuseDiam/2*y_space_fu;
            fuse_x = xStart_fu + fuse_w*tand(sweep_fu);
            fuseChord_vec = fuseLength + fuseLength*y_space_fu*(taper_fu -1);
            
            cmac_fu   = fuseLength*2/3*((1+taper_fu+taper_fu^2)/(1+taper_fu));
            w_cmac_fu = fuseDiam/6*(1+2*taper_fu)/(1+taper_fu);
            x_ac_fu   = w_cmac_fu*tand(sweep_fu)+0.25*cmac_fu;
            
            fuse_z_hori = zStart_fu*ones(1,noFuSpanEle);
            fuse_z_vert = zStart_fu + fuse_w;
            fuse_z_vert2 = zStart_fu + flip(fuse_w) - fuseDiam/2;
            
            % Horizontal Panel [1]
            for side = 0:1
                RHS = 1;
                if side
                    RHS = -1;
                end
                for i = 1:noFuSpanEle
                    ii = side*noFuSpanEle + i;
                    tempFuse1.panel.section{ii}.wing_x = n2c(fuse_x(i));
                    tempFuse1.panel.section{ii}.wing_y = n2c(RHS*fuse_w(i));
                    tempFuse1.panel.section{ii}.wing_z = n2c(fuse_z_hori(i));
                    tempFuse1.panel.section{ii}.twist  = '0';
                    tempFuse1.panel.section{ii}.chord  = n2c(fuseChord_vec(i));
                end
            end
            
            % Vertical Panel [2]
            for BOT = 0:1
                for i = 1:noFuSpanEle
                    ii = BOT*noFuSpanEle + i;
                    tempFuse2.panel.section{ii}.wing_x = n2c(fuse_x(i));
                    tempFuse2.panel.section{ii}.wing_y = '0';
                    if BOT
                        tempFuse2.panel.section{ii}.wing_z = n2c(fuse_z_vert2(i)-0.001);
                    else
                        tempFuse2.panel.section{ii}.wing_z = n2c(fuse_z_vert(i)-0.001);
                    end
                    tempFuse2.panel.section{ii}.twist  = '0';
                    tempFuse2.panel.section{ii}.chord  = n2c(fuseChord_vec(i));
                end
            end
            
            % Sort panels from right tip to left tip (or top to bottom)
            tempFuse1.panel.section = tempFuse1.panel.section([2*noFuSpanEle:-1:(noFuSpanEle+1),2:noFuSpanEle]);
            tempFuse2.panel.section = tempFuse2.panel.section([2*noFuSpanEle:-1:(noFuSpanEle+1),2:noFuSpanEle]);
            
            if strcmp(verticalTailType, 'noVerticalTail')
                if strcmp(horizontalTailType,'noHorizontalTail')
                    inputS.VAP.vehicle.wing = {tempWing, tempFuse1, tempFuse2};
                else
                    inputS.VAP.vehicle.wing = {tempWing, tempHori, tempFuse1, tempFuse2};
                end
            else
                if wingtips 
                    if strcmp(horizontalTailType,'noHorizontalTail')
                        inputS.VAP.vehicle.wing = {tempWing, tempVert, tempVert2 tempFuse1, tempFuse2};
                    else
                        inputS.VAP.vehicle.wing = {tempWing, tempHori, tempVert, tempVert2 tempFuse1, tempFuse2};
                    end
                else
                    if strcmp(horizontalTailType,'noHorizontalTail')
                        inputS.VAP.vehicle.wing = {tempWing, tempVert, tempFuse1, tempFuse2};
                    else
                        inputS.VAP.vehicle.wing = {tempWing, tempHori, tempVert, tempFuse1, tempFuse2};
                    end
                end
            end
            
            [Fu_area, FupanelArea] = fcnsWINGAREA(fuse_w,fuseChord_vec);
            
            GEOMETRY.fuseDiam    = fuseDiam;
            GEOMETRY.Fu_area     = Fu_area;
            GEOMETRY.FupanelArea = FupanelArea;
            GEOMETRY.x_ac_fu     = x_ac_fu;
            GEOMETRY.cmac_fu     = cmac_fu; 
            
            GEOMETRY.fuse_w     = fuse_w;
            GEOMETRY.fuse_x     = fuse_x;
            GEOMETRY.fuse_z_hori   = fuse_z_hori; 
            GEOMETRY.fuse_z_vert   = fuse_z_vert;
            GEOMETRY.fuseChord_vec = fuseChord_vec;
            
        otherwise
            error('Fuselage Type not recognised - TR')
    end

    %% PROPELLER 
    if strcmp(propellerName,'noProp') || fixedLift
        inputS.VAP.vehicle.rotor = {};
        
    else
        tempRotor = inputS.VAP.vehicle.rotor;
        tempRotor.collective         = n2c(collective);
        tempRotor.ref_diam           = n2c(propDiam);
        tempRotor.rotation_direction = n2c(rotation_direction);
        tempRotor.blades             = n2c(blades);
        tempRotor.rpm                = n2c(rpm);
        tempRotor.veh_x_hub          = n2c(veh_x_hub);
        tempRotor.veh_z_hub          = n2c(zStart_fu); % Prop in line with fuse 
        
        propRad  = propDiam/2;
        propGeom = dlmread(strcat('propellers/',propellerName,'.dat'));
        radialRaw = propGeom(:,1);
        chordRaw  = propGeom(:,2);
        pitchRaw  = propGeom(:,3);
        leadingEdgeRaw = fcnPROPLE(propellerName);
        
        radial_vec      = linspace(radialRaw(1) ,0.98 ,noRotorSpanEle);
        chordr_vec      = interp1(radialRaw,chordRaw,radial_vec,'spline')*propRad;
        pitch_vec       = interp1(radialRaw,pitchRaw,radial_vec,'spline') + pitchShift;
        leadingedge_vec = interp1(leadingEdgeRaw(:,1),leadingEdgeRaw(:,2),radial_vec,'spline')*propRad;
        
        tempPanel.strip_airfoil     = n2c(propAirfoil);
        tempPanel.spanwise_elements = '1';
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
    
    %% WEIGHT
    GEOMETRY.ref_span   = ref_span;
    GEOMETRY.ref_area   = ref_area;
    GEOMETRY.ref_cmac   = ref_cmac;    
    GEOMETRY.aspect_ratio = ref_span^2./ref_area;
    GEOMETRY.panelArea  = panelArea;
    GEOMETRY.wing_y     = wing_y;
    GEOMETRY.wing_x     = wing_x;
    GEOMETRY.wing_z     = wing_z; 
    GEOMETRY.chord_vec  = chord_vec;
    GEOMETRY.x_ac       = x_ac;
    GEOMETRY.y_ac       = y_cmac;   
    GEOMETRY.cog_z      = 0; % Will need to change when fuselage is added 5/10/2021 
    
    % Define weight and CoG, accounting for manual overrides
    if isnan(weightExplicit) 
        [WEIGHT, GEOMETRY] = fcnWEIGHT(VARIABLES, GEOMETRY, SETTINGS);
        if ~isnan(cogExplicitx)
            GEOMETRY.cog_x = cogExplicitx;
        end
    else 
        WEIGHT = weightExplicit;
        if isnan(cogExplicitx)
           [~, GEOMETRY] = fcnWEIGHT(VARIABLES, GEOMETRY, SETTINGS); 
        else
            GEOMETRY.cog_x = cogExplicitx;
        end
    end
    
    if strcmp(speed,'nan') || fixedLift
        inputS.VAP.vehicle.weight = n2c(WEIGHT*9.81); % Weight force    
    end
        
    %% OUTPUT
    inputFilename = strcat('tempVAPINPUT_',num2str(ID));  % Do I need to include the ID tag here? 
    %inputFilename = strcat(name,'.vap');
    fcnSTRUCT2XML(inputS, inputFilename)
end

%% Subroutines
function [ref_area, panelArea] = fcnsWINGAREA(wing_y,chord_vec)
% Calculates the wing area based on the sum of trapezedoil sections 

    panelArea = zeros(1,length(wing_y)-1);

    for i = 1:length(wing_y)-1
        dy = wing_y(i+1) - wing_y(i);
        panelArea(i) = 0.5*dy*(chord_vec(i+1) + chord_vec(i));
    end
      
    ref_area = 2*sum(panelArea); % Both sides of wing
end