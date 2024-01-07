function STABILITY = fcnSTABILITY(inputFilename, GEOMETRY, SETTINGS, VAP_IN, OUTP_first, SURF_first, COND_first, maxTimeStab, Id)
% This function numerically calculates the stability derivatives 
% This is to be used in the fcnMAINSOLVER_TOM function. See
% fcnSTABILITY_ANALYSIS to produce the validation plots
% 
% CMx, CMy produced by VAP, coordinates defined by input file?
% Shift the vehicle coordinates so that global coordinates line up with
% CoG?
%
% Vehicle origin is shifted so that global origin cooincides with COG
% 
% Static Criteria: Cma < 0, Cnb > 0, Clb < 0
% 
% Author: Tom Ryan
% 7/09/2021
%
% Status: This is the bain of my existence 12/10/2021

    % Extract variables 
    fcnSTRUCT2VARS(GEOMETRY)
    
    % Check if there is a predefined CG 
    if ~isnan(SETTINGS.cogExplicitx)
        cog_x = SETTINGS.cogExplicitx;
        cog_z = SETTINGS.cogExplicitz;
    else 
        cog_z = 0;
    end
    
    %% SET UP 
    % MAIN WING 
    % Difference vectors 
    dChord = [diff(chord_vec)/2, 0];
    dx     = [diff(wing_x)/2, 0];
    dy     = [diff(wing_y)/2, 0];
    dz     = [diff(wing_z)/2, 0];
    
    % Centre of wing panel geometry values 
    wingc_x     = wing_x(1:end-1) + dx(1:end-1);
    wingtemp_y     = wing_y(1:end-1) + dy(1:end-1);
    wingc_z     = wing_z(1:end-1) + dz(1:end-1);
    chordc_vec  = chord_vec(1:end-1) + dChord(1:end-1);
    
    % Reorganise panel area vector (current vectors are for one side only)
    fcnWING_REORDER(panelArea, wingc_x, wingc_z, chordc_vec)
    wingc_y = [-flip(wingtemp_y), wingtemp_y];
    
    % 1/4 Chord coordinates 
    x_025 = wingc_x + 0.25*chordc_vec;
    y_025 = wingc_y;
    z_025 = wingc_z;
    
    % Wing index
    wing_index = 1:SETTINGS.noWingDVE;
    
    % HORIZONTAL TAIL  
    if ~strcmp(SETTINGS.horizontalTailType,'noHorizontalTail')
        
        % Difference vectors 
        dChord_ht = [diff(horChord_vec)/2, 0];
        dx_ht     = [diff(hori_x)/2, 0];
        dy_ht     = [diff(hori_y)/2, 0];

        % Centre of wing panel geometry values 
        horic_x     = hori_x(1:end-1) + dx_ht(1:end-1);
        horic_y     = hori_y(1:end-1) + dy_ht(1:end-1);
        horChordc_vec  = horChord_vec(1:end-1) + dChord_ht(1:end-1);

        % Reorganise panel area vector (current vectors are for one side only)
        fcnWING_REORDER(HTpanelArea, horic_x, horic_y, horChordc_vec)

        % 1/4 Chord locations 
        x_025_ht = horic_x + 0.25*horChordc_vec;    
        y_025_ht = horic_y; 
        
        % HT Index 
        start_ht = SETTINGS.noWingDVE + 1;
        end_ht   = SETTINGS.noWingDVE + SETTINGS.noHTDVE;
        ht_index = start_ht:end_ht;
    end
    
    % VERTICAL TAIL
    if ~strcmp(SETTINGS.verticalTailType,'noVerticalTail')
        
        % Difference vectors 
        dChord_vt = [diff(vertChord_vec)/2, 0];
        dz_vt     = [diff(vert_z)/2, 0];
        dx_vt     = [diff(vert_x)/2, 0];

        % Centre of wing panel geometry values 
        vertc_x     = vert_x(1:end-1) + dx_vt(1:end-1);
        vertc_z     = vert_z(1:end-1) + dz_vt(1:end-1);
        vertChordc_vec  = vertChord_vec(1:end-1) + dChord_vt(1:end-1);

        % 1/4 Chord locations 
        z_025_vt = vertc_z; 
        x_025_vt = vertc_x + 0.25*vertChordc_vec;
       
        % VT Index
        end_ht   = SETTINGS.noWingDVE + SETTINGS.noHTDVE;
        end_vt = SETTINGS.noWingDVE + SETTINGS.noHTDVE + SETTINGS.noVTDVE;
        if strcmp(SETTINGS.horizontalTailType,'noHorizontalTail')
            end_ht = SETTINGS.noWingDVE;  
            end_vt = SETTINGS.noWingDVE + SETTINGS.noVTDVE;
        end
        start_vt = end_ht + 1;
        vt_index = start_vt:end_vt;
        
        if SETTINGS.wingtips
            vt2_index = vt_index + SETTINGS.noVTDVE;
            
            y_025_vt1 = wing_y(end);
            y_025_vt2 = -wing_y(end);
        end
    end
    
    % FUSELAGE 
    % [1] Hori Panel , [2] Vert Panel
    if ~strcmp(SETTINGS.fuselageType,'noFuselage')
        noFuseDVE = 2*(SETTINGS.noFuSpanEle - 1);
        
        dw_fu = [diff(fuse_w)/2,0];
        dx_fu = [diff(fuse_x)/2,0];
        dChord_fu = [diff(fuseChord_vec)/2,0];
        
        if isnan(SETTINGS.zStart_fu) 
            fu1c_z = (-fuseDiam/2 - 0.001)*ones(noFuseDVE,1);
        else
            fu1c_z = SETTINGS.zStart_fu*ones(noFuseDVE,1);
        end
        fu1temp_y = fuse_w(1:end-1) + dw_fu(1:end-1);
        fu1c_y = [-flip(fu1temp_y), fu1temp_y];
        
        fuseTemp_z = fuse_w(1:end-1) + dw_fu(1:end-1);
        fuseTemp_x = fuse_x(1:end-1) + dx_fu(1:end-1);
        fuseChordc_vec = fuseChord_vec(1:end-1) + dChord_fu(1:end-1);
        if isnan(SETTINGS.zStart_fu) 
            fu2c_z = [  -fuseDiam/2 - 0.001 + fuseTemp_z - fuseDiam/2, -fuseDiam/2 - 0.001 + fuseTemp_z];
        else
            fu2c_z = [ SETTINGS.zStart_fu + fuseTemp_z - fuseDiam/2, SETTINGS.zStart_fu + fuseTemp_z];
        end
        
        fcnWING_REORDER(fuseChordc_vec, FupanelArea, fuseTemp_x)
        
        x_025_fu1 = fuseTemp_x;
        x_025_fu2 = x_025_fu1;
        y_025_fu1 = fu1c_y;
        z_025_fu1 = fu1c_z;
        z_025_fu2 = fu2c_z;
        
        % Fuselage is always the last two panels
        % Indexing done inside derivative calculation 
        
    end
    
    %% GENERAL
    % Set up 
    StabilityInput = strcat('TEMPinStabiltiy_',num2str(Id));
    dalpha = [0 -1];
    
    % Load TEMP VAP input file generated by fcnGENERATE_INPUT() in a STRUCT
    inputS = fcnXML2STRUCT(inputFilename);
    
    % Read initial beta and alpha variabels
    alpha_0 = str2num(inputS.VAP.vehicle.alpha.Text);
    beta_0  = str2num(inputS.VAP.vehicle.beta.Text);
    inputS.VAP.settings.maxtime = n2c(maxTimeStab);
    
    
    %% Cma 
    if SETTINGS.stabSelection(1)
        First = 1;     
        for i = 1:length(dalpha)

            % Change Angle of attack 
            if First 
                OUTP = OUTP_first;
                COND = COND_first;
                SURF = SURF_first;
                First = 0;
            else
                inputS.VAP.vehicle.alpha = n2c(alpha_0 + dalpha(i));
                fcnSTRUCT2XML(inputS, StabilityInput)
                [OUTP, COND, ~, ~, ~, SURF, ~, ~, ~] = fcnVAP_MAIN(StabilityInput, VAP_IN);
            end

            Qinf = 0.5*OUTP.valDENSITY*COND.vecVEHVINF^2;

            % MAIN WING
            % Extract Moment about 1/4 chord for for given AoA (area normalised)
            M_025 = OUTP.vecCMDIST(wing_index) .*cos(SURF.vecDVEROLL(wing_index)) .*chordc_vec' .*SURF.vecDVEAREA(wing_index) .*Qinf; % *SURF.vecDVEAREA(wing_index)

            % Moment generated by wing normal force
            Mcg_wing  = OUTP.vecCNDIST(wing_index) .*cos(SURF.vecDVEROLL(wing_index)) .*SURF.vecDVEAREA(wing_index) .*Qinf .*(cog_x- x_025'); %
            
            % Moment from drag
            Mcg_wing_D = OUTP.vecCDPDIST(wing_index) .*cos(COND.vecVEHALPHA) .*SURF.vecDVEAREA(wing_index) .*Qinf .*(z_025' - cog_z);

            % Weighted average based on panel area 
            Cmcg(i) = sum(M_025 + Mcg_wing + Mcg_wing_D)/(Qinf *ref_area *ref_cmac);
            
            % VERTICAL TAIl 
            if ~strcmp(SETTINGS.verticalTailType,'noVerticalTail')
            
                % Moment from drag 
                M_cg_vt_Da = OUTP.vecCDPDIST(vt_index) .*cos(COND.vecVEHALPHA) .* SURF.vecDVEAREA(vt_index) .* Qinf .* (z_025_vt' - cog_z);
                M_cg_vt_Db = OUTP.vecCDPDIST(vt_index) .*sin(COND.vecVEHALPHA) .* SURF.vecDVEAREA(vt_index) .* Qinf .* (cog_x - x_025_vt');
                
                % Add on VT effect 
                Cmcg(i) = Cmcg(i) + sum(M_cg_vt_Da + M_cg_vt_Db)/(Qinf *ref_area *ref_cmac);        

                if SETTINGS.wingtips
                    % Moment from drag 
                    M_cg_vt2_Da = OUTP.vecCDPDIST(vt2_index) .*cos(COND.vecVEHALPHA) .* SURF.vecDVEAREA(vt2_index) .* Qinf .* (z_025_vt' - cog_z);
                    M_cg_vt2_Db = OUTP.vecCDPDIST(vt2_index) .*sin(COND.vecVEHALPHA) .* SURF.vecDVEAREA(vt2_index) .* Qinf .* (cog_x - x_025_vt');

                    % Add on VT effect 
                    Cmcg(i) = Cmcg(i) + sum(M_cg_vt2_Da + M_cg_vt2_Db)/(Qinf *ref_area *ref_cmac);
                end
            end

            % HORIZONTAL STABILISER
            if ~strcmp(SETTINGS.horizontalTailType,'noHorizontalTail')

                % Moment about 1/4 chord
                M_025_ht = OUTP.vecCMDIST(ht_index) .*cos(SURF.vecDVEROLL(ht_index)) .*SURF.vecDVEAREA(ht_index) .*HTpanelArea' *Qinf;

                % Moment due to lift moment 
                M_cg_ht = OUTP.vecCNDIST(ht_index) .*cos(SURF.vecDVEROLL(ht_index)) .*SURF.vecDVEAREA(ht_index) .*Qinf .*(cog_x-x_025_ht');

                % Add on HT effect 
                Cmcg(i)  = Cmcg(i) + sum(M_025_ht + M_cg_ht)/(Qinf *ref_area *ref_cmac);
            end
            
            % FUSELAGE
            if ~strcmp(SETTINGS.fuselageType,'noFuselage')
                fu1_start = length(OUTP.vecCNDIST) - 2*noFuseDVE + 1;
                fu2_start = fu1_start + noFuseDVE;
                fu1_index = fu1_start:(fu2_start-1);
                fu2_index = fu2_start:length(OUTP.vecCNDIST);
                
                % Fuse [1] Hori Panel
                M_cg_fu1_Da = OUTP.vecCDPDIST(fu1_index) .* cos(SURF.vecDVEPITCH(fu1_index)) .*SURF.vecDVEAREA(fu1_index) .*Qinf .*(z_025_fu1 - cog_z);
                M_cg_fu1_Db = OUTP.vecCDPDIST(fu1_index) .* sin(SURF.vecDVEPITCH(fu1_index)) .*SURF.vecDVEAREA(fu1_index) .*Qinf .*(cog_x - x_025_fu1');
                M_025_fu1   = OUTP.vecCMDIST(fu1_index)  .* SURF.vecDVEAREA(fu1_index) .* fuseChordc_vec' .*Qinf;
                
                % Fuse [2] Vert Panel
                M_cg_fu2_Da = OUTP.vecCDPDIST(fu2_index) .* cos(SURF.vecDVEPITCH(fu2_index)) .*SURF.vecDVEAREA(fu2_index) .*Qinf .*(z_025_fu2' - cog_z);
                M_cg_fu2_Db = OUTP.vecCDPDIST(fu2_index) .* sin(SURF.vecDVEPITCH(fu2_index)) .*SURF.vecDVEAREA(fu2_index) .*Qinf .*(cog_x - x_025_fu2');
                
                Cmcg(i) = Cmcg(i) + sum(M_cg_fu1_Da + M_cg_fu1_Db + M_025_fu1 + M_cg_fu2_Da + M_cg_fu2_Db)/(Qinf *ref_area *ref_cmac); 
            end
        end   

        % p = polyfit(dalpha,Cmcg,1);
        p = polyfit([alpha_0, alpha_0 + dalpha(2)],Cmcg,1);
        STABILITY.Cma = p(1)*180/pi;
        STABILITY.Cmcg_vec = Cmcg;

        % Revert alpha 
        inputS.VAP.vehicle.alpha = n2c(alpha_0);
    end
    

    First = 1; 
    for i = 1:length(dalpha) %% Beta loop            
       
       
        % Change sideslip angle
        if First 
            OUTP = OUTP_first;
            COND = COND_first;
            SURF = SURF_first;
            First = 0;
        else
            inputS.VAP.vehicle.beta = n2c(beta_0 + dalpha(i));
            fcnSTRUCT2XML(inputS, StabilityInput)
            [OUTP, COND, ~, ~, ~, SURF, ~, ~, ~] = fcnVAP_MAIN(StabilityInput, VAP_IN);
        end
        
        Qinf = 0.5*OUTP.valDENSITY.*COND.vecVEHVINF^2;
        
        %% Clb
        if SETTINGS.stabSelection(2)

            % MAIN WING 
            l_cg_z = OUTP.vecCNDIST(wing_index) .*cos(SURF.vecDVEPITCH(wing_index)) .*sin(SURF.vecDVEROLL(wing_index)) .*SURF.vecDVEAREA(wing_index) .*Qinf .*(cog_z- z_025');
            
            l_cg_y = OUTP.vecCNDIST(wing_index) .*cos(SURF.vecDVEPITCH(wing_index)) .*cos(SURF.vecDVEROLL(wing_index)) .*SURF.vecDVEAREA(wing_index) .*Qinf .* -y_025'; %cog_y = 0

            l_025 = OUTP.vecCMDIST(wing_index) .*sin(SURF.vecDVEROLL(wing_index)) .*SURF.vecDVEAREA(wing_index) .*chordc_vec' .*Qinf;

            Clcg(i) = sum(l_cg_z + l_cg_y + l_025)/(Qinf *ref_span *ref_area);

            % VERTICAL STABILIZER 
            if ~strcmp(SETTINGS.verticalTailType,'noVerticalStabilizer')

               % Moment due to normal force 
               l_cg_vt = OUTP.vecCNDIST(vt_index) .*SURF.vecDVEAREA(vt_index) .*Qinf .*( cog_z - z_025_vt');

               Clcg(i) = Clcg(i) + sum(l_cg_vt)/(Qinf*ref_span *ref_area);
               
               if SETTINGS.wingtips 
                   
                    % Moment from drag (wintips only 
                    l_cg_vt_d = OUTP.vecCDPDIST(vt_index) .*cos(SURF.vecDVEYAW(vt_index)) .*sin(SURF.vecDVEPITCH(vt_index)) .*SURF.vecDVEAREA(vt_index) .* Qinf .*-y_025_vt1;
                    l_vg_vt2_d = OUTP.vecCDPDIST(vt2_index) .*cos(SURF.vecDVEYAW(vt2_index)) .*sin(SURF.vecDVEPITCH(vt2_index)) .*SURF.vecDVEAREA(vt2_index) .*Qinf .*-y_025_vt2;
                   
                    % Moment due to normal force 
                    l_cg_vt2 = OUTP.vecCNDIST(vt2_index) .*SURF.vecDVEAREA(vt2_index) .*Qinf .*( cog_z - z_025_vt');

                    Clcg(i) = Clcg(i) + sum(l_cg_vt_d + l_vg_vt2_d + l_cg_vt2)/(Qinf*ref_span *ref_area);
               end
            end
            
            % FUSELAGE 
            % [1] Hori Panel, [2] Vert Panel
            if ~strcmp(SETTINGS.fuselageType,'noFuselage')
                fu1_start = length(OUTP.vecCNDIST) - 2*noFuseDVE + 1;
                fu2_start = fu1_start + noFuseDVE;
                fu1_index = fu1_start:(fu2_start-1);
                fu2_index = fu2_start:length(OUTP.vecCNDIST);
                
                l_cg_y_fu1 = OUTP.vecCNDIST(fu1_index) .*SURF.vecDVEAREA(fu1_index) .* Qinf .*  - y_025_fu1';
                
                l_cg_fu2 = OUTP.vecCNDIST(fu2_index) .*SURF.vecDVEAREA(fu2_index) .*Qinf .*( cog_z - z_025_fu2');
                
                Clcg(i) = Clcg(i) + sum(l_cg_y_fu1 + l_cg_fu2)/(Qinf*ref_span *ref_area);
            end  
        end
        
        %% CNb  
        if SETTINGS.stabSelection(3)

            % MAIN WING 
            N_cg_wing_Na = OUTP.vecCNDIST(wing_index) .* sin(SURF.vecDVEPITCH(wing_index)) .* cos(SURF.vecDVEYAW(wing_index)) .* SURF.vecDVEAREA(wing_index) .* Qinf .* -y_025';
            N_cg_wing_Nb = OUTP.vecCNDIST(wing_index) .* sin(SURF.vecDVEPITCH(wing_index)) .* sin(SURF.vecDVEYAW(wing_index)) .* SURF.vecDVEAREA(wing_index) .* Qinf .*(x_025' - cog_x);

            N_cg_wing_Da = OUTP.vecCDPDIST(wing_index) .* cos(SURF.vecDVEPITCH(wing_index)).* cos(SURF.vecDVEYAW(wing_index)) .* SURF.vecDVEAREA(wing_index) .* Qinf .* -y_025';
            N_cg_wing_Db = OUTP.vecCDPDIST(wing_index) .* cos(SURF.vecDVEPITCH(wing_index)).* sin(SURF.vecDVEYAW(wing_index)) .* SURF.vecDVEAREA(wing_index) .* Qinf .*(x_025' - cog_x);

            CNcg(i) = sum(N_cg_wing_Na + N_cg_wing_Nb + N_cg_wing_Da + N_cg_wing_Db) / (Qinf * ref_area * ref_span); 

            % VERTICAL STABILIZER 
            if ~strcmp(SETTINGS.verticalTailType, 'noVerticalTail') 

                N_025_vt = OUTP.vecCMDIST(vt_index) .* SURF.vecDVEAREA(vt_index) .* vertChordc_vec' .* Qinf;

                N_cg_vt = OUTP.vecCNDIST(vt_index) .* SURF.vecDVEAREA(vt_index) .*Qinf .*(x_025_vt' - cog_x);

                CNcg(i) = CNcg(i) + sum(N_025_vt + N_cg_vt)/(Qinf *ref_area * ref_span);
                
                if SETTINGS.wingtips 
                     
                    % Moment due to drag 
                    N_025_vt_da = OUTP.vecCDPDIST(vt_index) .* cos(SURF.vecDVEPITCH(vt_index)) .* cos(SURF.vecDVEYAW(vt_index)).* SURF.vecDVEAREA(vt_index) .* Qinf .* -y_025_vt1;
                    N_025_vt_db = OUTP.vecCDPDIST(vt_index) .* cos(SURF.vecDVEPITCH(vt_index)) .* sin(SURF.vecDVEYAW(vt_index)).* SURF.vecDVEAREA(vt_index) .* Qinf .* (x_025_vt' - cog_x);
                    
                    N_025_vt2_da = OUTP.vecCDPDIST(vt2_index) .* cos(SURF.vecDVEPITCH(vt2_index)) .* cos(SURF.vecDVEYAW(vt2_index)).* SURF.vecDVEAREA(vt2_index) .* Qinf .* -y_025_vt1;
                    N_025_vt2_db = OUTP.vecCDPDIST(vt2_index) .* cos(SURF.vecDVEPITCH(vt2_index)) .* sin(SURF.vecDVEYAW(vt2_index)).* SURF.vecDVEAREA(vt2_index) .* Qinf .* (x_025_vt' - cog_x);
                    
                    N_025_vt2 = OUTP.vecCMDIST(vt2_index) .* SURF.vecDVEAREA(vt2_index) .* vertChordc_vec' .* Qinf;

                    N_cg_vt2 = OUTP.vecCNDIST(vt2_index) .* SURF.vecDVEAREA(vt2_index) .*Qinf .*(x_025_vt' - cog_x);
                    
                    CNcg(i) = CNcg(i) + sum(N_cg_vt2 + N_025_vt2 + N_025_vt2_da + N_025_vt2_db + N_025_vt_da + N_025_vt_db)/(Qinf *ref_area * ref_span);
                end
            end
            % FUSELAGE 
            if ~strcmp(SETTINGS.fuselageType, 'noFuselage') 
                fu1_start = length(OUTP.vecCNDIST) - 2*noFuseDVE + 1;
                fu2_start = fu1_start + noFuseDVE;
                fu1_index = fu1_start:(fu2_start-1);
                fu2_index = fu2_start:length(OUTP.vecCNDIST);

                % Fuse [1] Hori Panel
                N_cg_fu1_Da = OUTP.vecCDPDIST(fu1_index) .* cos(SURF.vecDVEPITCH(fu1_index)) .* cos(SURF.vecDVEYAW(fu1_index)) .* SURF.vecDVEAREA(fu1_index) .* Qinf .* -y_025_fu1';
                N_cg_fu1_Db = OUTP.vecCDPDIST(fu1_index) .* cos(SURF.vecDVEPITCH(fu1_index)) .* sin(SURF.vecDVEYAW(fu1_index)) .* SURF.vecDVEAREA(fu1_index) .* Qinf .* (x_025_fu1' - cog_x);

                % Fuse [2] Vert Panel
                N_025_fu2 = OUTP.vecCMDIST(fu2_index) .*SURF.vecDVEAREA(fu2_index) .* fuseChordc_vec' .*Qinf; 
                N_cg_fu2 = OUTP.vecCNDIST(fu2_index) .* SURF.vecDVEAREA(fu2_index) .* Qinf .* (x_025_fu2' - cog_x);

                CNcg(i) = CNcg(i) + sum(N_cg_fu1_Da + N_cg_fu1_Db + N_025_fu2 + N_cg_fu2)/(Qinf * ref_area * ref_span);
            end

            % HORIZONTAL TAIL 
            if ~strcmp(SETTINGS.horizontalTailType, 'noHorizontalTail') 

                %N_cg_ht_Na = OUTP.vecCNDIST(ht_index) .* sin(SURF.vecDVEPITCH(ht_index)) .* cos(SURF.vecDVEYAW(ht_index)) .* SURF.vecDVEAREA(ht_index) .* Qinf .* -y_025_ht';
                %N_cg_ht_Nb = OUTP.vecCNDIST(ht_index) .* sin(SURF.vecDVEPITCH(ht_index)) .* sin(SURF.vecDVEYAW(ht_index)) .* SURF.vecDVEAREA(ht_index) .* Qinf .*(x_025_ht' - cog_x);

                N_cg_ht_Da = OUTP.vecCDPDIST(ht_index) .* cos(SURF.vecDVEPITCH(ht_index)).* cos(SURF.vecDVEYAW(ht_index)) .* SURF.vecDVEAREA(ht_index) .* Qinf .* -y_025_ht';
                N_cg_ht_Db = OUTP.vecCDPDIST(ht_index) .* cos(SURF.vecDVEPITCH(ht_index)).* cos(SURF.vecDVEYAW(ht_index)) .* SURF.vecDVEAREA(ht_index) .* Qinf .*(x_025_ht' - cog_x);

                CNcg(i) = CNcg(i) + sum(N_cg_ht_Da + N_cg_ht_Db) / (Qinf * ref_area * ref_span); 
            end
        end
    end
    
    if SETTINGS.stabSelection(2)
        % Clb
        p_2 = polyfit(dalpha,Clcg,1);
        STABILITY.Clb = p_2(1)*180/pi;
        STABILITY.Clcg_vec = Clcg;
    end 
    if SETTINGS.stabSelection(3) 
        % CNb
        p_3 = polyfit(dalpha,CNcg,1);
        STABILITY.CNb = p_3(1)*180/pi;
        STABILITY.CNcg_vec = CNcg;
    end
end