% This script creates of the optimisation output data 
% 
% Author: Tom Ryan

clear; clc;
close all
InitialisePlotting();

% Workspace 
cd('X:\Users\Tom\Documents\Uni Work\Y5S1\Thesis\aMainFiles')
curr_dir = pwd;
cd './VAP3point5'
addpath(curr_dir)


%% Load optimisation data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%loadFilename = 'genmav_NoProp_45span_DRIES_RESULTS'
loadFilename = 'WingTips_GWS4x4_30span_DRIES_RESULTS'
% loadFilename = 'WingTips_noprop_30span_DRIES_RESULTS'
%loadFilename = 'Flyingwing_noprop_45span_DRIES_RESULTS'
% loadFilename = 'Flyingwing_noprop_30span_DRIES_RESULTS'
%loadFilename = 'FlyingWing_GWS4x4_30span_DRIES_RESULTS'
%loadFilename = 'FlyingWing_GWS5x43_45span_DRIES_RESULTS'
%loadFilename = 'Genmav_GWS5x43_45span_DRIES_RESULTS'
% loadFilename = 'Genmav_APC7x6_45span_DRIES_RESULTS'
%loadFilename = 'FlyingWing_APC7x6_45span_DRIES_RESULTS'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% load(strcat(curr_dir,'\OUTPUT\',loadFilename,'.mat'))
load(strcat(curr_dir,'\OUTPUT\RESULTS\',loadFilename,'.mat'))
clear pwd;
curr_dir = pwd;
% SETTINGS.wingtips = 0;
% SETTINGS.systemWeight = 0.22;

% Axes 
Xtic_vec = 7:1:14;
Ytic_vec = 0.1:0.01:0.18;
ax_vec = [Xtic_vec(1),Xtic_vec(end),Ytic_vec(1),Ytic_vec(end)];

MaxIt = max(conv_histCost(:,end)) % Manual override if run does not go to complete MaxIt

% Filter out models that do not meet constraints 
constraintIndex = conv_histConstraint(:,1:end-1) < 0; 
constraintIndex_row = prod(constraintIndex,2);

conv_histCost_filt          = conv_histCost(constraintIndex_row > 0,:);
conv_histConstraint_filt    = conv_histConstraint(constraintIndex_row > 0,:);
conv_histPos_filt           = conv_histPos(constraintIndex_row > 0,:);

% Extract stability derivatives 
if SETTINGS.AnalyseStability
	conv_histStability.Cma = conv_histConstraint_filt(:,end-7) + SETTINGS.Cma_upper;
    conv_histStability.Clb = conv_histConstraint_filt(:,end-5) + SETTINGS.Clb_upper;
    conv_histStability.Cnb = conv_histConstraint_filt(:,end-3) + SETTINGS.CNb_upper;
end



% Iteration plot 
figure(1) 
    scatter(reshape(-conv_histCost_filt(:,1),[],1),reshape(conv_histCost_filt(:,2),[],1),60,conv_histCost_filt(:,3),'filled')
    hold on
    grid on
    c1 = colorbar;
    c1.Label.String = 'Iteration';
    %scatter(-CostPlot(:,1),CostPlot(:,2),'rx')
    xlabel('Lift to Drag Ratio')
    ylabel('Mass [kg]')
    xticks(Xtic_vec)
    yticks(Ytic_vec)
    axis(ax_vec)
%     axis tight


% Plots for each variable
for i = 1:SETTINGS.nVar
    figure(i+1)
    [VAR_denorm, label] = fcnDENORM_VAR(conv_histPos(:,i), SETTINGS.variableNames, SETTINGS.varMin, SETTINGS.varMax, i);
    VAR_denorm = VAR_denorm(constraintIndex_row > 0);
    scatter(reshape(-conv_histCost_filt(:,1),[],1),reshape(conv_histCost_filt(:,2),[],1),60,VAR_denorm,'filled')
    c_temp = colorbar;
    c_temp.Label.String = label;
    grid on
    xlabel('Lift to Drag Ratio')
    ylabel('Mass [kg]')
    xticks(Xtic_vec)
    yticks(Ytic_vec)
    axis(ax_vec)
%     axis tight
end

% Stability derivative Plots 
label_vec = {'C_{M_\alpha}','C_{l_\beta}','C_{n_\beta}'};
fields = {'Cma','Clb','Cnb'};
for i = 1:3 
    figure()
    scatter(reshape(-conv_histCost_filt(:,1),[],1),reshape(conv_histCost_filt(:,2),[],1),60,conv_histStability.(fields{i}),'filled')
    c_temp = colorbar;
    c_temp.Label.String = label_vec{i};
    grid on 
    xlabel('Lift to Drag Ratio')
    ylabel('Mass [kg]')
    xticks(Xtic_vec)
    yticks(Ytic_vec)
    axis(ax_vec)
%     axis tight
end

%  AR plot for wingtips 
for i = 1:SETTINGS.nVar
    [conv_histPos_filt_denorm_full(:,i), ~] = fcnDENORM_VAR(conv_histPos_filt(:,i), SETTINGS.variableNames, SETTINGS.varMin, SETTINGS.varMax, i);
end
b_vec = conv_histPos_filt_denorm_full(:,9);
cr_vec = conv_histPos_filt_denorm_full(:,10);
taper_vec = conv_histPos_filt_denorm_full(:,12);
S_vec = b_vec.*(cr_vec.*(1+taper_vec)./2);
AR_vec = b_vec.^2./S_vec;
figure() 
    scatter(reshape(-conv_histCost_filt(:,1),[],1),reshape(conv_histCost_filt(:,2),[],1),60,AR_vec,'filled')
    hold on
    grid on
    c1 = colorbar;
    c1.Label.String = 'Aspect Ratio';
    caxis([1,4])
    xlabel('Lift to Drag Ratio')
    ylabel('Mass [kg]')
    xticks(Xtic_vec)
    yticks(Ytic_vec)
    axis(ax_vec)

% Extract specific geometries from the final iteration 

% final_cost = conv_histCost_filt(conv_histCost_filt(:,end) == (MaxIt), 1:2);
% final_pos = conv_histPos(conv_histPos(:,end) == (MaxIt), 1:nVar);
% [~, index] = min(final_cost);
% ind_maxLtoD   = index(1);
% ind_minWeight = index(2);
% ind_bestComp  = fcnBESTCOMP(final_cost);
% 
% VARS_maxLtoD    = final_pos(ind_maxLtoD,:);
% VARS_minWeight  = final_pos(ind_minWeight,:);
% VARS_bestComp   = final_pos(ind_bestComp,:);
% 
% hFig(1) = fcnPLOT_GEOM(VARS_maxLtoD, SETTINGS);
% hFig(1).Name = 'Max LtoD';
% hFig(2) = fcnPLOT_GEOM(VARS_minWeight, SETTINGS);
% hFig(2).Name = 'Min Weight';
% hFig(3) = fcnPLOT_GEOM(VARS_bestComp, SETTINGS);
% hFig(3).Name = 'Best Compromise';

%% Correlation Plot 
lastIts = 5;
%selectVars = [3,7,9,10,11,12]; % WingTips 30 cm
%selectVars = [4,7,9,10,11]; % Flying Wing 45 prop
selectVars = [2, 4,6,8,9,10]; % Genmav 7 inprop
%selectVars = [3, 5, 7, 9,10,11]; % Genmav noprop
%selectVars = [3,4,7,9,10,11]; % Flying Wing 30cm
%corrNames = SETTINGS.variableNames(selectVars);
for i = 1:SETTINGS.nVar
    [conv_histPos_filt_denorm(:,i), ~] = fcnDENORM_VAR(conv_histPos_filt((end-nPop*lastIts):end,i), SETTINGS.variableNames, SETTINGS.varMin, SETTINGS.varMax, i);
end
%corrNames = {'VT Vol','Dihed [deg]','Span [m]','Root Chord [m]','Sweep [deg]','Taper','L/D','Mass [kg]','Cma','Clb','Cnb'}; % Wingtips 30 cm
%corrNames = {'VT Vol','Dihed [deg]','AR','Root Chord [m]','Sweep [deg]','L/D','Mass [kg]','Cma','Clb','Cnb'};% Flying Wing 45 prop
%corrNames = {'VT Vol','VT2 Vol','Dihed [deg]','AR','Root Chord [m]','Sweep [deg]','L/D','Mass [kg]','Cma','Clb','Cnb'};% Flying Wing 45 prop
corrNames = {'VT Vol','HT Vol','Dihed [deg]','AR','Root Chord [m]','Sweep [deg]','L/D','Mass [kg]','Cma','Clb','Cnb'}; % Genmav 7in
corrData  = [conv_histPos_filt_denorm(:,selectVars), -conv_histCost_filt((end-nPop*lastIts):end,1), conv_histCost_filt((end-nPop*lastIts):end,2), conv_histStability.Cma((end-nPop*lastIts):end), conv_histStability.Clb((end-nPop*lastIts):end), conv_histStability.Cnb((end-nPop*lastIts):end)];
figure
[R,PValue, h_corr] = corrplot2(corrData, 'varNames', corrNames);
set(gcf,'color','w');
set(findall(gcf,'-property','FontSize'),'FontSize',16)
% corr(corrData)

%% Extra Plots 
% figure % Vt vs Cma
%     x = conv_histPos_filt(:,3);
%     y = conv_histPos_filt(:,4);
%     z = conv_histStability.Cma;
%     x(isnan(x)) = [];
%     y(isnan(y)) = [];
%     z(isnan(z)) = [];
%     f = fit([x,y], z, 'lowess');
%     box on 
%     grid on 
%     hold on 
%     plot(f, [x,y], z)
%     scatter3(x, y, z, 50, z, 'filled')
%     alpha 0.5
%     xlabel('V_t Top')
%     ylabel('V_t Bot')
%     zlabel('C_{m_\alpha}')
%     xlim([0,0.8])
%     ylim([0,0.4])

% figure % Sweep vs AR vs Cma 
%     x = conv_histPos_filt(:,11);
%     y = conv_histPos_filt(:,9);
%     z = conv_histStability.Cma;
%     x(isnan(x)) = [];
%     y(isnan(y)) = [];
%     z(isnan(z)) = [];
%     f = fit([x,y], z, 'lowess');
%     box on 
%     grid on 
%     hold on 
%     plot(f, [x,y], z)
%     scatter3(x, y, z, 50, z, 'filled')
%     alpha 0.5
%     xlabel('Sweep [deg]')
%     ylabel('AR ')
%     zlabel('C_{m_\alpha}')
%     xlim([0,0.8])
%     ylim([0.1,0.3])

% figure % Vt vs Clb
    % x = conv_histPos_filt(:,3);
    % y = conv_histPos_filt(:,4);
    % z = conv_histStability.Clb;
    % x(isnan(x)) = [];
    % y(isnan(y)) = [];
    % z(isnan(z)) = [];
    % f = fit([x,y], z, 'lowess');
    % box on 
    % grid on 
    % hold on 
    % plot(f, [x,y], z)
    % scatter3(x, y, z, 50, z, 'filled')
    % alpha 0.5
    % xlabel('V_t Top')
    % ylabel('V_t Bot')
    % zlabel('C_{l_\beta}')
    % xlim([0,0.8])
    % ylim([0,0.4])

% figure% Vt vs Clb
    % x = conv_histPos_filt(:,3);
    % y = conv_histPos_filt(:,4);
    % z = conv_histStability.Cnb;
    % x(isnan(x)) = [];
    % y(isnan(y)) = [];
    % z(isnan(z)) = [];
    % f = fit([x,y], z, 'lowess');
    % box on 
    % grid on 
    % hold on 
    % plot(f, [x,y], z)
    % scatter3(x, y, z, 50, z, 'filled')
    % alpha 0.5
    % xlabel('V_t Top')
    % ylabel('V_t Bot')
    % zlabel('C_{n_\beta}')
    % xlim([0,0.8])
    % ylim([0,0.4])

toSave = 1; % This part won't work on Linux 
if toSave 
    folderName = strcat('X:\Users\Tom\Documents\Uni Work\Y5S1\Thesis\Optimisation Results\Results\',loadFilename);
    mkdir(folderName);

    for i = (1:SETTINGS.nVar + 1)
        if i == 1
            saveas(figure(i), strcat(folderName, '\Iterations.png'))
        else 
            saveas(figure(i), strcat(folderName, '\', SETTINGS.variableNames{i-1}, '.png'))
        end
    end
    ii = 1;
%     for i = (SETTINGS.nVar+2):(SETTINGS.nVar+4)
%         saveas(figure(i), strcat(folderName, '\', fields{ii}),'.png')
%         ii = ii + 1;
%     end

    saveas(hFig(1), strcat(folderName, '\MaxLtoD.png'))
    saveas(hFig(2), strcat(folderName, '\MinWeight.png'))
    saveas(hFig(3), strcat(folderName, '\BestCompromise.png'))
end

% Revert directory
% cd(curr_dir)
