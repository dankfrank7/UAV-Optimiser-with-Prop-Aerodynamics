function hFig = fcnPLOT_OUTPUT(conv_histCost, conv_histPos, conv_histConstraint ,CostPlot, SETTINGS, MaxIt, saveFilename, toSave)
% This is the function version of ascript creates of the optimisation output data 
% 
% Author: Tom Ryan
    InitialisePlotting();
    
    % Axes 
    tic_vec = 8:2:16;
    ax_vec = [tic_vec(1),tic_vec(end),0.26,0.40];


    % Filter out models that do not meet constraints 
    conv_histCost_filt = conv_histCost;
    [~, cols] = size(conv_histCost);
    [rows, ~] = size(conv_histConstraint);
    constraintIndex = conv_histConstraint(:,1:end-1) < 0; 
    for r = 1:rows
        if prod(constraintIndex(r,:)) == 0
            conv_histCost_filt(r,:) = nan(1,cols);
        end
    end


    % Iteration plot 
    figure(1) 
        scatter(reshape(-conv_histCost(:,1),[],1),reshape(conv_histCost(:,2),[],1),60,conv_histCost(:,3),'filled')
        hold on
        grid on
        c1 = colorbar;
        c1.Label.String = 'Iteration';
        scatter(-CostPlot(:,1),CostPlot(:,2),'rx')
        xlabel('Lift to Drag Ratio')
        ylabel('Mass [kg]')


    % Plots for each variable
    for i = 1:SETTINGS.nVar
        figure(i+1)
        [VAR_denorm, label] = fcnDENORM_VAR(conv_histPos(:,i), SETTINGS.variableNames, SETTINGS.varMin, SETTINGS.varMax, i);
        scatter(reshape(-conv_histCost(:,1),[],1),reshape(conv_histCost(:,2),[],1),60,VAR_denorm,'filled')
        c_temp = colorbar;
        c_temp.Label.String = label;
        grid on
        xlabel('Lift to Drag Ratio')
        ylabel('Mass [kg]')
    end

    % Extract specific geometries from the final iteration 
    final_cost = conv_histCost(conv_histCost(:,end) == (MaxIt+1), 1:2);
    final_pos = conv_histPos(conv_histPos(:,end) == (MaxIt + 1), 1:SETTINGS.nVar);
    [~, index] = min(final_cost);
    ind_maxLtoD   = index(1);
    ind_minWeight = index(2);
    ind_bestComp  = fcnBESTCOMP(final_cost);

    VARS_maxLtoD    = final_pos(ind_maxLtoD,:);
    VARS_minWeight  = final_pos(ind_minWeight,:);
    VARS_bestComp   = final_pos(ind_bestComp,:);

    hFig(1) = fcnPLOT_GEOM(VARS_maxLtoD, SETTINGS);
    hFig(1).Name = 'Max LtoD';
    hFig(2) = fcnPLOT_GEOM(VARS_minWeight, SETTINGS);
    hFig(2).Name = 'Min Weight';
    hFig(3) = fcnPLOT_GEOM(VARS_bestComp, SETTINGS);
    hFig(3).Name = 'Best Compromise';

    if toSave 
        folderName = strcat('X:\Users\Tom\Documents\Uni Work\Y5S1\Thesis\Optimisation Results\Initial\',saveFilename);
        mkdir(folderName);
        
        for i = (1:SETTINGS.nVar + 1)
            if i == 1
                saveas(figure(i), strcat(folderName, '\Iterations.png'))
            else 
                saveas(figure(i), strcat(folderName, '\', SETTINGS.variableNames{i-1}, '.png'))
            end
        end
        
        saveas(hFig(1), strcat(folderName, '\MaxLtoD.png'))
        saveas(hFig(2), strcat(folderName, '\MinWeight.png'))
        saveas(hFig(3), strcat(folderName, '\BestCompromise.png'))
    end
end