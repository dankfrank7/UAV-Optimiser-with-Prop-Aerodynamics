function ind_bestComp = fcnBESTCOMP(final_cost)
% Functiond finds the distance of every cost position from [-maxL/D, 0] and
% outputs a vector of distances. This is then used to find the optimal
% compromise design. 
% 
% Author: Tom Ryan
    
    [maxLD, ~] = min(final_cost(:,1));
    
    % Define datum point
    x0 = 1;
    y0 = 0;

    % Distance of each point in the final cost vec to the datum point
    dist_vec = sqrt((final_cost(:,1)./maxLD-x0).^2 + (final_cost(:,2)-y0).^2);
    
    % Best compromise is the closest point to datum
    [~, ind_bestComp] = min(dist_vec);
    
%     scatter(-final_cost(:,1),final_cost(:,2),60,dist_vec,'filled')
%     dist_vec = colorbar;
    
end