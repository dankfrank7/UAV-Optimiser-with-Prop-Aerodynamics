function thicknessAtMAC = fcnAIRFOIL_THICKNESS(wingAirfoil, ref_cmac)
% This function looks up the airfoil coordinate .dat file and estimates the
% thickness at the MAC
% 
% Author: Tom Ryan
% 7/09/2021

    % Load airfoil coordinates 
    
    filename = strcat('airfoils\',wingAirfoil, '.mat');
    load(filename,'coord')
    
    % Seperate into two halfs 
    top_index = coord(:,2) > 0; 
    top = coord(top_index,:);
    
    bot_index = coord(:,2) < 0;
    bot = coord(bot_index,:);
    
    new_x = linspace(0,1,100);
    
    % Redistribute coordinate points 
    top_new = interp1(top(:,1),top(:,2),new_x,'spline');
    bot_new = interp1(bot(:,1),bot(:,2),new_x,'spline');
    
    thick_vec = top_new + abs(bot_new);
    
    thicknessAtMAC = max(thick_vec)*ref_cmac;
end