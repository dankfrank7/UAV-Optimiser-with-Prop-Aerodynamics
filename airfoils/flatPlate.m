% This script generates the .mat file for a flate plate airfoil to be used
% in the propeller validation file 
% Author: Tom Ryan 
% 24/08/2021

% Start with a clear workspace
clear; clc;
close all 

airfoil_name = 'flateplate';
coord = table2array(readtable(strcat(airfoil_name,'.dat')));
Re_range = [10, 100, 1000]*1e3;
Alpha_range = -20:.25:30;
pol = zeros(length(Alpha_range),9,length(Re_range));
for i = 1:length(Re_range)
    pol(:,1,i) = Alpha_range';                              % 1st col is a
    pol(:,8,i) = ones(length(Alpha_range),1)*Re_range(i);   % 8th col is Re
    pol(:,6,i) = ones(length(Alpha_range),1);               % 6th col is Top_Xtr - not sure if this is necessary 
end
    
% Save file 
saveName = strcat(airfoil_name,'.mat');
save(saveName,'airfoil_name','coord','Re_range','pol')
