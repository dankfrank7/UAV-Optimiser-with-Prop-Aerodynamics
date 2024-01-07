clc
clear

% From a dat file
airfoil_name = 'S5010';
coord = dlmread([airfoil_name, '.dat'],'',1,0);

Reynolds = [10:10:100,125:25:250,300:50:500,600:100:1000]*1e3;

Alpha = [-20:0.25:30];

VAP3_airfoil_gen(airfoil_name, coord, Reynolds, Alpha)

