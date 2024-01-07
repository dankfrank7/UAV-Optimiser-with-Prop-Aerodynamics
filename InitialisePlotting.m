function InitialisePlotting()
myred           = [216 30 49]/255;
myblue          = [27 99 157]/255;
myblack         = [0 0 0]/255;
mygreen         = [0 128 0]/255;
mycyan          = [2 169 226]/255;
myyellow        = [251 194 13]/255;
mygray          = [89 89 89]/255;
set(groot,'defaultAxesColorOrder',[myblack;myblue;myred;mygreen;myyellow;mycyan;mygray]);
%% define some general plot parameters
alw             = 1;                        % AxesLineWidth
fsz             = 14;                       % Fontsize
lw              = 2.5;                        % LineWidth
msz             = 12;                       % MarkerSize

set(0,'defaultLineLineWidth',lw);           % set the default line width to lw
set(0,'defaultLineMarkerSize',msz);         % set the default line marker size to msz
set(0,'defaultLineLineWidth',lw);           % set the default line width to lw
set(0,'defaultLineMarkerSize',msz);         % set the default line marker size to msz
set(0,'defaultAxesFontSize',fsz)             % set the default axes font size
set(0,'defaultLegendFontSize',18)
return