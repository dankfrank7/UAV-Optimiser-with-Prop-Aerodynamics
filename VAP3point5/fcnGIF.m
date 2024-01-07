function [hFig3] = fcnGIF(valTIMESTEP, FLAG, SURF, VISC, WAKE, COND, INPU, case_num)

hFig = fcnPLOTPKG(valTIMESTEP, FLAG, SURF, VISC, WAKE, COND, INPU);

% These two lines added by tom to force white background and make figure
% large - 18/10/2021
%hFig.Position = [1200 300 1000 1000]; 
set(gcf,'color','w');

% view([33 22])
view([-45 20])
hFig3 = figure(3);
frame = getframe(hFig3);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);

% Write to the GIF File

gif_str = ['GIF/output_',num2str(case_num),'.gif'];

if valTIMESTEP == 1
    imwrite(imind,cm, gif_str,'gif', 'Loopcount',inf);
else
    imwrite(imind,cm, gif_str,'gif','WriteMode','append', 'DelayTime', 0.2);
end