function [hFig2] = fcnPLOTBODY_2(verbose, valNELE, matDVE, matVLST, matCENTER, sym)
% New plot body function, plots to the next avaliable figure NOT just
% figure 3
% 
% Author: Tom Ryan

hFig2 = figure(); % This is the only line changed 
% clf(3)

patch('Faces',matDVE,'Vertices',matVLST,'FaceColor',[255 90 90]./255)
if sym == true
    patch('Faces',matDVE,'Vertices',[matVLST(:,1) matVLST(:,2).*-1 matVLST(:,3)],'FaceColor',[255 90 90]./255)
end
hold on


% alpha(0.5)

if verbose == 1
    for ii = 1:valNELE
        str = sprintf('%d',ii);
        text(matCENTER(ii,1),matCENTER(ii,2),matCENTER(ii,3),str,'Color','k','FontSize',15);
    end
    
    for ii = 1:length(matVLST(:,1))
        str = sprintf('%d',ii);
        text(matVLST(ii,1),matVLST(ii,2),matVLST(ii,3),str,'Color','g','FontSize',15);
    end
    
end

hold off


axis equal
axis tight
set(gcf,'color','w');

% With axes
box on
grid on
xlabel('X-Dir','FontSize',15);
ylabel('Y-Dir','FontSize',15);
zlabel('Z-Dir','FontSize',15);

% % % Clear background
% set(gca,'xtick',[]);
% set(gca,'ytick',[]);
% set(gca,'ztick',[]);
view(3)