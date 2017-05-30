function [FigHandle] = PlotModelsByWeight( Hc1,Dc1,Hc2,Dc2,dobs_c1,WNorm,...
    FontSize)
%PlotModelsByWeight Plots scatter plot of Hc vs Dc for posterior samples
%colored by weights given WNorm.
%
% Plots scatter plot of Hc vs Dc for posterior samples colored by weights given
% in WNorm.
%
% Inputs:
%   Hc1: Canonical coeffcients for prior historical
%   Dc1: Canonical coeffcients for prior forecasts
%   Hc2: Canonical coeffcients for posterior historical
%   Hc1: Canonical coeffcients for posterior forecasts
%   dobs_c1: Canonical coeffcients for reference observed historical
%   WNorm: Normalized weights of posterior models
%
% Author: Lewis Li (lewisli@stanford.edu)
% Date: March 4th 2016

FigHandle = figure('Position', [100, 100, 1000, 800]);


% Check if WNorm is unique
if length(unique(WNorm)) > 1
    scatter(Dc1(:,1),Hc1(:,1),100,[0.5 0.5 0.5],'filled'); hold on;
    scatter(Dc2(:,1),Hc2(:,1),100,WNorm,'filled');
    h = colorbar;    
    ylabel(h, 'Weights','FontSize',FontSize);
else
    h1 = scatter(Dc1(:,1),Hc1(:,1),100,'k','filled'); hold on;
    h2 = scatter(Dc2(:,1),Hc2(:,1),100,'r','filled');
    hlegend = legend([h1(1), h2(1)], 'Prior', 'Resampled Proposal');
    set(hlegend,'Location','southeast');
end


h = plot([dobs_c1(1) dobs_c1(1)],[min(Hc1(:,1)) ...
    max(Hc1(:,1))],'b','LineWidth',2);
xlabel('d_c','FontSize',FontSize);
ylabel('h_c','FontSize',FontSize);
set(gcf,'color','w');
set(gca,'FontSize',FontSize-4);
grid on;
axis tight; axis square;
tt1=text(dobs_c1(1)+0.05,min(Hc1(:,1))+0.5,'d_{obs}',...
    'Fontweight','b','FontSize',FontSize);
set(tt1,'color','b','FontSize',FontSize);

end

