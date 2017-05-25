function [ h ] = PlotLowDimModels(D,H,DObs,Type,FontSize)
%PlotLowDimModels Plot a low dimensional projection of the reservoir models
%   Plot the models in some low dimension space according to data and
%   forecast
% Inputs:
%   D: Data
%   H: Forecast
%   DObs: Observed data
%   Type: c or f for canonical space or functional space
% Return:
%   h: handle to figure
%
% Author: Lewis Li (lewisli@stanford.edu)
% Date: May 2nd 2016

if (nargin < 4)
    FontSize=12;
end

ScatterSize=100;
ObservedLineThickness=3;
NumPlots=3;
MaxPlotPerRow=3;
NumDimensions=size(H,2);

for i=1:NumDimensions

    
    if mod(i,NumPlots)==1
        h=figure('Units', 'normalized', 'Position', [0,0,1,1]);
    end
    
    subplot(NumPlots/MaxPlotPerRow,MaxPlotPerRow,mod(i-1,NumPlots)+1);
    hold on;
    scatter(D(:,i),H(:,i),ScatterSize,'filled')
    plot([DObs(i),DObs(i)],[min(H(:,i)),max(H(:,i))],'r-',...
        'LineWidth',ObservedLineThickness);
    
    text(DObs(i) + abs(DObs(i))*0.25,min(H(:,i)) + ...
        abs(min(H(:,i)))*0.25,'d_{obs}','Fontweight','b','FontSize',FontSize);
    xlabel(['d_',num2str(i),'^' Type],'FontSize',FontSize);
    ylabel(['h_',num2str(i),'^' Type],'FontSize',FontSize);
    rho = corrcoef(D(:,i),H(:,i));
    title(['\rho = ' num2str(rho(1,2))],'FontSize',FontSize);
   
    set(gca,'FontSize',FontSize);
    axis square; axis tight;
    set(gcf,'color','w');
    
end


end

