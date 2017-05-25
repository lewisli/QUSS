function [ h ] = PlotLowDimScatter(D,H,DObs,TruthRealization,Type,FontSize)
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


h=figure('Units', 'normalized', 'Position', [0,0,1,1]);
for i=1:3
    subplot(1,3,i);
    hold on;
    scatter(D(:,i),H(:,i),ScatterSize,'filled')
    plot([DObs(i),DObs(i)],[min(H(:,i)),max(H(:,i))],'r-',...
        'LineWidth',ObservedLineThickness);
    text(DObs(i) + abs(DObs(i))*0.25,min(H(:,i)) + ...
        abs(min(H(:,i)))*0.25,'d_{obs}','Fontweight','b','FontSize',FontSize);
    xlabel(['d_',num2str(i),'^' Type],'FontSize',FontSize);
    ylabel(['h_',num2str(i),'^' Type],'FontSize',FontSize);
    coeff = corrcoef(D(:,i),H(:,i));

    
    if (TruthRealization~=0)
       plot([min(D(:,i)),max(D(:,i))],[(H(TruthRealization,i)),(H(TruthRealization,i))],'g-',...
        'LineWidth',ObservedLineThickness);
    end
    
    set(gca,'FontSize',FontSize);
    axis square; axis tight;
    set(gcf,'color','w');
    
end



% coeff = corrcoef(D(:,1),H(:,1));
% % title(['d^' Type '_1 vs h^' Type '_1. \rho = ' ...
% %     num2str(coeff(2)) ],'FontSize',FontSize);
% set(gca,'FontSize',FontSize);
% axis square; axis tight;
% 

% hold on;
% scatter(D(:,2),H(:,2),ScatterSize,'filled')
% plot([DObs(2),DObs(2)],[min(H(:,2)),max(H(:,2))],'r-',...
%     'LineWidth',ObservedLineThickness);
% text(DObs(2) + abs(DObs(2))*0.25,min(H(:,2)) + ...
%     abs(min(H(:,2)))*0.25,'d_{obs}','Fontweight','b','FontSize',FontSize);
% xlabel(['d_',num2str(2),'^' Type],'FontSize',FontSize);
% ylabel(['h_',num2str(2),'^' Type],'FontSize',FontSize);
% coeff = corrcoef(D(:,2),H(:,2));
% % title(['d^' Type '_2 vs h^' Type '_2. \rho = ' ...
% %     num2str(coeff(2)) ],'FontSize',FontSize);
% set(gca,'FontSize',FontSize);
% axis square; axis tight;
% 
% subplot(133);



end

