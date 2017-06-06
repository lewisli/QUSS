function [hFigures] = PlotUpdatedQuantiles(PriorPrediction,...
    ProposalPrediction,TruePrediction,HcProposal,Weights,CurrentTime,FontSize)
%PlotUpdatedQuantiles Plot updated quantiles
%
% Plots prior vs posterior prediction quantiles as well as the proposal
% predictions colored by weight.
%
% Inputs:
%   PriorPrediction: Struct containing prior prediction
%   ProposalPrediction: Struct containing proposal prediction
%   TruePrediction: Struct containing true prediction
%   HcProposal: Canonical coeffcients for proposal prediction
%   Weights: Normalized weights of proposal models
%   CurrentTime: Current time step (days)
%   FontSize: Font siez for plots
%
% Outputs:
%   hFigures: Figure handles
%
% Author: Lewis Li (lewisli@stanford.edu)
% Date: March 4th 2016

ColorMap = parula; ColorMap = ColorMap(1:end,:);
ProposalForecasts = ProposalPrediction.data;

% We need to plot the lower weight ones first, or else nothing shows up
[SortedWeights,I] = sort(Weights);

hold on;
% Plot prior
for i = 1:length(HcProposal)
    h1 = plot(PriorPrediction.time, PriorPrediction.data(i,:),...
        'color',[0.5 0.5 0.5]);
end

% Plot posterior by weight
for i = 1:length(HcProposal)
    ColorIndex = min(round(SortedWeights(i)/max(SortedWeights)*...
        length(ColorMap))+1,length(ColorMap));
    ColorIndex = max(ColorIndex,1);
    Color = ColorMap(ColorIndex,:);
    h2 = plot(PriorPrediction.time, ProposalForecasts(I(i),:),'color',...
        Color,'LineWidth',2);
end
caxis([0 max(SortedWeights)]);

h3 = plot(PriorPrediction.time,TruePrediction.data,...
     'r','LineWidth',3);

xlabel('Time (days)','FontSize',FontSize);
ylabel(['Forecasted: ' PriorPrediction.name],'FontSize',FontSize);
hlegend = legend([h1,h2,h3],'Prior','Posterior','Reference');
set(hlegend,'FontSize',FontSize);
title({'Updated Forecasts Using ', [num2str(CurrentTime) ...
    ' days of data']},'FontSize',FontSize);
hcolor = colorbar;  ylabel(hcolor,'Weight','FontSize',FontSize-4);
set(gcf,'color','w'); axis square;
set(gca,'FontSize',FontSize-4); axis square;

% Format axis labels
YMin = min(PriorPrediction.data(:));
YMax = max(PriorPrediction.data(:));
ylim([YMin YMax]); set(gcf,'color','w'); axis square; 
ax = gca; 
ax.XTick = roundn(linspace(PriorPrediction.time(1),...
    PriorPrediction.time(end),5),2);
set(gca,'XTickLabel',sprintf('%3.0f\n',ax.XTick))

PriorForecasts = PriorPrediction.data;
ResampledForecasts = ProposalPrediction.data;

[ PriorQuantiles,PosteriorQuantiles ] = ComputeQuantiles(PriorForecasts,...
    ResampledForecasts );

hFigures(2) = figure('Position', [100, 100, 1000, 800]);
hold on;
h0 = plot(PriorPrediction.time,PriorQuantiles','color',[0.5 0.5 0.5],...
    'LineWidth',3);
h2 = plot(PriorPrediction.time,TruePrediction.data,...
     'r','LineWidth',3);
h3 = plot(PriorPrediction.time,PosteriorQuantiles,'b--','LineWidth',3);
legend([h0(1), h2(1),h3(1)],'Prior','Reference','SIR-CFCA');

xlabel('Time (days)');ylabel(['Forecasted: ' PriorPrediction.name]);
title({'Updated Quantiles using ', [num2str(CurrentTime) ...
    ' days of data']},'FontSize',FontSize);
set(gca,'FontSize',FontSize-4);

% Format axis labels
YMin = min([PriorQuantiles(:);PosteriorQuantiles(:)]);
YMax = max([PriorQuantiles(:);PosteriorQuantiles(:)]);
ylim([YMin YMax]); set(gcf,'color','w'); axis square; 
ax = gca; 
ax.XTick = roundn(linspace(PriorPrediction.time(1),...
    PriorPrediction.time(end),5),2);
set(gca,'XTickLabel',sprintf('%3.0f\n',ax.XTick))

end

