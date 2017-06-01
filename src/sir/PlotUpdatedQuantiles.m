function [hFigures] = PlotUpdatedQuantiles(ForecastStruct,...
    ForecastProposal,ForecastTruth,HcProposal,Weights,CurrentTime,FontSize)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% Plot updated forecasts by weights

ColorMap = parula; ColorMap = ColorMap(1:end,:);
ProposalForecasts = ForecastProposal.data;

% We need to plot the lower weight ones first, or else nothing shows up
[SortedWeights,I] = sort(Weights);

%hFigures(1) = figure('Position', [100, 100, 1000, 800]);

hold on;
% Plot prior
for i = 1:length(HcProposal)
    h1 = plot(ForecastStruct.time, ForecastStruct.data(i,:),...
        'color',[0.5 0.5 0.5]);
end

% Plot posterior by weight
for i = 1:length(HcProposal)
    ColorIndex = min(round(SortedWeights(i)/max(SortedWeights)*...
        length(ColorMap))+1,length(ColorMap));
    ColorIndex = max(ColorIndex,1);
    Color = ColorMap(ColorIndex,:);
    h2 = plot(ForecastStruct.time, ProposalForecasts(I(i),:),'color',...
        Color,'LineWidth',2);
end
caxis([0 max(SortedWeights)]);

h3 = plot(ForecastStruct.time,ForecastTruth.data,...
     'r','LineWidth',3);

xlabel('Time (days)','FontSize',FontSize);
ylabel(['Forecasted: ' ForecastStruct.name],'FontSize',FontSize);
hlegend = legend([h1,h2,h3],'Prior','Posterior','Reference');
set(hlegend,'FontSize',FontSize);
title({'Updated Forecasts Using ', [num2str(CurrentTime) ...
    ' days of data']},'FontSize',FontSize);
hcolor = colorbar;  ylabel(hcolor,'Weight','FontSize',FontSize-4);
set(gcf,'color','w'); axis square;
set(gca,'FontSize',FontSize-4); axis square;

% Format axis labels
YMin = min(ForecastStruct.data(:));
YMax = max(ForecastStruct.data(:));
ylim([YMin YMax]); set(gcf,'color','w'); axis square; 
ax = gca; 
ax.XTick = roundn(linspace(ForecastStruct.time(1),...
    ForecastStruct.time(end),5),2);
set(gca,'XTickLabel',sprintf('%3.0f\n',ax.XTick))

PriorForecasts = ForecastStruct.data;
ResampledForecasts = ForecastProposal.data;

[ PriorQuantiles,PosteriorQuantiles ] = ComputeQuantiles(PriorForecasts,...
    ResampledForecasts );

hFigures(2) = figure('Position', [100, 100, 1000, 800]);
hold on;
h0 = plot(ForecastStruct.time,PriorQuantiles','color',[0.5 0.5 0.5],...
    'LineWidth',3);
h2 = plot(ForecastStruct.time,ForecastTruth.data,...
     'r','LineWidth',3);
h3 = plot(ForecastStruct.time,PosteriorQuantiles,'b--','LineWidth',3);
legend([h0(1), h2(1),h3(1)],'Prior','Reference','SIR-CFCA');

xlabel('Time (days)');ylabel(['Forecasted: ' ForecastStruct.name]);
title({'Updated Quantiles using ', [num2str(CurrentTime) ...
    ' days of data']},'FontSize',FontSize);
set(gca,'FontSize',FontSize-4);

% Format axis labels
YMin = min([PriorQuantiles(:);PosteriorQuantiles(:)]);
YMax = max([PriorQuantiles(:);PosteriorQuantiles(:)]);
ylim([YMin YMax]); set(gcf,'color','w'); axis square; 
ax = gca; 
ax.XTick = roundn(linspace(ForecastStruct.time(1),...
    ForecastStruct.time(end),5),2);
set(gca,'XTickLabel',sprintf('%3.0f\n',ax.XTick))

end

