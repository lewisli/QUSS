function [ ] = PlotPosteriorSamplesAndQuantiles( ForecastStruct,...
    TruthRealization, SampledPosteriorRealizations, PriorQuantiles,...
    PosteriorQuantiles,SavePath)
%PlotPosteriorSamplesAndQuantiles Plots CFCA posterior realizations and
%quantiles
%
% Author: Lewis Li (lewisli@stanford.edu)
% Date: March 4th 2016

if (nargin < 6)
    SaveOn = false;
else
    SaveOn = true;
end



% Plot quantiles
FigHandle = figure('Position', [100, 100, 1049, 895]);
subplot(121);
hold on;
h0 = plot(ForecastStruct.time,PriorQuantiles','color',[0.5 0.5 0.5],...
    'LineWidth',3);
h2 = plot(ForecastStruct.time,ForecastStruct.data(TruthRealization,:),...
    'r','LineWidth',3);
h3 = plot(ForecastStruct.time,PosteriorQuantiles,'b--','LineWidth',3);

legend([h0(1), h3(1),h2(1)],'Prior','Posterior','Reference');
xlabel('t(days)');ylabel(['Forecasted: ' ForecastStruct.name]);axis square;
title('Quantiles');
set(gca,'FontSize',18);
axis tight;


% Plot posterior samples
subplot(122);
hold on;

%% Plot prior
for i = 1:size(ForecastStruct.data,2)
    h1=plot(ForecastStruct.time,ForecastStruct.data(i,:),...
        'color',[0.5 0.5 0.5]);
end


NumPosteriorSamples = size(SampledPosteriorRealizations,1);
for i = 1:NumPosteriorSamples
    h2=plot(ForecastStruct.time,SampledPosteriorRealizations(i,:),...
        'color','b');
end
h3=plot(ForecastStruct.time,ForecastStruct.data(TruthRealization,:),...
    'color','r','LineWidth',3); axis square;
set(gca,'FontSize',18);

legend([h1, h2,h3],'Prior','Posterior','Reference');
xlabel('t(days)');ylabel(['Forecasted: ' ForecastStruct.name]);axis square;
title('Posterior Samples');
axis tight;

set(gcf,'color','w');

if SaveOn == true
    export_fig([SavePath '_Posterior'], '-png','-m3');
end


end

