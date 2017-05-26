function [ ] = PlotPosteriorSamplesAndQuantiles( PredictionStruct,...
    TruthRealization, SampledPosteriorRealizations, PriorQuantiles,...
    PosteriorQuantiles,SavePath)
%PLOTPOSTERIORSAMPLESANDQUANTILES Plots posterior realizations and
%quantiles
%
% Inputs:
%   PredictionStruct: Struct containing prior prediction variable
%   TruthRealization: Index of realization that is used as "Truth"
%   SampledPosteriorRealizations: Posterior realizations in time domain
%   PriorQuantiles: Prior prediction quantiles
%   PosteriorQuantiles: Posterior prediction quantiles
%   SavePath: Path to save figures [optional]
%
% Outputs:
%   None:
%
% Author: Lewis Li (lewisli@stanford.edu)
% Date: March 4th 2016

if (nargin < 6)
    SaveOn = false;
else
    SaveOn = true;
end

% Plot quantiles
figure('Position', [100, 100, 1049, 895]);
subplot(121);
hold on;
h0 = plot(PredictionStruct.time,PriorQuantiles','color',[0.5 0.5 0.5],...
    'LineWidth',3);
h2 = plot(PredictionStruct.time,PredictionStruct.data(TruthRealization,:),...
    'r','LineWidth',3);
h3 = plot(PredictionStruct.time,PosteriorQuantiles,'b--','LineWidth',3);

legend([h0(1), h3(1),h2(1)],'Prior','Posterior','Reference');
xlabel('t(days)');ylabel(['Forecasted: ' PredictionStruct.name]);axis square;
title('Quantiles');
set(gca,'FontSize',18);
axis tight;


% Plot posterior samples
subplot(122);
hold on;

%% Plot prior
for i = 1:size(PredictionStruct.data,2)
    h1=plot(PredictionStruct.time,PredictionStruct.data(i,:),...
        'color',[0.5 0.5 0.5]);
end


NumPosteriorSamples = size(SampledPosteriorRealizations,1);
for i = 1:NumPosteriorSamples
    h2=plot(PredictionStruct.time,SampledPosteriorRealizations(i,:),...
        'color','b');
end
h3=plot(PredictionStruct.time,PredictionStruct.data(TruthRealization,:),...
    'color','r','LineWidth',3); axis square;
set(gca,'FontSize',18);

legend([h1, h2,h3],'Prior','Posterior','Reference');
xlabel('t(days)');ylabel(['Forecasted: ' PredictionStruct.name]);axis square;
title('Posterior Samples');
axis tight;

set(gcf,'color','w');

if SaveOn == true
    export_fig([SavePath '_Posterior'], '-png','-m3');
end


end

