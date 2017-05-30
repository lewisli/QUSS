function [ Hc2, Dc2 ] = ComputeCanonicalOnExternalBasis( HistoricalPrior, ...
    ForecastPrior,HistoricalProposal,ForecastProposal,EigenvalueTolerance,...
    HistoricalTruth,OutlierPercentile,NumHistEig,NumPredEig)
%ComputeCanonicalOnExternalBasis Computes canonical coeffcients of
%posterior samples using same eigenfunctions as prior to get them on the
%same support
%

% Author: Lewis Li
% Date: March 4th 2016

NumSamples = size(HistoricalProposal.data,1);
Dc2 = zeros(NumSamples,NumHistEig);
Hc2 = zeros(NumSamples,NumPredEig);,

h = waitbar(0,'Computing Proposal Hc and Dc...');
for i = 1:NumSamples
    
    % Add truth to proposal for calculating mu
    HistoricalPrior.data = cat(1,HistoricalPrior.data,...
        HistoricalProposal.data(i,:,:));
    ForecastPrior.data = [ForecastPrior.data; ForecastProposal.data(i,:)];
    
    [ ~, ~, Dc1,~,Hc1,~ ,~,~] = ...
        ComputeCFCAPosterior( HistoricalPrior, ForecastPrior, ...
        HistoricalTruth, EigenvalueTolerance,OutlierPercentile);
    
    if (size(Dc1(end,:),2) >=  NumHistEig)
        Dc2(i,:) = Dc1(end,1:NumHistEig);
    end
    
    if (size(Hc1(end,:),2) >=  NumPredEig)
        Hc2(i,:) = Hc1(end,1:NumPredEig);
    end
    
    close all;
    
    ForecastPrior.data(end,:)= [];
    HistoricalPrior.data(end,:,:)= [];
    
    display(['Computing score for realization ' num2str(i)]);
    waitbar(i  / NumSamples);
end
close(h);
end

