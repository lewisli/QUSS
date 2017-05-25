function [ C_f ] = EstimateFunctionalErrorCovariance( HistoricalStruct, EigenTolerance,epsilon)
%EstimateFunctionalErrorCovariance Uses a Monte Carlo approach to compute
%   Detailed explanation goes here
addpath('../cfca');
% Step 1: Extract

NummRealizations = size(HistoricalStruct.data,1);
NumTimeSteps = size(HistoricalStruct.data,2);
NumHistoricalResponses = size(HistoricalStruct.data,3);

histPCA = ComputeHarmonicScores(HistoricalStruct,0);
[noiseless_scores, ~] = MixedPCA(histPCA,0,EigenTolerance);

% 
df_error = zeros(size(noiseless_scores));

for i = 1:NummRealizations
     fprintf('Working on realization %d\n',i)
     NoisyStruct = HistoricalStruct;
     
     NoisyStruct.data(i,:,:) = HistoricalStruct.data(i,:,:) + ...
         randn(1,NumTimeSteps,NumHistoricalResponses)*epsilon;
     
     noisyPCA = ComputeHarmonicScores(NoisyStruct,0);
     
     [noisy_score, ~] = MixedPCA(noisyPCA,0,size(noiseless_scores,2));
     
     df_error(i,:) = noisy_score(i,:) - noiseless_scores(i,:);

end
 
C_f = cov(df_error);

end

