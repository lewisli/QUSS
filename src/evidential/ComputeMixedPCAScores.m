function [ score ] = ComputeMixedPCAScores( HistoricalStruct,EigenTolerance )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Compute d_f for noiseless
NummRealizations = size(HistoricalStruct.data,1);
NumTimeSteps = size(HistoricalStruct.data,2);
NumHistoricalResponses = size(HistoricalStruct.data,3);
PlotLevelPCA = 0;
MinEigenValues = 2;


histPCA = ComputeHarmonicScores(HistoricalStruct,PlotLevelPCA);


% Get number of required harmonics required for forecasts
nHarmHist = 0;
for i = 1:NumHistoricalResponses
    responseNumHarms = GetNumHarmonics(histPCA{i},...
        MinEigenValues,EigenTolerance);
    nHarmHist = max(nHarmHist,responseNumHarms);
end

% Generate matrix of historical harmonic scores
Df = zeros(NummRealizations,nHarmHist*NumHistoricalResponses);

% Iterate over each historical response (ex: P1, P2)
for i = 1:NumHistoricalResponses
    harmscrhist=histPCA{i}.harmscr;
    
    % Need to re-arrange harmonic scores into Df such that the first
    % eigenvalues are placed in first
    for j = 1:nHarmHist
        Df(:,(i-1)*nHarmHist + j) = harmscrhist(:,j);
    end
end


% Run PCA On Df
% Poor naming choice by FDA package... we need to remove it from the path
% in order to run PCA
rmpath('../../common//fda_matlab');
DfStar = Df;
[~,score,~] = pca(DfStar);



end

