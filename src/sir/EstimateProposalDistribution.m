function [ ProposalParameterVal,ProposalPDFs] = EstimateProposalDistribution(...
    PriorData,ObservedData,PriorParameters)
%ESTIMATEPROPOSALDISTRIBUTION Estimate the proposal distribution using
%currently observed data
%   Detailed explanation goes here

addpath('../src/thirdparty/likelihood_continuous/');

% Compute FPCA on prior data + observed data
data_FPCA = ComputeHarmonicScores(PriorData,ObservedData,0);

% ComputeHarmonicScore stacks the observed data as the last row of the
% resulting scores
obs_realization = size(PriorData.data,1) + 1;

% Perform Mixed PCA
eigentolerance = 0.9;
rmpath('../src/thirdparty/fda_matlab/');
[mpca_scores, mpca_obs] = MixedPCA(data_FPCA,obs_realization,...
    eigentolerance);

NumParameters = size(PriorParameters,2);
NumPDFPoints = 100;
ProposalParameterVal = zeros(NumPDFPoints,NumParameters);
ProposalPDFs = zeros(NumPDFPoints,NumParameters);

% Stack mpca_scores and mpca_obs to compute probability
mpca_stacked = [mpca_scores; mpca_obs];

for i = 1:NumParameters
    ParameterValues = PriorParameters(:,i);
    
    % Kernel Density Estimation
    [ProposalParameterVal(:,i),ProposalPDFs(:,i)] = ...
        UpdateProbabilityContinuous(mpca_stacked,ParameterValues);
end

end

