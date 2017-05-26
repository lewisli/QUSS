function [ PriorQuantiles,PosteriorQuantiles ] = ComputeQuantiles(Prior,...
    Posterior)
%COMPUTEQUANTILES Computes P10-P50-P90 for prior and posterior samples
%
% Inputs:
%   Prior: (NPriorReals x NTimeSteps) prior runs
%   Posterior: (NPostReals x NTimeSteps) posterior runs
%
% Outputs:
%   PriorQuantiles: (3 x NTimeSteps) prior quantiles
%   PosteriorQuantiles: (3 x NTimeSteps) posterior quantiles
%
% Author: Lewis Li (lewisli@stanford.edu)
% Date:    Feburary 7th 2016 

% P10-P50-P90
quantiles = [.1,.5,.9];

PriorQuantiles = quantile(Prior,quantiles);
PosteriorQuantiles = quantile(Posterior,quantiles);

end

