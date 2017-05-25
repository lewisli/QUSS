function [ PriorQuantiles,PosteriorQuantiles ] = ComputeQuantiles(Prior,...
    Posterior )
%ComputeQuantiles Compute P10P50P90
%   Computes P10P50P90 for prior and posterior samples
%
% Inputs:
%   Prior: prior runs
%   Posterior: posterior runs
%
% Outputs:
%   PriorQuantiles: prior quantiles
%   PosteriorQuantiles: posterior quantiles
%
% 

quantiles = [.1,.5,.9];

PriorQuantiles = quantile(Prior,quantiles);
PosteriorQuantiles = quantile(Posterior,quantiles);

end

