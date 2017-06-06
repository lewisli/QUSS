function [ResampledModels, PosteriorQuantiles,PriorQuantiles] = ...
    UpdateSIRPosterior(PriorData, PriorPrediction, ProposalData, ...
    ProposalPrediction, ObservedData, CurrentTime, LoadCanonicalFromSave)
%UpdateSIRPosterior Computes updated prediction posterior using SIR.
%
% Given the prior, proposal and observed data, we can compute a weight for
% each proposal model, and accordingly compute the updated posterior.
%
% Inputs:
%   PriorData: Struct containing prior data
%   PriorPrediction: Struct containing prior prediction
%   ProposalData: Struct containing proposal data
%   ProposalPrediction: Struct containing proposal data
%   ObservedData: Struct containing observed data
%   CurrentTime: Current time step (days)
%   LoadCanonicalFromSave: Flag to indicate if canonical coordinates have
%   been precomputed
%
% Outputs:
%   ResampledModels: Index corresponding to proposal models kept after
%   resampling
%   PosteriorQuantiles: Quantiles corresponding to posterior predictions
%   PriorQuantiles: Quantiles corresponding to prior predictions
%
% Author: Lewis Li (lewisli@stanford.edu)
% Date: March 4th 2016

EigenTolerance = 0.95;
C_D  = 0;
PlotLevel = 0;

% First step is to compute f(h|d_{obs}) using Evidential Learning on prior
[mu_Target,C_Target,~,HcPrior, ~]=ComputePosteriorPrediction(...
    PriorData, PriorPrediction, ObservedData, EigenTolerance,C_D,PlotLevel);

% We next project the proposal models into the canonical domain

DataDir=['../data/SIR/3dslresults/Time' num2str(CurrentTime) '/'];
% Since the canonical projection takes a while to run, we will run it once
% and save it for later use.
if (LoadCanonicalFromSave == 0)
    [DcProposal,HcProposal] = ComputeCanonicalOnExternalBasis(PriorData,...
        PriorPrediction,ProposalData,ProposalPrediction,EigenTolerance,...
        ObservedData);
    save([DataDir 'ProposalCanonical.mat'], 'DcProposal','HcProposal');
else
    load([DataDir 'ProposalCanonical.mat']);
end

% Get density of f(h^\star)
[~,xPrior] = ksdensity(HcPrior(:,1));

% Apply kernel density estimation on H_c2 to get density of
% \hat{f}(h^\star|d_obs)
[fProposal,xProposal] = ksdensity(HcProposal(:,1));

% Get density of f(h^\star|d_obs)
NumPointsForDensity = 1000;
xTarget = linspace(min(xPrior),max(xPrior),NumPointsForDensity);
fTarget = normpdf(xTarget,mu_Target(1),(C_Target(1,1)));

% Compute weights of models used to calculate \hat{f}(h*|d_obs)
[WNorm] = ComputeModelWeights( HcProposal, fProposal, xProposal,...
    fTarget, xTarget );

% Perform resampling
R = [DcProposal(:,1) HcProposal(:,1)];
[~, ~, IndexAfterResampling] = ...
    SystematicResampling(R,WNorm);
ResampledModels = unique(IndexAfterResampling);

% Compute weighted posterior
[PriorQuantiles,PosteriorQuantiles] = ComputeQuantiles(PriorPrediction.data,...
    ProposalPrediction.data);

end