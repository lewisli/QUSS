function [dc_proposal, hc_proposal] = ComputeCanonicalOnExternalBasis(PriorData, ...
    PriorPrediction,ProposalData,ProposalPrediction,EigenTolerance,...
    Observed)
%ComputeCanonicalOnExternalBasis Computes canonical coeffcients of
%posterior samples using same eigenfunctions as prior to get them on the
%same support
%

% Author: Lewis Li
% Date: March 4th 2016

NumProposalModels = size(ProposalData.data,1);

% Perform dimension reduction on prior data/prediction to get dimensions
addpath('../src/thirdparty/fda_matlab/');
PlotLevelPCA = 0;

% If Observed is a struct
if (isstruct(Observed))
    data_FPCA = ComputeHarmonicScores(PriorData,Observed,PlotLevelPCA);
    obs_realization = size(PriorData.data,1) + 1;
else
    data_FPCA = ComputeHarmonicScores(PriorData,[],PlotLevelPCA);
    obs_realization = Observed;
end

% Perform Mixed PCA
rmpath('../src/thirdparty/fda_matlab/');
[mpca_scores, ~] = MixedPCA(data_FPCA,obs_realization,...
    EigenTolerance);
addpath('../src/thirdparty/fda_matlab/');

% Perform dimension reduction on prediction variable
pred_FPCA = ComputeHarmonicScores(PriorPrediction,[],PlotLevelPCA);

% Plot Eigenvalue Functions
MinEigenValues = 3;

% Get number of required harmonics required for forecasts
NumPredEig = GetNumHarmonics(pred_FPCA{1}, MinEigenValues,EigenTolerance);
NumDataEig = size(mpca_scores,2);

dc_proposal = zeros(NumProposalModels, NumDataEig);
hc_proposal = zeros(NumProposalModels, NumPredEig);

for i = 1:NumProposalModels
    
    % Add truth to proposal for calculating mu
    PriorData.data = cat(1,PriorData.data,...
        ProposalData.data(i,:,:));
    PriorPrediction.data = [PriorPrediction.data; ProposalPrediction.data(i,:)];
    
    [ ~, ~, d_c, h_c] = ...
        ComputePosteriorPrediction(PriorData, PriorPrediction, ...
        Observed, EigenTolerance,0,0);
    
    if (size(d_c(end,:),2) >=  NumDataEig)
        dc_proposal(i,:) = d_c(end,1:NumDataEig);
    else
        dc_proposal(i,1:MinEigenValues) = d_c(end,1:MinEigenValues);
    end
    
    if (size(h_c(end,:),2) >=  NumPredEig)
        hc_proposal(i,:) = h_c(end,1:NumPredEig);
    else
        hc_proposal(i,1:MinEigenValues) = h_c(end,1:MinEigenValues);
    end
    
   
    PriorPrediction.data(end,:)= [];
    PriorData.data(end,:,:)= [];
    
    if (mod(i,10) == 0)
        display(['Working on ' num2str(i)]);
    end
end
end

