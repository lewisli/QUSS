% function [h_reconstructed] = UndoCanonicalFunctional(h_c_post,B, h_f,Time,...
% PriorPredictionFPCA)
%
% UndoCanonicalFunctional Undo the canonical and functional data analysis
% transformations
%   Posterior samples are generated in the canonical space, we need to undo
%   the canonical transform as well as the functional data transformations
%   to return the samples back into the time domain.
%
% Inputs:
%   h_c_post: (NReals x NCanonicalDim) posterior samples
%   B: (NOriginalDim x NCanonicalDim) Transformation matrix from CCA
%   Time: (NTimeSteps) Vector of time steps that we need to project back onto
%   PriorPredictionFPCA: FPCA struct that contains the eigenfunctions need
%   to undo that transformation
%
% Outputs:
%   h_reconstructed: (NReals X NOriginalDim) Reconstructed posterior 
%   realizations

function [h_reconstructed] = UndoCanonicalFunctional(h_c_post,B, h_f,Time,...
    PriorPredictionFPCA)

% Get number of posterior samples
NumPosteriorSamples = length(h_c_post);

% Undo the canonical correlation analysis
HpostCoef = h_c_post'*pinv(B)+repmat(mean(h_f,1)',...
        1,NumPosteriorSamples)';
    
% Undo the FPCA (mean_FDA + sum(HpostCoef*PhiH))
numPredCoeffs = size(h_f,2);
PhiH = eval_fd(Time,PriorPredictionFPCA{1}.harmfd);

% Undo the basis projection of FDA
h_reconstructed = repmat(eval_fd(Time,PriorPredictionFPCA{1}.meanfd),...
        1,NumPosteriorSamples) + PhiH(:,1:numPredCoeffs)*...
        HpostCoef(:,1:numPredCoeffs)';

end