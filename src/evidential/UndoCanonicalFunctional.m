function [ h_reconstructed] = UndoCanonicalFunctional(h_c_post,B, h_f,Time,...
    PriorPredictionFPCA)
%UndoCanonicalFunctional Undo the canonical and functional data analysis
%transformations
%   Posterior samples are generated in the canonical space, we need to undo
%   the canonical transform as well as the functional data transformations
%   to return the samples back into the time domain.
%
% Inputs:
%   h_c_post: posterior samples
%   B: Transformation matrix from CCA
%   Time: Vector of time steps that we need to project back onto
%   PriorPredictionFPCA: FPCA struct that contains the eigenfunctions need
%   to undo that transformation
%
% Outputs:
%   h_reconstructed: Reconstructed posterior realizations

NumPosteriorSamples = length(h_c_post);
HpostCoef = h_c_post'*pinv(B)+repmat(mean(h_f,1)',...
        1,NumPosteriorSamples)';
    
% Finally, we reconstruct the time series (mean_FDA + sum(HpostCoef*PhiH))
numPredCoeffs = size(h_f,2);
    
% Principal components for H
PhiH = eval_fd(Time,PriorPredictionFPCA{1}.harmfd);

% Re-construct time series
h_reconstructed = repmat(eval_fd(Time,PriorPredictionFPCA{1}.meanfd),...
        1,NumPosteriorSamples) + PhiH(:,1:numPredCoeffs)*...
        HpostCoef(:,1:numPredCoeffs)';
end