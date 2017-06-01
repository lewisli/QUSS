function [mpca_scores, mpca_obs] = MixedPCA(FunctionalStruct,truth_real,...
    eigentolerance)
%MIXEDPCA Computes Mixed PCA of a FPCA Object
%
% Inputs:
%   FunctionalStruct: Structure containing harmscr for each response
%   variable
%   truth_real[Optional]: Whether to set aside a realization as d_obs
%   eigentolerance: Number of eigenvalues to keep to maintain this variance
% Outputs:
%   mpca_scores: Scores of response variables with 99% of variance kept
%   mpca_obs: Score of observed data
%
% Author: Lewis Li (lewisli@stanford.edu)
% Date: May 24th 2016

% This is the number of variables
num_wells = length(FunctionalStruct);

% Concenated normalized scores
norm_scores = [];

for i = 1:num_wells
    % Perform regular PCA on each well
    [~,~,latent] = pca(FunctionalStruct{i}.harmscr);

    % Normalize the PCA scores by the first singular value, which is the
    norm_score = FunctionalStruct{i}.harmscr/sqrt(latent(1));

    % Concanate the norm_score
    norm_scores = [norm_scores norm_score];
end

% Perform PCA on concatenated matrix
[~,mpca_scores,~,~,explained] = pca(norm_scores);

% Compute explained variance
explained = cumsum(explained)/sum(explained);

% Check number of components to keep
eigenToKeep = 3;

if (eigentolerance<1)
    ix = max(find(explained > eigentolerance, 1, 'first'),eigenToKeep);
else
    ix = eigentolerance;
end

% Whether we set aside a truth realization
if truth_real==0
    mpca_scores = mpca_scores(:,1:ix);
    mpca_obs = 0;
else
    avail_real = setdiff(1:size(mpca_scores,1),truth_real);
    mpca_obs = mpca_scores(truth_real,1:ix);
    mpca_scores = mpca_scores(avail_real,1:ix);
end

end

