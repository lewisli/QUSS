function [h_c_post] = SampleCanonicalPosterior( ...
    mu_posterior, C_posterior, NumPosteriorSamples, h_c )
%SampleCanonicalPosterior Generate samples of the canonical posterior (that
%is sampled from f(h|d_obs), and reconstructs the time series
%
%   Samples from a multivariate Gaussian determined by input mean and
%   covariance, then converts samples from canonical spcae back into time
%   domain.
%
% Inputs:
%   mu_posterior: posterior mean
%   C_posterior: posterior covariance
%   NumPosteriorSamples: number of posterior samples to generate
%   h_c: prior models in canonical space (used to undo normal score
%   transform)
%
% Outputs:
%   h_c_post: Posterior samples in canonical coefficients (NST is undone)
%
% Author: Lewis Li (lewisli@stanford.edu)
% Date: March 4th 2017

FontSize = 12;
PosteriorSamples = mvnrnd(mu_posterior',C_posterior,NumPosteriorSamples)';
PosteriorSamplesTransformed = BackTransform(PosteriorSamples,h_c);

hold on;
h1=scatter(h_c(:,1),h_c(:,2),60,[0.5 0.5 0.5]);
h2=scatter(PosteriorSamplesTransformed(1,:),...
PosteriorSamplesTransformed(2,:),60,'b','filled');
legend([h1(1),h2(1)],'Prior','Posterior');
ylabel('h^c_2','FontSize',FontSize); xlabel('h^c_1','FontSize',FontSize)
set(gca,'FontSize',FontSize); axis tight; set(gcf,'color','w');

h_c_post = PosteriorSamplesTransformed;
end

