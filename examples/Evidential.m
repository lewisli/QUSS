% Evidential.m
% Tutorial file on how to use Evidential Learning for forecasting future
% prediction rates from historical production rates.
%
% Author: Lewis Li (lewisli@stanford.edu)
% Original Date: December 30th 2016
% Last Modified: Janurary 30th 2017

close all; clear all; clc;
addpath('../src/evidential/');
addpath('../src/thirdparty/fda_matlab/');
prior_path = '../data/evidential/Prior.mat';
load(prior_path);

%% Plot input responses
FontSize = 20;


% Set aside a realization to use as the "truth"
TruthRealization = 12; 
NumPriorRealizations=length(PriorData.data);
AvailableRealizations = setdiff(1:NumPriorRealizations,TruthRealization);

PlotPriorResponses(PriorData,TruthRealization,FontSize);
%PlotPriorResponses(PriorPrediction,TruthRealization,FontSize);

%% Dimension Reduction On Both Data and Prediction Variables

% Minimum numbeer of eigenvalues and % of variance to keep after dim
% reduction
MinEigenValues = 3; EigenTolerance = 0.97;

% We first perform FPCA on both d and h
PriorData.spline=[3 40]; % 3rd order spline with 40 knots
PriorDataFPCA = ComputeHarmonicScores(PriorData,4);
PriorPrediction.spline=[3 20]; % 3rd order spline with 20 knots
PriorPredictionFPCA = ComputeHarmonicScores(PriorPrediction,0);

% Perform Mixed PCA on FPCA components for d
rmpath('../src/thirdparty/fda_matlab/');
[d_f, dobs_fpca] = MixedPCA(PriorDataFPCA,TruthRealization,EigenTolerance);
addpath('../src/thirdparty/fda_matlab/');

% Get number of FPCA components to keep for h
nHarmPred = GetNumHarmonics(PriorPredictionFPCA{1}, MinEigenValues, ...
    EigenTolerance);
h_f = PriorPredictionFPCA{1}.harmscr(AvailableRealizations,1:nHarmPred);

% Plot prior models in functional space
PlotLowDimModels(d_f,h_f,dobs_fpca,'f',FontSize);

%% Apply CCA transformation and compute posterior distribution in canonical space

% Compute canonical transformation on the functional components
[A, B, ~, d_c,h_c] = canoncorr(d_f,h_f);
dobs_c=(dobs_fpca-mean(d_f))*A;

% Plot prior models in canonical space
PlotLowDimModels(d_c,h_c,dobs_c,'c',FontSize);

% Apply a normal score transform to the h_c
h_c_gauss = NormalScoreTransform(h_c,0);

% Find best linear fit between Dc and h_c_gauss
G = d_c'/h_c_gauss';

% Compute misfit covariance
DDiff= d_c'-G*h_c_gauss';
C_T = DDiff*DDiff'/length(d_c);

% Compute posterior using Linear Gaussian Regression Equations
C_H = cov(h_c_gauss);
mu_prior = mean(h_c_gauss)';
 mu_posterior = mu_prior + C_H*G'*pinv(G*C_H*G' + C_T)*(dobs_c'-G*mu_prior);
C_posterior = inv(G'*pinv(C_T)*G + inv(C_H));

%% Generate samples from the posterior in canonical space
addpath('../src/thirdparty/fda_matlab/');
NumPosteriorSamples = 100;
h_c_post = SampleCanonicalPosterior(mu_posterior,C_posterior,...
    NumPosteriorSamples,h_c);

% Undo the CCA and FPCA transformations
h_reconstructed = UndoCanonicalFunctional(h_c_post, B, h_f,...
    PriorPrediction.time, PriorPredictionFPCA);

%% Compute and plot posterior responses and quantiles
[PriorQuantiles, PosteriorQuantiles] = ComputeQuantiles(...
    PriorPrediction.data, h_reconstructed');
PlotPosteriorSamplesAndQuantiles(PriorPrediction,TruthRealization, ...
    h_reconstructed',PriorQuantiles,PosteriorQuantiles);


