function [mu_posterior, C_posterior, d_c, h_c, dobs_c] = ...
    ComputePosteriorPrediction(DataStruct, PredictionStruct, ...
    Observed, EigenTolerance,C_D,PlotLevel)
%ComputePosteriorPrediction Computes posterior distribution of prediction
%using Linear Gaussian Regression
%   Performs Linear Gaussian regression
%
% Inputs:
%   DataStruct: Struct containing data variable
%   PredictionStruct: Struct containing forecast responses
%   Observed: Struct or index corresponding to bserved 
%   EigenTolerance: % of variance we will use to pick how many eigenvalues
%   are kept
%   C_D: Error covariance in time domain
%   PlotLevel: Plot canonical and functional results

%
% Outputs:
%   mu_posterior: posterior mean
%   C_posterior: posterior covariance
%   Dc: historical in canonical coordinates
%   Df: historical in functional coordinates
%   Hc: forecast in canonical coordinates
%   Hf: forecast in functional coordinates
%   B: tranformation matrix from CCA
%   dobs_c: observed data in canonical coordinates
%
% Author: Lewis Li (lewisli@stanford.edu)
% Date:    Feburary 5th 2016

% By default do not plot
if (nargin < 6)
    PlotLevel = 0;
end

FontSize=12;

% Depends on if the truth realization is part of the DataStruct
PlotLevelPCA = 0;
addpath('../src/thirdparty/fda_matlab/');
if (isstruct(Observed))
    AvailableRealizations = 1:size(DataStruct.data,1);
    data_FPCA = ComputeHarmonicScores(DataStruct,Observed,PlotLevelPCA);
    obs_realization = size(DataStruct.data,1) + 1;
else
    AvailableRealizations = setdiff(1:size(DataStruct.data,1),...
        Observed);
    data_FPCA = ComputeHarmonicScores(DataStruct,[],PlotLevelPCA);
    obs_realization = Observed;
end

% Perform Mixed PCA
rmpath('../src/thirdparty/fda_matlab/');
[mpca_scores, mpca_obs] = MixedPCA(data_FPCA,obs_realization,...
        EigenTolerance);
addpath('../src/thirdparty/fda_matlab/');
    
% Perform dimension reduction on prediction variable
pred_FPCA = ComputeHarmonicScores(PredictionStruct,[],PlotLevelPCA);

% Plot Eigenvalue Functions
MinEigenValues = 3;
 
% Get number of required harmonics required for forecasts
nHarmPred = GetNumHarmonics(pred_FPCA{1}, MinEigenValues,EigenTolerance);

% Generate matrix of prediction functional components
harmscrpred=pred_FPCA{1}.harmscr;
h_f = harmscrpred(AvailableRealizations,1:nHarmPred);

% Generate matrix of data functional components
d_f = mpca_scores;
dobs_f = mpca_obs;

% Perform CCA
[A, B, ~, d_c,h_c] = canoncorr(d_f,h_f);
 
% Project dobs_f into canonical space
dobs_c=(dobs_f-mean(d_f))*A;
 
if (PlotLevel == 1)
    PlotLowDimModels(d_c,h_c,dobs_c,'c',FontSize);
    PlotLowDimModels(d_f,h_f,dobs_f,'f',FontSize);
end
 
% Perform a normal score transform
Hc_gauss = NormalScoreTransform(h_c,0);
C_H = cov(Hc_gauss);
H_CG_Mean = mean(Hc_gauss)';
 
% Find best linear bit between Dc and Hc_gauss
G = d_c'/Hc_gauss';
DDiff= d_c'-G*Hc_gauss';
C_T = DDiff*DDiff'/length(d_c);

if C_D==0
    C_Dc = zeros(size(C_T));
else
    C_Df = EstimateFunctionalErrorCovariance(DataStruct,...
        EigenTolerance,C_D);
    C_Dc = A'*C_Df*A;
end

% Perform Gaussian Regression
mu_posterior = H_CG_Mean + C_H*G'*pinv(G*C_H*G' + C_T+C_Dc)*(...
    dobs_c'-G*H_CG_Mean);
C_posterior = C_H - C_H*G'*inv(G*C_H*G' + C_T+C_Dc)*G*C_H;

end

