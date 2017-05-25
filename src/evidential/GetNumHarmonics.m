function [ nHarm ] = GetNumHarmonics(PCAObj, MinEigenValues, EigenTolerance)
%GetNumHarmonics Compute the number of functional components is required to
% capture EigenTolerance percentage of the variance.
%   The PCAObj contains the harmonic scores of a response curve projected
%   into functional space. As PCA is performed on those coefficients, most
%   of the variation should (hopefully) be contained in the first few
%   eigenvalues. 
%
% Inputs:
%   PCAObj: PCA object from FPCA
%   MinEigenValues: The min number of eigenvalues to keep
%   EigenTolerance: The percentage of variation we would like to keep
%
% Outputs:
%   nHarm: Number of functional components
%
% Author: Lewis Li (lewisli@stanford.edu)
% Date: March 24th 2016

nHarm = max(MinEigenValues,sum(cumsum(PCAObj.varprop)<EigenTolerance));

end

