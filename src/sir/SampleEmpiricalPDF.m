function [ Samples ] = SampleEmpiricalPDF( Values, pdf, NumSamples)
% SAMPLEEMPIRICALPDF Sample an empirical PDF
%
% Generates samples from a PDF that has been determined empirically
%
% Inputs:
%   Values: Values corresponding to variable being sampled
%   pdf: Density of those values
%   NumSamples: Number of samples to generate
%
% Outputs:
%   Samples: Number of samples to generate
%
% Author: Lewis Li (lewisli@stanford.edu)
% Date: March 4th 2016


SumPDF = cumsum(pdf);
MaxCDF = max(SumPDF);

Samples = zeros(NumSamples,1);

% Put in a shuffle because calling MATLAB from ssh causes issues
% otherwise...
rng('shuffle');
r = (MaxCDF).*rand(NumSamples,1);

for i = 1:NumSamples
    tmp = abs(SumPDF-r(i));
    [~, idx] = min(tmp);
    Samples(i) = Values(idx);
   
end

end

