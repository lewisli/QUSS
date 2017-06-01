function [ Samples ] = SampleEmpiricalPDF( SampleValues, pdf, NumSamples)
%SampleEmpiricalPDF Generates samples from some empirical PDF

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
    Samples(i) = SampleValues(idx);
   
end

end

