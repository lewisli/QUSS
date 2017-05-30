function [ WNorm ] = ComputeModelWeights( Hc2, f2, xi2, f3, xi3 )
%ComputeModelWeights Computes a weight given proposal defined by [f2,xi2]
% [f3, xi3]
%   Detailed explanation goes here
%
% Author: Lewis Li
% Date: Feburary 27th 2016

Weights = zeros(length(Hc2),1);

for i = 1:length(Hc2)
    % value of h^(*j)
    tmp = abs(Hc2(i,1)-xi2);
    [~, idx] = min(tmp); %index of closest value
    ProposalValue = f2(idx);
    tmp = abs(Hc2(i,1)-xi3);
    [~, idx] = min(tmp); %index of closest value
    TargetValue = f3(idx);
    Weights(i) = TargetValue/ProposalValue;
end

% Normalize weights between 0 and 1
WNorm = Weights./sum(Weights);

end

