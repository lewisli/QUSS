function [ResampledStates, NewWeights, Index ] = SystematicResampling(State, ...
    InputWeights )
%SystematicResampling Perform systematic resampling on States given
%InputWeights
%
% SystematicResampling implemented as described in 
% Douc, Randal, and Olivier Cappé. 
% "Comparison of resampling schemes for particle filtering." Image and 
% Signal Processing and Analysis, 2005. ISPA 2005. 
%
% Author: Lewis Li (lewisli@stanford.edu)
% Date: Feburary 27th 2016

NumParticles = length(InputWeights);
u1 = rand/NumParticles;
u = u1;

i = 1;
j = 1;

ResampledStates = State;
NewWeights = InputWeights;
WeightsCDF = cumsum(InputWeights);
Index = zeros(length(ResampledStates),1);
while (true)

    % Construct CDF of particles
    if (WeightsCDF(i) > u)
        ResampledStates(j,:) = State(i,:);
        Index(j) = i;
        NewWeights(j) = 1/NumParticles;
        j = j+1;
        u = u1 + (j-1)/NumParticles;
    else
        i = i+1;
    end
    
    if (j > NumParticles)
        break;
    end
end

end

