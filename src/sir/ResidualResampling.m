function [ResampledStates, NewWeights] = ResidualResampling(State, ...
    InputWeights)
%ResidualResampling Perform residual resampling on States given
%InputWeights
%
% Residual Resampling implemented as described in 
% Douc, Randal, and Olivier Cappé. 
% "Comparison of resampling schemes for particle filtering." Image and 
% Signal Processing and Analysis, 2005. ISPA 2005. 
%
% Author: Lewis Li (lewisli@stanford.edu)
% Date: Feburary 27th 2016


% Allocate n_i = N_particles * w_i copies of particle x_i
NumParticles = length(InputWeights);
NumCopies = floor(NumParticles*InputWeights);
ResampledStates = zeros(size(State));

% Fill in deterministic portion of particles
n = 1;
for i = 1:NumParticles
    
    for j = 1:NumCopies(i)
        ResampledStates(n,:) = State(i,:);
        n = n+1;
    end
end

% Fill in remainder
m = NumParticles - sum(NumCopies);

% Residual weights
ResidualWeights = (NumParticles*InputWeights - NumCopies)/m;

Q = cumsum(ResidualWeights);

while (n <= NumParticles)
    sampl = rand;  %(0,1]
    j = 1;
    
    while(Q(j) < sampl)
        j = j+1;
    end
    
    ResampledStates(n,:) = State(j,:);
    
    n = n +1;
end

NewWeights = ones(length(InputWeights),1)/length(InputWeights);
end

