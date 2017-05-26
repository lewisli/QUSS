%
% This function computes the L1-norm distance between the class-conditional 
% distribution of a single parameter and the class-conditional distribution
% of that parameter additionally conditioned to a second parameter. 

% By default, the distance is computed for all parameters, except the conditional parameter
 
% Author: Celine Scheidt
% Date: August 2012

function InteractionSensitivity = L1normInteractions(ParametersValues,CondIdx,Clustering,NbBins,WhichParam)

%% Input Parameters
%   - ParametersValues: matrix (NbModels x NbParams) of the parameter values
%   - CondIdx: index of the conditional parameter (e.g. y in x|y)
%   - Clustering: Clustering results
%   - NbBins: number of bins for the conditional parameter
%   - WhichParam (optional).  To compute only the distance for one parameter
%                             (ParametersValues(:,WhichParam)|ParametersValues(:,CondIdx))

%% Output Parameters 
%   - InteractionSensitivity: Array ((NbParams-1) x Nbcluster x NbBins) containing the L1-norm for each
%   parameter (i.e x|y, z|y, t|y, ..), 1st dim), each cluster (2nd dim) and for each level (3rd dim).


NbClusters = length(Clustering.medoids);
NbParams = size(ParametersValues,2);

%% definition of the levels (bin) for the conditional parameter

if length(unique(ParametersValues(:,CondIdx))) == NbBins
    levels = sort(unique(ParametersValues(:,CondIdx)));
else 
    levels = quantile(ParametersValues(:,CondIdx),(1:NbBins-1)/NbBins);
end


%%  Compute the distance
if nargin < 5
    InteractionSensitivity = zeros(NbParams,NbClusters,NbBins);
else
    InteractionSensitivity = zeros(NbClusters,NbBins);
end

for j = 1:NbClusters
    idx_c =  find(Clustering.T == j); % find points in the cluster
    q_prior = quantile(ParametersValues(idx_c,:),(1:1:99)./100,1);  % distribution of parameters in the entire cluster: F(p|c)
    
    % Bin the conditional parameter and store the indices in a vector
    for l = 1:NbBins
        if l == 1
            idx_cl = idx_c(ParametersValues(idx_c,CondIdx) <= levels(l));
        elseif l == NbBins
            idx_cl = idx_c(ParametersValues(idx_c,CondIdx) > levels(l-1));
        else
            idx_cl = idx_c(all(horzcat(ParametersValues(idx_c,CondIdx) <= levels(l),ParametersValues(idx_c,CondIdx) > levels(l-1)),2));
        end

        % Compute the L1norm
        if nargin < 5 % compute distance for all parameters  
            q_inter = quantile(ParametersValues(idx_cl,:),(1:1:99)./100,1); % distribution of parameters conditioned to pj in bin l: F(p|i(pi,tl),c)
            for i =1:NbParams 
                if i == CondIdx
                    continue
                else
                    InteractionSensitivity(i,j,l) = norm(q_prior(:,i)-q_inter(:,i),1);    % L1-norm
                end
            end
        else % compute only the given interaction
            q_prior = quantile(ParametersValues(idx_c,WhichParam),(1:1:99)./100);  % distribution of parameter p in the entire cluster: F(p|c)
            q_inter = quantile(ParametersValues(idx_cl,WhichParam),(1:1:99)./100); % distribution of parameter p conditioned to pj in bin l:F(p|i(pi,tl),c)
            InteractionSensitivity(j,l) = norm(q_prior-q_inter,1);                 % L1-norm
        end       
    end
end

if nargin < 5
    InteractionSensitivity(CondIdx,:,:) = [];
end
    
end

