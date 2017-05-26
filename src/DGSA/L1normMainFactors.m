%
% This function computes the L1-norm distance between the prior cdf and the
% cdf of each cluster for each parameter
 
% Author: Celine Scheidt
% Date: August 2012, updated 2014

function L1norm = L1normMainFactors(Clustering,ParametersValues)

%% Input Parameters
%   - Clustering: Clustering results
%   - ParametersValues: matrix (NbModels x NbParams) of the parameter values

%% Output Parameters 
%   - L1norm: Matrix (NbParams x NbCluster)containing the L1norm for each parameter (row) and each cluster (column).

NbParameters = size(ParametersValues,2);
NbClusters = length(Clustering.medoids);
L1norm = zeros(NbParameters,NbClusters);

q_prior = quantile(ParametersValues,(1:1:99)./100,1);  % prior distribution    
q_class = struct([]);

for j = 1:NbClusters
    q_class{j} = quantile(ParametersValues(Clustering.T == j,:),(1:1:99)./100,1);  % distribution per class
end

% Numerical calculation of the L1 norm for each parameter
for i = 1:NbParameters 
    for j = 1:NbClusters
        L1norm(i,j) = norm(q_prior(:,i) - q_class{j}(:,i),1);
    end
end

end