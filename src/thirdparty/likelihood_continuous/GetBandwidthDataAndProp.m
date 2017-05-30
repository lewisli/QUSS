%
% This function defines the optimal (anisotropic) bandwidth for the data X and the continuous property 'Prop' using
% Silverman's rule of thumb based on k-medoid clustering
 
% Author: Celine Scheidt
% Date: June 2014


function [sigDataProp,Clustering] = GetBandwidthDataAndProp(X,Prop,Clustering)

%% Input Parameters
%   - X: matrix of coordinates in nD space.  Each row is a point, the
%        columns represent the dimension n of the space. The last row is
%        the location where the density must be evaluated
%   - Prop: vector containing the values of the continuous property
%   - Clustering (optional): There are two possibility for this variable:
%                 1. structure containing the clustering results (see kmedoid function). The cluster containing the observed data should be the 1st cluster.  
%                 If not provided, clustering is performed.  
%              OR:   
%                 2. nbclusters: number of clusters to construct.  If not provided, it is defined based on the silhouette index

%% Output Parameters 
%   - sigDataProp: vector containing the bandwidth for the data and the
%                  continuous property
%   - Clustering: Results of the clustering, which contains:
%        - label: vector of length Nbmodels, which contains the cluster index each model belongs to
%        - medoids: vector of length nbcluster containing the medoids
%        - weights: vector of length nbcluster containing the number of models in each cluster.
%

% References: 
% 1. Silverman, B.W.: Density Estimation for Statistics and Data Analysis. Chapman and Hall, London (1986)
% 2. Scheidt et al.: Updating joint uncertainty in trend and depositional scenario for exploration and early appraisal,
%    submitted to Computational Geosciences, 2014 


%% Evaluate the bandwidth for the Data
if nargin == 3   
    [sigData,Clustering] = GetBandwidthData(X,Clustering); 
else
    [sigData,Clustering] = GetBandwidthData(X); 
end

%% Evaluate the bandwidth for the continuous property based on the clustering (using only models from cluster 1)
ModelsSelectedForBandwidth = (Clustering.label ==1); % take only models belonging to cluster 1 (i.e. containing the observed data)
sigProp = std(Prop(ModelsSelectedForBandwidth(1:end-1))) * (4/(3*Clustering.weights(1)-1))^(1/5); % Silverman's rule of thumb (remove medoid, as Prop is not known)

%% Combine both bandwidth
sigDataProp = [sigData(1,:),sigProp];

end

    