%
% This function defines the optimal (anisotropic) bandwidth for the data X using
% Silverman's rule of thumb.  Clustering is performed to use points close
% to the estimation point.
 
% Author: Celine Scheidt
% Date: June 2014


function [sigX,Clustering] = GetBandwidthData(X,Clustering)

%% Input Parameters
%   - X: matrix (N+1 x n) containing the coordinates of the N prior models projected in nD space.  
%        Each row is a point, the columns represent the dimension n of the space. The last row is
%        the location where the density must be evaluated
%   - Clustering (optional): There are two possibility for this variable:
%                 1. structure containing the clustering results (see kmedoid function). The cluster containing the observed data should be the 1st cluster.  
%                 If not provided, clustering is performed.  
%              OR:   
%                 2. nbclusters: number of clusters to construct.  If not provided, it is defined based on the silhouette index


%% Output Parameters 
%   - sigX: vector containing the bandwidth for the data
%   - Clustering: Results of the clustering, which contains:
%        - label: vector of length N + 1, which contains the cluster index that each model belongs to
%        - medoids: vector of length nbcluster containing the index of the medoids
%        - weights: vector of length nbcluster containing the number of models in each cluster.
%

% References: 
% 1. Silverman, B.W.: Density Estimation for Statistics and Data Analysis. Chapman and Hall, London (1986)
% 2. Scheidt et al.: Updating joint uncertainty in trend and depositional scenario for exploration and early appraisal,
%    submitted to Computational Geosciences, 2014 

if nargin == 2 && isstruct(Clustering) %if clustering results provided
    
    nbclusters = length(Clustering.medoids);

else    % if clustering results not provided
    
    %% Clustering with fixed medoid at the point where the density must be estimated.

    % Evaluate the distance between points
    D = pdist(X);

    % If less than 15 points, no need to apply clustering.
    if size(X,1) < 15
        nbclusters = 1;
        Clustering.label = ones(size(X,1),1);
        Clustering.weights = size(X,1);
    else
        if nargin == 1 % if the number of clusters is not provided, find the number of clusters automatically based on the silhouette index
            % Evaluate the silouhette index for 2 to 6 clusters
            s = zeros(1,5);
            for k = 2:6
                Clustering = kmedoidFixedMedoidData(D,k,20);  
                sil = silhouette(X,Clustering.label);
                s(k-1) = mean(sil);
            end
            % Optimal number of clusters corresponds to the first non-increasing
            % value of the silhoutte index
            k = 1;
            while (s(k) < s(k+1) && k < 4), k = k+1; end
            nbclusters = k+1;
        else
            nbclusters = Clustering;
        end
    %    Do the clustering with the optimal number of clusters
        Clustering = kmedoidFixedMedoidData(D,nbclusters,20);
    end

    % If the cluster containing the point where the density must be estimated
    % contains less than 3 points, group cluster with the closest one.

    if Clustering.weights(1) <= 3  
        nbclusters = nbclusters-1;
        Dmedoids = dist(X(Clustering.medoids,:)'); % distance between the medoids
        [~,idx] = min(Dmedoids(2:end,1)); idx = idx +1;
        Clustering.label(Clustering.label == idx) = 1;
        Clustering.weights(1) = Clustering.weights(1) + Clustering.weights(idx);
    end
end


% For each cluster and each dimension, the bandwidth is defined using Silverman's rule of thumb.
d = size(X,2);
sigX = zeros(nbclusters,d);
for k = 1:nbclusters
        sigX(k,:) = std(X(Clustering.label == k,:))*(4/((d+2)*Clustering.weights(k)))^(1/(d+4));
end

% If sigma is zeros, then put a small value
if any(any(sigX<=0)), sigX(sigX <=0) = .01*max(sigX(:)); end

end

    