%
% This function applies a k-medoid algorithm in metric space and fixes the history (observed data) as a medoid
 
% Author: Celine Scheidt
% Date: March 2014


function Clustering = kmedoidFixedMedoidData(D,nbclusters,nbrep)

%% Input Parameters
%   - D: distance matrix (N+1) x (N+1) of all the models (N prior models and the location of the observed data).  
%        The last line and column correspond to the observed data
%   - nbclusters: number of clusters to construct
%   - nbrep: (optional). Number of clustering performed. The best cluster configuration is returned.

%% Output Parameters 
%   - Clustering: Results of the clustering, which contains:
%        - label: is a vector of length N+1, which contains the index of the cluster that each model belongs to
%        - medoids: is a vector of length nbcluster containing the index of the medoids
%        - weights: is a vector of length nbcluster containing the number of models in each cluster.
%


% Reference: http://en.wikipedia.org/wiki/K-medoids
%            Kaufman, L. and Rousseeuw, P.J. (1987), Clustering by means of Medoids, 
%            in Statistical Data Analysis Based on the –Norm and Related Methods, 
%            edited by Y. Dodge, North-Holland, 405–416.


% Ensure D is a square matrix
if any(size(D) ==1), D = squareform(D,'tomatrix'); end;

maxIterations = 50;
npoints = size(D,1);
minDistBest = Inf;

if nargin < 3
    nbrep = 1;
end

for iter = 1:nbrep % nbrep clustering are performed, the best is returned
    
    % 1. Initalize: randonly select nbclusters of the npoints data points
    % as the medoids (the first medoid is the history)
    initMedoids = [npoints;randsample(npoints-1,nbclusters-1)];

    
    % 2. Associate each data point to the closest medoid
    [minDistInit, label] = min(D(initMedoids,:));

    currentMedoids = initMedoids;
    minDistCurrent = minDistInit;
    
    label_prev = NaN(1,npoints);
    nbIter = 0;
    
    while any(label ~= label_prev) && nbIter < maxIterations % while cluster configuration is changing and maxIteration not reached
              
        label_prev = label;
        
        % 3. For each medoid m
        for m = 2:nbclusters
            
            NoMedoid = setdiff(1:npoints,currentMedoids);
            NewMedoids = currentMedoids;
            
             % For each non-medoid data point o
            for o = 1:length(NoMedoid)               
                % Swap m and o and compute the cost of the configuration
                NewMedoids(m) = NoMedoid(o);
                [minDist, label] = min(D(NewMedoids,:));
                cost = sum(minDist) - sum(minDistCurrent);
                
                if cost < 0  % 4. Select the configuration with the lowest cost
                    currentMedoids(m) = NoMedoid(o);
                    [minDistCurrent, label] = min(D(currentMedoids,:));
                end
            end          
        end
        
        nbIter = nbIter+1;
    end
    
    [minDist, label] = min(D(currentMedoids,:));

    
    % Return the best clustering configuration among the nbrep tested
     if sum(minDist) < sum(minDistBest) 
        minDistBest = minDist;
        labelBest = label;
        currentMedoidsBest = currentMedoids;
    end
end

%% Once the medoids are defined, store the outputs
weights = zeros(nbclusters,1);
for i = 1:nbclusters
    weights(i) = sum(labelBest == i);
end
    
Clustering.label = labelBest;  
Clustering.medoids = currentMedoidsBest;  
Clustering.weights = weights;

end
