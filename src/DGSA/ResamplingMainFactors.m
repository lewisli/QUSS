%
% This function performs a resampling procedure and returns the
% ASL-based sensitivity value and the alpha-percentile value for each parameter and each cluster.
 
% Author: Celine Scheidt
% Date: August 2012, updated 2017

function [ASLvalue,FactorResampledQuantile] = ResamplingMainFactors(ParametersValues,Clustering,inputStruct)

%% Input Parameters
%   - ParametersValues: matrix (NbModels x NbParams) of the parameter values
%   - Clustering: Clustering results
%   - inputStruct: Structure containing any of the fields 'alpha','NbDraw'and 'ObservedL1Norm'
%                         - alpha: alpha value to return alpha-percentile from the resampling (default is 0.95)
%                         - NbDraw: number of samples drawn during the
%                                   resampling procedure (default is 3000)
%                         - ObservedL1Norm: Matrix (NbParams x NbClusters) of L1norm per class resulting from classification

%% Output Parameters 
%   - ASLvalue: Matrix (NbParams x NbClusters) containing the ASL-based sensitivities for each
%               parameter (row) and each cluster (column). If observed L1-norm not
%               provided, return an empty matrix
%   - FactorResampledQuantile: Matrix (NbParams x NbClusters) containing the alpha-percentile for each
%                             parameter (row) and each cluster (column).

    nbParameters = size(ParametersValues,2);
    nbclusters = length(Clustering.medoids);
    FactorResampledQuantile = zeros(nbParameters,nbclusters);
    ASLvalue = zeros(nbParameters,nbclusters);
    
    % default values
    NbDraw = 3000;alpha = 0.95;ObservedL1Norm = [];
    if nargin == 3
        PossibleFields = {'NbDraw'; 'alpha';'ObservedL1Norm'};
        fnames = fieldnames(inputStruct);
        for i = 1:length(fnames); if ~any(strcmp(fnames(i),PossibleFields)); warning('%s is not a possible input\n',fnames{i}); end;end;
        if isfield(inputStruct,'NbDraw'); NbDraw = inputStruct.NbDraw; end
        if isfield(inputStruct,'ObservedL1Norm'); ObservedL1Norm = inputStruct.ObservedL1Norm; end  
        if isfield(inputStruct,'alpha'); alpha = inputStruct.alpha; end    
    end
    
    
     q_prior = quantile(ParametersValues,(1:1:99)./100,1);  % prior distribution    

    % Do the resampling procedure for each parameter and each cluster
    for i = 1:nbParameters
        ParametersValuesi = ParametersValues(:,i);
        for  j = 1:nbclusters
            ResampledL1norm = zeros(NbDraw,1);
            x_redraw = zeros(Clustering.weights(j),NbDraw);
           
            % Apply resampling
            for iter = 1:NbDraw 
                x_redraw(:,iter) = ParametersValuesi(randsample(size(ParametersValues,1),Clustering.weights(j)));  
            end
            q = quantile(x_redraw,(1:1:99)./100,1);  
            for iter = 1:NbDraw
                ResampledL1norm(iter) = norm(q_prior(:,i)-q(:,iter),1); % L1-norm between prior and resampled CDFs
            end
            
            % Evaluate the alpha-percentile of the new samples
            FactorResampledQuantile(i,j) = quantile(ResampledL1norm,alpha);
            
            %% evaluate the ASL if observed L1-norm provided
            if ~isempty(ObservedL1Norm)
                [ResampledL1normSorted,~] = sort(ResampledL1norm);
                [~,IdxMin] = min(abs(ResampledL1normSorted - ObservedL1Norm(i,j)));
                ASLvalue(i,j) = IdxMin/NbDraw;

                if IdxMin == NbDraw % if observed L1-norm larger than all samples, do extrapolation to define ASL value
                     ResampledL1norm(ResampledL1norm==0) = eps;
                     [D, ~] = allfitdist(ResampledL1norm); % fit a distribution on the resampled L1-norms
                     BestFitIdx = 1; GoNextBestFit = 1;
                    while (~isempty(D(BestFitIdx).DistName) && GoNextBestFit == 1)
                     GoNextBestFit = 0;
                     switch D(BestFitIdx).DistName
                         case 'lognormal'
                             ASLvalue(i,j) = logncdf(ObservedL1Norm(i,j),D(1).Params(1),D(1).Params(2));
                         case 'gamma'
                             ASLvalue(i,j) = gamcdf(ObservedL1Norm(i,j),D(1).Params(1),D(1).Params(2));
                         case 'weibull'
                             ASLvalue(i,j) = wblcdf(ObservedL1Norm(i,j),D(1).Params(1),D(1).Params(2));
                         case 'beta'
                             ASLvalue(i,j) = betacdf(ObservedL1Norm(i,j),D(1).Params(1),D(1).Params(2));
                         otherwise
                             GoNextBestFit = 1;BestFitIdx = BestFitIdx+1;
                     end
                    end

                end
            end
        end
    end
    if ~isempty(ObservedL1Norm)
        ASLvalue = ASLvalue.*100;
    else
        ASLvalue = [];
    end
end

