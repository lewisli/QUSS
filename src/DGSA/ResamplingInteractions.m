%
% This function performs a resampling procedure for the interactions and
% return the ASL-based sensitivity and the alpha-percentile of the L1-norm 
% for each parameter given the conditional parameter, each cluster and each bin.
 
% Author: Celine Scheidt
% Date: August 2012, updated 2017

function [ASLvalue,ResampledL1normInteractions]  = ResamplingInteractions(ParametersValues,CondIdx,Clustering,NbBins,inputStruct)

%% Input Parameters
%   - ParametersValues: matrix (NbModels x NbParams) of the parameter values
%   - CondIdx: index of the conditional parameter (e.g y in x|y)
%   - Clustering: Clustering results
%   - NbBins: number of bins for the conditional parameter
%   - inputStruct: Structure containing any of the fields 'alpha','NbDraw'and 'ObservedL1Norm'
%                         - alpha: alpha-percentile for the resampling (default is 0.95)
%                         - NbDraw: number of samples drawns for the
%                                   resampling procedure (default is 3000)
%                         - ObservedL1Norm: L1norm resulting from classification


%% Output Parameters
%   - ASLvalue: Array ((NbParams-1) x Nbcluster x NbBins) containing the ASL-based sensitivities for each
%     parameter (1st dim), each cluster (2nd dim) and for each level (3rd dim).
%   - ResampledL1normInteractions: Array ((NbParams-1) x Nbcluster x NbBins) containing the alpha-percentile for each
%     parameter (1st dim), each cluster (2nd dim) and for each level (3rd dim).


NbClusters = length(Clustering.medoids);
NbParams = size(ParametersValues,2);
ASLvalue = zeros(NbParams,NbClusters,NbBins);
ResampledL1normInteractions = zeros(NbParams,NbClusters,NbBins);

% default values
NbDraw = 3000;alpha = 0.95;ObservedL1Norm = [];
if nargin == 5
    PossibleFields = {'NbDraw'; 'alpha';'ObservedL1Norm'};
    fnames = fieldnames(inputStruct);
    for i = 1:length(fnames); if ~any(strcmp(fnames(i),PossibleFields)); warning('%s is not a possible input\n',fnames{i}); end;end;
    if isfield(inputStruct,'NbDraw'); NbDraw = inputStruct.NbDraw; end
    if isfield(inputStruct,'ObservedL1Norm'); ObservedL1Norm = inputStruct.ObservedL1Norm; end  
    if isfield(inputStruct,'alpha'); alpha = inputStruct.alpha; end    
end
    
%% Definition of levels (bin) for the conditional parameter

if length(unique(ParametersValues(:,CondIdx))) == NbBins
    levels = sort(unique(ParametersValues(:,CondIdx)));
else 
    levels = quantile(ParametersValues(:,CondIdx),(1:NbBins-1)/NbBins);
end


%% Resampling procedure 

for j = 1:NbClusters
        idx_c =  find(Clustering.T == j); % find points in the cluster
        q_prior = quantile(ParametersValues(idx_c,:),(1:1:99)./100,1);% distribution of parameter p in the entire cluster: F(p|c)
         % Bin the conditional parameter and store the indices in a vector
         for l = 1:NbBins
            if l == 1
                idx_cl = idx_c(ParametersValues(idx_c,CondIdx) <= levels(l));
            elseif l == NbBins
                idx_cl = idx_c(ParametersValues(idx_c,CondIdx) > levels(l-1));
            else
                idx_cl = idx_c(all(horzcat(ParametersValues(idx_c,CondIdx) <= levels(l),ParametersValues(idx_c,CondIdx) > levels(l-1)),2));
            end

            if  length(idx_cl) <= 3  % too few values in the bin -> nan
                ASLvalue(:,j,l) = NaN(NbParams,1);
                ResampledL1normInteractions(:,j,l) = NaN(NbParams,1);
            else
                for i =1:NbParams 
                    ResampledL1norm = zeros(NbDraw,1);
                    ParametersValuesi = ParametersValues(:,i);
                    x_redraw = zeros(length(idx_cl),NbDraw);
                    for iter = 1:NbDraw
                        x_redraw(:,iter) = ParametersValuesi(idx_c(randsample(Clustering.weights(j),length(idx_cl))));  % sample points
                    end
                    q = quantile(x_redraw,(1:1:99)./100,1);  % Resampled distribution
                    for iter = 1:NbDraw % resampling procedure
                        ResampledL1norm(iter) = norm(q_prior(:,i)-q(:,iter),1); % L1-norm
                    end
                   
                    % Evaluate the alpha-percentile of the new samples       
                    ResampledL1normInteractions(i,j,l) = quantile(ResampledL1norm,alpha);
                    
                    if ~isempty(ObservedL1Norm)
                        [ResampledL1normvalSorted,~] = sort(ResampledL1norm);
                        if i > CondIdx, idx = i - 1; else idx = i; end % interaction with itself not computed
                        if CondIdx == NbParams && i == NbParams, idx = i-1;end

                        [~,IdxMin] = min(abs(ResampledL1normvalSorted - ObservedL1Norm(idx,j,l)));
                        ASLvalue(i,j,l) = IdxMin/NbDraw;

                        if IdxMin == NbDraw % if observed L1-norm is larger than L1-norm for all samples, do extrapolation to define ASL value (by fitting a distribution)
                             ResampledL1norm(ResampledL1norm==0) = eps;
                             [D, ~] = allfitdist(ResampledL1norm); % fit a distribution on the resampled L1-norms
                             BestFitIdx = 1; GoNextBestFit = 1;
                             
                            while (~isempty(D(BestFitIdx).DistName) && GoNextBestFit == 1)
                             GoNextBestFit = 0;
                             switch D(BestFitIdx).DistName
                                 case 'lognormal'
                                     ASLvalue(i,j) = logncdf(ObservedL1Norm(idx,j,l),D(1).Params(1),D(1).Params(2));
                                 case 'gamma'
                                     ASLvalue(i,j) = gamcdf(ObservedL1Norm(idx,j,l),D(1).Params(1),D(1).Params(2));
                                 case 'weibull'
                                     ASLvalue(i,j) = wblcdf(ObservedL1Norm(idx,j,l),D(1).Params(1),D(1).Params(2));
                                 case 'beta'
                                     ASLvalue(i,j) = betacdf(ObservedL1Norm(idx,j,l),D(1).Params(1),D(1).Params(2));
                                 otherwise
                                     GoNextBestFit = 1;BestFitIdx = BestFitIdx+1;
                             end
                            end
                        end
                    end
                end
             end
        end
end

ResampledL1normInteractions(CondIdx,:,:) = [];
 if ~isempty(ObservedL1Norm)
     ASLvalue(CondIdx,:,:) = [];
     ASLvalue = ASLvalue.*100;
 else
        ASLvalue = [];
 end
    
end