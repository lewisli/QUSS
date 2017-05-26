%
% Distance-based generalized sensitivity analysis (dGSA)
% Evaluation of sensitivity of the intercations
 
% Author: Celine Scheidt
% Date: August 2013, updated 2017

function [GlobalASLInteractions,GlobalL1Interactions,ASLvalueInteractions,L1SensitivityInteractions] = dGSA_Interactions(Clustering,ParametersValues,NbBins,ParametersNames,inputs)

%% Input Parameters
%   - Clustering: Clustering results
%   - ParametersValues: matrix (NbModels x NbParams) of the parameter values
%   - NbBins: Vector containing the number of bins per parameter
%   - ParametersNames: List containing the interaction names to be displayed on the y-axis
%   - inputs (optional). Structure containing any of the fields
%                       'alpha','NbDraw','PlotType', 'ASLMainFactors' and
%                       'L1MainFactors'
%                         - alpha: alpha-percentile for the resampling (default is 0.95)
%                         - NbDraw: number of samples drawns for the
%                                   resampling procedure (default is 3000)
%                         - PlotType: type of plot to display. Possible
%                                     choices are 'ASL','L1norm','ASLandL1norm' or
%                                     'None'.  Default is 'ASL'
%                         - ASLMainFactors: Vector (NbParams x 1) or Matrix (NbParams x NbClusters) 
%                                           containing the ASL-based global sensitivities for each input parameters. (for plotting only)
%                         - L1MainFactors: Matrix (NbParameter x NbClass) containing the normalized L1norm for each parameter (row) and
%                                           each cluster (column).(for plotting only)

%% Output Parameters
%   - GlobalASLInteractions: Matrix (NbParams x NbParams) of the global
%                            ASL-based sensitivities for each interactions
%   - GlobalL1Interactions:  Matrix (NbParams x NbParams) of the global
%                            standardized L1-norm sensitivities for each interactions
%   - ASLvalueInteractions: Array ((NbParams-1) x Nbcluster x NbBins) containing the ASL-based sensitivities for each
%                            parameter (1st dim), each cluster (2nd dim) and for each level (3rd dim).
%   - L1SensitivityInteractions: Array ((NbParams-1) x Nbcluster x NbBins) containing the standardized L1-norm sensitivities for each
%                            parameter (1st dim), each cluster (2nd dim) and for each level (3rd dim).


NbParams = size(ParametersValues,2);
NbClusters = size(Clustering.medoids,2);

% Set default values
NbDraw = 3000; alpha = .95; PlotType = 'ASL'; ASLMainFactor = NaN(NbParams,1); L1MainFactors =  NaN(NbParams,1);
if nargin == 5; 
    PossibleFields = {'NbDraw'; 'PlotType'; 'alpha';'ASLMainFactor';'L1MainFactors'};
    PossibleFieldsPlot = {'ASL';'L1norm';'ASLandL1norm'; 'None'};
    fnames = fieldnames(inputs);
    for i = 1:length(fnames); if ~any(strcmp(fnames(i),PossibleFields)); warning('%s is not a possible input\n',fnames{i}); end;end;
    if isfield(inputs,'NbDraw'); NbDraw = inputs.NbDraw; end
    if isfield(inputs,'PlotType'); if (any(strcmp(inputs.PlotType,PossibleFieldsPlot))),PlotType = inputs.PlotType; end  ; end
    if isfield(inputs,'alpha'); alpha = inputs.alpha; end    
    if isfield(inputs,'ASLMainFactor'); ASLMainFactor = inputs.ASLMainFactor; end    
    if isfield(inputs,'L1MainFactors'); L1MainFactors = inputs.L1MainFactors; end    
end

% Evaluate the conditionnal interactions for each parameter, each class and each bin
L1Interactions = NaN(NbParams,NbParams-1,NbClusters,max(NbBins));  % array containing all the Interactions 
ResampledAlphaQuantile = NaN(NbParams,NbParams-1,NbClusters,max(NbBins));  
ASLvalueInteractions = NaN(NbParams,NbParams-1,NbClusters,max(NbBins));  

for params = 1:NbParams
    % compute L1 norm
    L1InteractionsParams = L1normInteractions(ParametersValues,params,Clustering,NbBins(params)); 
    L1Interactions(params,:,:,1:NbBins(params)) = L1InteractionsParams(:,:,1:NbBins(params));
    
    % do resampling
    InputsResampling = struct('NbDraw',NbDraw,'alpha',alpha,'ObservedL1Norm',L1InteractionsParams);
    [ASLvalue,ResampledL1norm] = ResamplingInteractions(ParametersValues,params,Clustering,NbBins(params),InputsResampling);
    ResampledAlphaQuantile(params,:,:,1:NbBins(params)) = ResampledL1norm;
    ASLvalueInteractions(params,:,:,1:NbBins(params)) = ASLvalue;
end


% Compute sensitivity values
L1SensitivityInteractions = L1Interactions./ResampledAlphaQuantile;  % standardized measure of sensitivity 


%Compute and display the Global sensitivities
GlobalASLInteractions = ComputeGlobalSA(ASLMainFactor,ASLvalueInteractions,'ASL');
GlobalL1Interactions = ComputeGlobalSA(L1MainFactors,L1SensitivityInteractions,'L1norm');

if strcmp(PlotType,'ASLandL1norm'); 
    DisplayInteractionsASL(GlobalASLInteractions,ParametersNames)
    DisplayInteractionsL1(L1MainFactors,L1SensitivityInteractions,ParametersNames)
elseif strcmp(PlotType,'L1norm'); 
    DisplayInteractionsL1(L1MainFactors,L1SensitivityInteractions,ParametersNames)
elseif strcmp(PlotType,'ASL'); 
    DisplayInteractionsASL(GlobalASLInteractions,ParametersNames)
end


end



