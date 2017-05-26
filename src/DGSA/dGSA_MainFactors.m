%
% Distance-based generalized sensitivity analysis (DGSA)
% Evaluation of sensitivity of the main factors only
% Pareto plots are used to display the resuls
 
% Author: Celine Scheidt
% Date: August 2012, updated March 2017

function [GlobalASLvalue,StandardizedSensitivity,ASLvalueperClass,SensitivityMainFactors] = dGSA_MainFactors(Clustering,ParametersValues,ParametersNames,inputs)

%% Input Parameters
%   - Clustering: Clustering results, obtained from the kmedoids function
%   - ParametersValues: matrix (NbModels x NbParams) of the parameter values
%     (numerical values should be provided, even for discrete parameters)
%   - ParametersNames: list containing the names of the parameters 
%   - inputs (optional). Structure containing any of the fields 'alpha','NbDraw' and 'PlotType'
%                         - alpha: alpha value to return alpha-percentile from the resampling (default is 0.95)
%                         - NbDraw: number of samples drawn during the
%                                   resampling procedure (default is 3000)
%                         - PlotType: type of plot to display. Possible
%                         choices are 'ASL','L1norm','ASLandL1norm' or 'None'



%% Output Parameters 
%   - GlobalASLvalue: Vector containing the global ASL-based sensitivities for each input parameters. 
%   - StandardizedSensitivity: Vector containing the global standardized measure of sensitivity for each parameter
%   - ASLvalueperClass: Matrix (NbParameter x NbClass) containing the ASL-based sensitivities for each parameter (row) and for each
%                       class (column)
%   - SensitivityMainFactors: Matrix (NbParameter x NbClass) containing the normalized L1norm for each parameter (row) and
%                             each cluster (column).


% Set default values
NbDraw = 3000; alpha = .95; PlotType = 'ASL'; 
if nargin == 4; 
    PossibleFields = {'NbDraw'; 'PlotType'; 'alpha'};
    PossibleFieldsPlot = {'ASL';'L1norm';'ASLandL1norm'; 'None'};
    fnames = fieldnames(inputs);
    for i = 1:length(fnames); if ~any(strcmp(fnames(i),PossibleFields)); warning('%s is not a possible input\n',fnames{i}); end;end;
    if isfield(inputs,'NbDraw'); NbDraw = inputs.NbDraw; end
    if isfield(inputs,'PlotType'); if (any(strcmp(inputs.PlotType,PossibleFieldsPlot))),PlotType = inputs.PlotType; end  ; end
    if isfield(inputs,'alpha'); alpha = inputs.alpha; end    
end

if any(Clustering.weights < 5); warning('Too small cluster size - the DGSA results may be inacurate\n'); IDX = (Clustering.weights < 5); else IDX = []; end
if length(Clustering.T) ~= size(ParametersValues,1); error('The number of input parameter combination is different than the number of models in the clustering\n'); end

% Evaluate L1-norm in each class
L1MainFactors = L1normMainFactors(Clustering,ParametersValues);

% Resampling
InputsResampling = struct('NbDraw',NbDraw,'alpha',alpha,'ObservedL1Norm',L1MainFactors);
[ASLvalueperClass,ResampledQuantile] = ResamplingMainFactors(ParametersValues,Clustering,InputsResampling);

% Remove cluster with too few data (if any)
ASLvalueperClass(:,IDX) = [];
ResampledQuantile(:,IDX) = [];


% Compute sensitivity values
SensitivityMainFactors = L1MainFactors./ResampledQuantile;  % standardized measure of sensitivity 
StandardizedSensitivity = mean(SensitivityMainFactors,2); % global standardized sensitivty measures
GlobalASLvalue = max(ASLvalueperClass,[],2);

% Display the Global sensitivities
if strcmp(PlotType,'ASLandL1norm'); 
    Pareto_GlobalSensitivityASL(GlobalASLvalue,ParametersNames,alpha);
    Pareto_GlobalSensitivityL1(SensitivityMainFactors,ParametersNames);
elseif strcmp(PlotType,'L1norm'); 
    Pareto_GlobalSensitivityL1(SensitivityMainFactors,ParametersNames);
elseif strcmp(PlotType,'ASL'); 
    Pareto_GlobalSensitivityASL(GlobalASLvalue,ParametersNames,alpha);
end

end