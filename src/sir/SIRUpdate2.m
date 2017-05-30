function [ ResampledModels,N_eff,hFigures ] = SIRUpdate(CurrentTime,Debug)
%SIRIteration
%   Detailed explanation goes here
%
% Inputs:
%   CurrentTime: Time at current iteration in days


% We need to know

%---------------------------Include Directories---------------------------------
addpath('../cfca');
addpath('../sir');
addpath('../../thirdparty/fda_matlab');
%-------------------------------------------------------------------------------

%---------------------------Default Parameters----------------------------------
ForecastColumn = 4;     % Oil rate (stb/days)
HistoricalColumn = 4;   % Oil rate (stb/days)
TimeColumn = 2;         % Time (days)
ForecastStart = 133;    % Corresponds to 8000 days, when new well is drilled
TotalNumTimeSteps = 200;% Number of time steps for interpolating simulator time
TotalNumberDays = 12000;% Total number of days in simulation
SplineParameters = [4 8];

EigenvalueTolerance = 0.995;
OutlierPercentile = 100;

FontSize = 40;

ObservedDataPath = '../../data/simulationresults/Truth.mat';
PriorDataPath = '../../data/simulationresults/Prior.mat';

% Well we are trying to forecast
ForecastObjectName = {'PNEW'};

% Well(s) we are using as historical data
HistoricalObjectName = {'P1','P2','P3','P4','P5'};

%-------------------------------------------------------------------------------

% Find which index of interpolated time is closest to current time step
InterpolatedTime = linspace(1,TotalNumberDays,TotalNumTimeSteps);
[~,HistoricalEnd] = min(abs(InterpolatedTime-CurrentTime));

% Load observed data
load(ObservedDataPath);
[HistoricalTruth,ForecastTruth] = GenerateDataStructsInterpolated(Data,Names,...
    ForecastColumn,HistoricalColumn,TimeColumn,HistoricalEnd,ForecastStart,...
    TotalNumTimeSteps,SplineParameters,SplineParameters,...
    ForecastObjectName,HistoricalObjectName,...
    TotalNumberDays);

% Load Prior Data
load(PriorDataPath);
[HistoricalStruct,ForecastStruct] = GenerateDataStructsInterpolated(Data,Names,...
    ForecastColumn,HistoricalColumn,TimeColumn,HistoricalEnd,ForecastStart,...
    TotalNumTimeSteps,SplineParameters,SplineParameters,...
    ForecastObjectName,HistoricalObjectName,...
    TotalNumberDays);

% Load Proposal Model
ProposalDataPath = ['../../data/simulationresults/Time' ...
    num2str(CurrentTime) '.mat'];
if exist(ProposalDataPath, 'file') ~= 2
    display('Proposal model simulation results not found. Check path');
    return;
end
load(ProposalDataPath);
[HistoricalProposal,ForecastProposal] = ...
    GenerateDataStructsInterpolated(Data,Names,...
    ForecastColumn,HistoricalColumn,TimeColumn,HistoricalEnd,ForecastStart,...
    TotalNumTimeSteps,[4 8],[4 8],ForecastObjectName,HistoricalObjectName,...
    TotalNumberDays);

% Location to save canonical components
CoordinateSavePath = ['../../data/canonicalcoord/Time' ...
    num2str(CurrentTime) '.mat'];

if (Debug == 0)
    % Perform CFCA on prior to obtain target distribution as well as canonical
    % coordinates of prior models
    PlotLevel = 0;
    [ muTarget, CTarget, DcPrior, ~,HcPrior, ~, ~, dobs_c] = ...
        ComputeCFCAPosterior(...
        HistoricalStruct, ForecastStruct, HistoricalTruth, EigenvalueTolerance,...
        OutlierPercentile,PlotLevel,FontSize);
    
    % Projects proposal models onto same eigenbasis as prior models
    NumPredEig = size(HcPrior,2);
    NumHistEig = size(DcPrior,2);
    [ HcProposal, DcProposal ] = ComputeCanonicalOnExternalBasis(HistoricalStruct, ...
        ForecastStruct,HistoricalProposal,ForecastProposal,EigenvalueTolerance,...
        HistoricalTruth,OutlierPercentile,NumHistEig,NumPredEig);
    
    % Save the canonical components
    save(CoordinateSavePath,'muTarget','CTarget','DcPrior','HcPrior',...
        'dobs_c','HcProposal','DcProposal');
else
    % Load the canonical components
    load(CoordinateSavePath);
end

% Find out which models are left after resampled
[IndexAfterResampling, Weights,N_eff,hFigures] = ...
    CalculateModelWeightsAndResample( DcPrior, HcPrior, ...
    DcProposal, HcProposal, dobs_c, muTarget, CTarget, CurrentTime,4);

ResampledModels = unique(IndexAfterResampling);

hQuantiles = PlotUpdatedQuantiles(ForecastStruct,...
    ForecastProposal,ForecastTruth,HcProposal,Weights,CurrentTime,FontSize);

hFigures = [hFigures hQuantiles];

ResampledModels = 0;
end

