function [ CurrentIterationParameters,hFigures ] = SIRPredict(CurrentTime,...
    NextTime, ResampledModels,Debug )
%SIRPredict Summary of this function goes here
%   Detailed explanation goes here

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

ObservedDataPath = '../../data/simulationresults/Truth.mat';
CurrentProposalDataPath = ['../../data/simulationresults/Time' ...
    num2str(CurrentTime) '.mat'];

CurrentProposalParameterPath=['../../data/parameters/TimeStep' ...
    num2str(CurrentTime) '/Sampled.mat'];

ParameterNamePath = '../../data/parameters/ParameterNames.mat';
TrueParameterPath = '../../data/parameters/Truth/Truth.mat';
%-------------------------------------------------------------------------------

% Load parameter names
load(ParameterNamePath);

% Find which index of interpolated time is closest to next time step
InterpolatedTime = linspace(1,TotalNumberDays,TotalNumTimeSteps);
[~,HistoricalEnd] = min(abs(InterpolatedTime-NextTime));

% Due to difficulties using KDE in high dimensions, instead of mixed PCA,
% we use the sum of the fields
ForecastObjectName = {'PNEW'};
HistoricalObjectName = {'Field'};

% Load Proposal Model from current time step
load(CurrentProposalDataPath);
[HistoricalProposal,~] = ...
    GenerateDataStructsInterpolated(Data,Names,...
    ForecastColumn,HistoricalColumn,TimeColumn,HistoricalEnd,ForecastStart,...
    TotalNumTimeSteps,SplineParameters,SplineParameters,...
    ForecastObjectName,HistoricalObjectName,...
    TotalNumberDays);

% Use only resampled models
if (isempty(ResampledModels)) % empty means we use all models
    ResampledModels = 1:length(HistoricalProposal.data);
end
HistoricalProposal.data = HistoricalProposal.data(ResampledModels,:);
HistoricalProposal.RunNames = HistoricalProposal.RunNames(ResampledModels);

% Load observed data
load(ObservedDataPath);
[HistoricalTruth,~] = GenerateDataStructsInterpolated(Data,Names,...
    ForecastColumn,HistoricalColumn,TimeColumn,HistoricalEnd,ForecastStart,...
    TotalNumTimeSteps,SplineParameters,SplineParameters,ForecastObjectName,...
    HistoricalObjectName,TotalNumberDays);

% Load last iteration proposal model parameters
load(CurrentProposalParameterPath);
CurrentIterationParameters = ...
    CurrentIterationParameters(HistoricalProposal.RunNames,:);

% Load true parameter value (for visualizing only)
load(TrueParameterPath);

ProposalEmpiricalPath = ['../../data/parameters/TimeStep' ...
    num2str(NextTime)];
if exist(ProposalEmpiricalPath,'dir') ~= 7
    mkdir(ProposalEmpiricalPath);
end

% Computes updated proposal distribution using KDE on last iteration
% proposal models simulated to current iteration's time step
if (Debug == 0)
    PlotLevel = 1;
    [SampleForParam, pdfParam] = ModelProposalDistribution(...
        HistoricalProposal, HistoricalTruth, CurrentIterationParameters,...
        ParameterNames,TrueParameters,...
        EigenvalueTolerance,PlotLevel);
    save([ProposalEmpiricalPath '/Empirical.mat'],...
        'SampleForParam','pdfParam');
else
    load([ProposalEmpiricalPath '/Empirical.mat']);
    
    hFigures = PlotUpdatedParameters(SampleForParam,pdfParam,...
        CurrentIterationParameters,TrueParameters,ParameterNames,...
        CurrentTime,NextTime);
end

CurrentIterationParameters = 0;

%% We can now sample from the current iteration parameter
NumProposalModels = 500; % Generate more models than we need for now
% NumProposalModels = length(ResampledModels);
PlotLevel = 0;
if (Debug == 0)
    % Sample proposal models for current iteration
    CurrentIterationParameters = SampleNewProposalModels(SampleForParam,...
        pdfParam,NumProposalModels,ParameterNames,PlotLevel);
    % Saves sampled parameters
    save([ProposalEmpiricalPath '/Sampled.mat'],'CurrentIterationParameters');
else
    load([ProposalEmpiricalPath '/Sampled.mat']);
end


end

