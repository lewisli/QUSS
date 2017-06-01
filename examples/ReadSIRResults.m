CurrentTime = 1000;
% Load Proposal Model
ProposalDataPath = ['../data/SIR/simulationresults/Time' ...
    num2str(CurrentTime) '.mat'];
if exist(ProposalDataPath, 'file') ~= 2
    display('Proposal model simulation results not found. Check path');
    return;
end
load(ProposalDataPath);

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

FontSize = 24;

% Well we are trying to forecast
ForecastObjectName = {'PNEW'};

% Well(s) we are using as historical data
HistoricalObjectName = {'P1','P2','P3','P4','P5'};

% Find which index of interpolated time is closest to current time step
InterpolatedTime = linspace(1,TotalNumberDays,TotalNumTimeSteps);
[~,HistoricalEnd] = min(abs(InterpolatedTime-CurrentTime));


[ProposalData,ProposalPrediction] = ...
    GenerateDataStructsWithInterpolation(Data,Names,...
    ForecastColumn,HistoricalColumn,TimeColumn,HistoricalEnd,ForecastStart,...
    TotalNumTimeSteps,[4 8],[4 8],ForecastObjectName,HistoricalObjectName,...
    TotalNumberDays);

%% TODO:
% Make a prior, proposal and observed for each time step