function [Historical, Forecast] = GenerateDataStructsWithInterpolation(Data, ...
    Names, ForecastColumn, HistoricalColumn, TimeColumn, ...
    HistoricalEnd,ForecastStart, NumTimeSteps, ForecastSpline, ...
    HistoricalSpline, ForecastObjName, HistoricalObjName, EndTime)

%GenerateDataStructsWithInterpolation Generates data structure for use with CFCA.
%   Converts 3DSL simulation data into seperate structures for historical
%   and forecast data to be used with PFA-CFCA. Applies an interpolation
%   step on the response curves to ensure that they are on the same time scale, 
%   as time steps for models from flow simulation may differ according to
%   convergence criteron.
%
% Inputs:
%   Data: Cell array of realizations and response curves. See --.m for
%   structure.
%   Names: Names of properties to measure 
%   ForecastColumn: Column in data array corresponding to forecast property
%   HistoricalColumn: Column in data array corresponding to historical property
%   ForecastStart: Time step at which forecast starts
%   NumTimeSteps: Total number of time steps to be used for interpolation
%   ForecastSpline: Order and number of knots for forecast spline
%   HistoricalSpline: Order and number of knots for historical spline
%   ForecastObjName: Object name for forecast (P1,P2,etc). Can only be 1
%       response
%   HistoricalObjName: Object name for historical (P1,P2,etc). Can be
%       multiple responses
%   EndTime: Final time step for simulation in days. Used for interpolating
%   the time
% Outputs:
%   Historical: Object containing historical data/properties
%   Forecast: Object containing forecast data/properties
%
%
% Author: Lewis Li (lewisli@stanford.edu)
% Date: March 4th 2016

% Number of realizations in Data
NumRealizations = sum(~cellfun('isempty',Data));

% Number of responses considered as Historical (ex: P1, P2)
NumHistoricalResponses = length(HistoricalObjName);

% Preallocate matrix for forecast and historical data
ForecastData = zeros(NumRealizations,NumTimeSteps);
HistoricalData = zeros(NumRealizations,NumTimeSteps,NumHistoricalResponses);

% Matrix to store run ids
RunNames = zeros(NumRealizations,1);

% Sampled Time Axis
Time = linspace(1,EndTime,NumTimeSteps);

% Perform CFCA to get plot of f(h|d_obs)
SuccessfulRuns = 1;

for i = 1:NumRealizations
    % Read in simulation time steps
    SimulatorTime = Data{i}.Field(:,TimeColumn);
    
    % See when simulation finished (in case some realizations did not run
    % properly on cluster)
    MaxTime = max(SimulatorTime);
    if (MaxTime >= EndTime)
        
        % Interpolate forecast response onto common time steps
        ForecastData(SuccessfulRuns,:) = interp1(SimulatorTime, Data{i}.(...
            ForecastObjName{1})(:,ForecastColumn),Time,'linear');
        
        % Interpolate historical response onto common time steps
        for j = 1:NumHistoricalResponses
            HistoricalData(SuccessfulRuns,:,j) = interp1(SimulatorTime, ...
                Data{i}.(HistoricalObjName{j})(:,HistoricalColumn),...
                Time,'linear');
        end
        
        
        % Get the IDs of the runs that did run successfully
        RunNames(SuccessfulRuns) = str2double(regexprep(Data{i}.Name,...
            '[^0-9]',''));
        SuccessfulRuns = SuccessfulRuns+1;
    else
        % Display a warning message about those that failed...
%         display(['Warning: Run ' num2str(i) ' ran only for ' ...
%             num2str(MaxTime) ' days']);
    end
    
end

RealizationIndices = 1:1:SuccessfulRuns-1;

% Generate structs for historical
Historical = struct('name',Names{HistoricalColumn});
Historical.('data') = HistoricalData(RealizationIndices,...
    1:HistoricalEnd,:);
Historical.('time') = Time(1:HistoricalEnd);
Historical.ObjNames = HistoricalObjName;
Historical.spline = HistoricalSpline;
Historical.type = 'Historical';
Historical.('RunNames') = RunNames;

% Generate structs for forecast
ForecastStart = ForecastStart + 1;
Forecast = struct('name',Names{ForecastColumn});
Forecast.('data') = ForecastData(RealizationIndices,...
    ForecastStart:end);
Forecast.('time') = Time(ForecastStart:end);
Forecast.spline = ForecastSpline;
Forecast.ObjNames = ForecastObjName;
Forecast.type = 'Forecast';
Forecast.('RunNames') = RunNames;

end

