function [ h ] = PlotResponses(DataStruct,TruthRealization,FontSize)
%PLOTRESPONSES: Plots input response for data and prediction variables
%
% Inputs:
%   HistoricalStruct: Struct containing historical data
%   ForecastStruct: Struct containing forecast data
%   TruthRealization: Index of realization taken to be the truth
%   FontSize: Font size in plots
%
% Author: Lewis Li (lewisli@stanford.edu)
% Original Date: March 4th 2016
% Last Updated: September 26th 2016

if (nargin < 3)
    FontSize=12;
end

NumRealizations = size(DataStruct.data,1);
NumResponses = size(DataStruct.data,3);

MaxCols = min(3,NumResponses);
NumRows = 2;

MaxFiguresPerPage = 6;

for i = 1:NumResponses
    if mod(i,MaxFiguresPerPage) == 1
        FigID =fix(i/MaxFiguresPerPage) + 1;
        h=figure('Units', 'normalized', 'Position', [0,0,1,1]);
        NumCols = min(NumResponses - MaxFiguresPerPage*(FigID-1),MaxCols);
    end
    
    subplot(NumRows,MaxCols,mod(i-1,MaxFiguresPerPage)+1);
    
    
    hold on; title([ DataStruct.type ': ' DataStruct.ObjNames{i}],...
        'FontSize',FontSize); axis square;
    
    for j = 1:NumRealizations
        h1 = plot(DataStruct.time, DataStruct.data(j,:,i)',...
            'color',[0.5 0.5 0.5]);
    end
    xlabel('Time (days)','FontSize',FontSize);
    ylabel(DataStruct.name,'FontSize',FontSize);
    
    
    if (isstruct(TruthRealization))
        h2 = plot(TruthRealization.time, ...
            TruthRealization.data(:,:,i)','r','LineWidth',3);
        hlegend = legend([h1,h2],'Prior','Observed');
        set(hlegend,'Location','southwest');
    elseif TruthRealization>0
        h2 = plot(DataStruct.time, ...
            DataStruct.data(TruthRealization,:,i)','r','LineWidth',3);
        hlegend = legend([h1,h2],'Prior','Observed');
        set(hlegend,'Location','southwest');
    end
    axis tight;
    set(gcf,'color','w');
    set(gca,'FontSize',FontSize);
    
end

end

