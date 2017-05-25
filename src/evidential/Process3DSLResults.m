function [ Data, PropertyNames ] = Process3DSLResults( Directory, ...
    NumRealizations, WellNames )
%Process3DSLResults Reads in results from 3DSL runs
%   Processes a set of 3DSL runs located in $Directory
%
%   Author: Lewis Li
%   Original Date: Feburary 22rd 2016

ResultNumHeaders = 4;
NumRunsRead = 0;
Data = cell(NumRealizations,1);
PropertyNames = {};
for i = 1:NumRealizations
    RunName = ['Run' num2str(i)];
    RunDirectory = [Directory RunName '/'];
    
    RealizationStruct = struct();
    
    RunComplete = 0;
    
    for w = 1:length(WellNames)
        WellName = WellNames{w};
        WellPath = [RunDirectory RunName '.' WellName '.wel'];
        
        % Checks if all well exists
        if exist(WellPath, 'file') == 2
            RunComplete = RunComplete + 1;
            
            ImportedDataStruct = importdata(WellPath,' ',ResultNumHeaders);
            PropertyNames = strsplit(char(ImportedDataStruct.textdata(...
                ResultNumHeaders)));
            PropertyNames = PropertyNames(2:end);
            
            RealizationStruct.(WellName) = ...
                ImportedDataStruct.data;
        end
    end
    
    if (RunComplete == length(WellNames))
        display([RunName ' is complete']);
        RealizationStruct.('Name') = RunName;
        
        RealizationStruct.('Field') = zeros(size(ImportedDataStruct.data));
        
        % Compute Field Production
        for w = 1:length(WellNames)
            WellName = WellNames{w};
            RealizationStruct.('Field') = RealizationStruct.('Field') + RealizationStruct.(WellName);
        end
        
        % The ts, time and bhp don't need to be summed
        RealizationStruct.('Field')(:,1:3) = RealizationStruct.('Field')(:,1:3)/length(WellNames);
        
        NumRunsRead = NumRunsRead + 1;
        Data{NumRunsRead} = RealizationStruct;
    else
        display([RunName ' is not complete... moving onto next one']);
    end
end

end