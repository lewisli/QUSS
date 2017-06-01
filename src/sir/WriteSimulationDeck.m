function WriteSimulationDeck(BaseCaseDatPath, ModelParameterValues, ...
    OutputFolder, TrialName )
%WriteSimulationDeck Generates an 3DSL simulation deck for a given base
%case and set of parameters values sampled using regular Monte Carlo
%
%   Generates a simulation deck for the Libyan Case. Assumes that
%   ModelParametersValues contains 12 uncertain parameters (Swir, Swor,
%   krw_end, kro_end,no,nw,FaultMulti1,FaultMulti2,FaultMulti3,FaultMulti4,
%   Viscosity, OWC).
% Inputs:
%   BaseCaseDatPath: Path to base line simulation deck
%   ModelParameterValues: Struct containing the sampled model parameters
%   OutputFolder: Path to write the output simulation decks
%   TrialName: Name of current set of runs
%
% Notes: No error checking is currently implemented, assumes each parameter
% has the same number of realizations and that all 12 parameters are
% present.
%
% Author: Lewis Li (lewisli@stanford.edu)
% Date: May 25th 2016

% Figure out how many realizations we have
ParameterNames = fieldnames(ModelParameterValues);
NbSimu = size(ModelParameterValues.(ParameterNames{1}),1);

% Compute relative permabilities of oil and water using Corey expressions
[ SW_corey_m, krw_model, kro_model ] = GenerateCoreyCurves(...
    ModelParameterValues,NbSimu,0);

% Get number of discretizations used for constructing rel perms
RelPermEntries = size(SW_corey_m,2);

% Allocate a cell array that we will use to store the baseline
s=cell(GetNumberOfLines(BaseCaseDatPath),1);

% Read from base case
fid = fopen(BaseCaseDatPath);
lineCt = 1;
tline = fgetl(fid);

while ischar(tline)
    s{lineCt} = (tline);
    lineCt = lineCt + 1;
    tline = fgetl(fid);
end

% Directory to store output simulation decks
OutputDirectory = [OutputFolder '/' TrialName '/'];

% Generate a seperate deck for each 
for k=1:NbSimu
    
    FolderNameIteration = [OutputDirectory 'Run', num2str(k)];
    
    %creating an new folder for this iteration
    %checking if there is already a folder with that name
    if exist(FolderNameIteration,'dir') ~= 7
        mkdir(FolderNameIteration);
    end
    
    % Create new file
    file_name = [FolderNameIteration '/Run', num2str(k), '.dat'];
    fileID = fopen(file_name,'w+');
    
    % Loading everything before the Faultmultiplier
    MULTFLTIndex = SearchCellArray('MULTFLT',s);
    for j=1:MULTFLTIndex
        fprintf(fileID,'%c',s{j});
        fprintf(fileID,'\n');
    end
    
    
    % Fault Multiplier
    F1Mult=['fault_1' blanks(1) ...
        num2str(ModelParameterValues.('FaultMulti1')(k)) blanks(1) '/'];
    F2Mult=['fault_2' blanks(1) ...
        num2str(ModelParameterValues.('FaultMulti2')(k)) blanks(1) '/'];
    F3Mult=['fault_3' blanks(1) ...
        num2str(ModelParameterValues.('FaultMulti3')(k)) blanks(1) '/'];
    F4Mult=['fault_4' blanks(1) ...
        num2str(ModelParameterValues.('FaultMulti4')(k)) blanks(1) '/'];
    
    fprintf(fileID,'%s',F1Mult);
    fprintf(fileID,'\n');
    fprintf(fileID,'%s',F2Mult);
    fprintf(fileID,'\n');
    fprintf(fileID,'%s',F3Mult);
    fprintf(fileID,'\n');
    fprintf(fileID,'%s',F4Mult);
    fprintf(fileID,'\n');
    fprintf(fileID,'/ \n');
    
    PVMULTIndex = SearchCellArray('PVMULT',s);
    CVISCOSITIESIndex = SearchCellArray('CVISCOSITIES',s);
    
    % Writing everything before Viscosity
    for j=PVMULTIndex:CVISCOSITIESIndex
        fprintf(fileID,'%c',s{j});
        fprintf(fileID,'\n');
    end
    
    % Printing the PVTs out - in this case only viscosity is changed
    Visc=[num2str(ModelParameterValues.('Viscosity')(k)) blanks(1) '0.1 0.4 /'];
    fprintf(fileID,'%s',Visc);
    
    KRWOIndex = SearchCellArray('KRWO',s);
    
    % Writing everything before Viscosity
    for j=CVISCOSITIESIndex+1:KRWOIndex
        fprintf(fileID,'%c',s{j});
        fprintf(fileID,'\n');
    end
    
    % Writing out the Relperms
    
    formatSpecRelPerm = '%4.4f %4.8f %4.8f %s\n';
    fprintf(fileID,'%s\n', '--    Sw        krw       kro      Pc');
    for j=1:RelPermEntries
        fprintf(fileID,formatSpecRelPerm,SW_corey_m(k,j),...
            krw_model(k,j),kro_model(k,j));
        fprintf(fileID,'\n');
    end
    fprintf(fileID, '/\n');
    fprintf(fileID,'%s\n', 'END RELPERMS');
    %check if the slash works in the simulation
    
    INITIALCOND = SearchCellArray('INITIALCOND',s);
    OWC = SearchCellArray('OWC',s);
    
    % Writing everything before OWC
    for j=INITIALCOND:OWC
        fprintf(fileID,'%c',s{j});
        fprintf(fileID,'\n');
    end
    
    % Writing out the OWC
    OWCValue=['-', num2str(ModelParameterValues.('OWC')(k))];
    fprintf(fileID,'%s',OWCValue);
    fprintf(fileID,'\n');
    fprintf(fileID, '/\n');
    
    EndInitialCondition=SearchCellArray('END INITIALCOND',s);
    
    % Writing out the rest
    for j=EndInitialCondition:GetNumberOfLines(BaseCaseDatPath)
        fprintf(fileID,'%c',s{j});
        fprintf(fileID,'\n');
    end
    
    fclose(fileID);
end

end

