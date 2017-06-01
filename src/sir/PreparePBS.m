function [  ] = PreparePBS( NumRealizations,CaseName,TrialName,...
    OutputDir,ClusterName)
%PreparePBS Prepares PBS Scripts for running 3DSL Files on Cluster
%   
% Author: Lewis Li
% Date: March 4th 2016
if (nargin < 4)
    ClusterName=2010;
end

BaselineRunDirectory = [OutputDir TrialName '/pbs/'];

%creating an new folder for this iteration
%checking if there is already a folder with that name
if exist(BaselineRunDirectory,'dir') ~= 7
    mkdir(BaselineRunDirectory);
end

JobNames = cell(NumRealizations,1);

ClusterExecutable='/opt/3dsl/bin/3dsl-linux-x64';

if (ClusterName == 2010)
    ClusterDataDir = ['/data/groups/scrf/users/lewisli/3DSLRuns/' CaseName];
else
    ClusterDataDir = '/home/lewisli/Data/3DSLRuns';
end

for i = 1:NumRealizations
    JobName = [TrialName '_Run' num2str(i)];
    JobNames{i} = JobName;
    FilePath = [BaselineRunDirectory JobName '.pbs'];
    ClusterDataFilePath = [ClusterDataDir '/' TrialName ...
        '/Run' num2str(i) '/Run' num2str(i) '.dat'];
    
    fileID = fopen(FilePath,'w+');
    
    fprintf(fileID,'#!/bin/bash\n');
    fprintf(fileID,'#PBS -l nodes=1:ppn=4\n');
    fprintf(fileID,['#PBS -o ' JobName '.out\n']);
    fprintf(fileID,['#PBS -N ' JobName '\n']);
    fprintf(fileID,'#PBS -j oe\n');
    fprintf(fileID,'#PBS -l walltime=2:00:00\n');
    fprintf(fileID,'#PBS -V\n');
    fprintf(fileID,'#PBS -q default\n');
    if (ClusterName == 2010)
        fprintf(fileID,'#PBS -W x="PARTITION:sw121"\n');
    end
    fprintf(fileID,'#PBS -M lewisli@stanford.edu\n');
    fprintf(fileID,'\n');
    
    fprintf(fileID,['echo -n ' '''Job is running on node ''; cat $PBS_NODEFILE\n']);
    if (ClusterName == 2010)
        fprintf(fileID,'export RLM_LICENSE=9090@10.1.2.11\n');
    else
        fprintf(fileID,'export RLM_LICENSE=28000@rcf-nfs\n');
    end
    fprintf(fileID,'echo Looking for RLMS License at $RLM_LICENSE\n');
    fprintf(fileID,[ClusterExecutable ' ' ClusterDataFilePath '\n']);
    fprintf(fileID, '\n');
    
    fclose(fileID);
    
end

%% Make a file with all job names in it
JobListPathDir = [OutputDir TrialName '/Logs/'];
JobListPathPath = [JobListPathDir 'JobList'];
%creating an new folder for this iteration
%checking if there is already a folder with that name
if exist(JobListPathDir,'dir') ~= 7
    mkdir(JobListPathDir);
end

fileID = fopen(JobListPathPath,'w+');
for i = 1:NumRealizations
    fprintf(fileID, ['../pbs/' JobNames{i} '.pbs\n']);
end
fclose(fileID);

end

