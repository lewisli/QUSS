%
% This function computes the global SA mesures for either the ASL-based
% senstitvities of standardized L1-norms

% Author: Celine Scheidt
% Date: August 2012, updated 2017


function GlobalSA = ComputeGlobalSA(SensitivityMainFactorperClass,SensitivityPerClassperBins,method)

%% Input Parameters
%   - SensitivityMainFactor: Matrix (NbParameter x NbClass) containing the standardized L1norm for each parameter (row) and
%                            each cluster (column).
%   - SensitivityPerInteractions: 4D array (NbParams x NbParams-1 x NbClusters x max(NbBins)) containing the sensitivity
%                                 values for each interaction, each class and each bin.
%   - method: 'ASL' or 'L1norm'

%% Output Parameters 
%   - GlobalSA: Vector containing the global sensitivities for each input parameters. 


NbParams = size(SensitivityPerClassperBins,1);
GlobalSA = NaN(NbParams,NbParams);

if strcmp(method,'ASL')

    SensitivityMainFactorMax = max(SensitivityMainFactorperClass,[],2);
    
    for j = 1:NbParams
        c = 0;
        for i = 1:NbParams
            if i ~=j
                GlobalSA(i,j) = max(max(SensitivityPerClassperBins(j,i-c,:,:)));
            else
                GlobalSA(i,j) = SensitivityMainFactorMax(i);
                c = 1;
            end
        end
    end
    
elseif strcmp(method,'L1norm')

    % Measure of conditional interaction sensitivity per class
    SensitivityPerClass = nanmean(SensitivityPerClassperBins,4);

    % Average measure of sensitivity for each parameter
    SensitivityPerInteractions = nanmean(SensitivityPerClass,3);
    for j = 1:NbParams
        c = 0;
        for i = 1:NbParams
            if i ~=j
                GlobalSA(i,j) = SensitivityPerInteractions(j,i-c);                
            else
                GlobalSA(i,j) = mean(SensitivityMainFactorperClass(i,:));
                c = 1;
            end
        end
    end
else
    warning('incorrect input')
end

end

