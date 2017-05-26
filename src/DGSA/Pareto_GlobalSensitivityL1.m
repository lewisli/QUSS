%
% This function creates a Pareto plot showing the standardized measure of
% sensitivities for each parameter.
% In red, the sensitive parameters, in blue the non-sensitive parameters 
 
% Author: Celine Scheidt
% Date: August 2012. Update Nov 2016

function Pareto_GlobalSensitivityL1(SensitivityValues,ParamsNames)

%% Input Parameters
%   - SensitivityValues: vector (NbParams x NbClusters) of the parameter sensitivities
%   - ParametersNames: list containing the names of the parameters 


alpha = .95;
NbParams = length(SensitivityValues);

if NbParams > 10
    TextSize = 8;
else
     TextSize = 8;
end  

MaxSensitivity = max(SensitivityValues,[],2);

% Sort from less sensitive to most sensitive
[~, SortedSA] = sort(mean(SensitivityValues,2),'ascend');  
SensitivityValues = mean(SensitivityValues(SortedSA,:),2);
ParamsNames = ParamsNames(SortedSA);
MaxSensitivity = MaxSensitivity(SortedSA);

figure
axes('FontSize',TextSize,'fontweight','b');  hold on;
for i = 1:NbParams
        C = GetColor(alpha,MaxSensitivity(i),.80,1.2);
        barh(i,SensitivityValues(i),'FaceColor',C,'BarWidth',0.8,'LineStyle','-','LineWidth',2); 
end

set(gca,'YTick',1:NbParams)
set(gca,'YTickLabel',ParamsNames)
box on; ylim([0 NbParams+1]);

Cmap = CreateColorMap();
colormap(Cmap);
colorbar('YTick',[1.15 1.5 1.85],'YTickLabel',{'Insensitive','Important','Critical'},'fontsize',TextSize,'fontweight','b');  % to update


end



function Cmap = CreateColorMap()
   
    Cmap = zeros(100,3); 
    Cmap(1:30,1) = zeros(1,30); %001 -> 111 
    Cmap(1:30,2) = zeros(1,30);
    Cmap(1:30,3) = ones(1,30);
    Cmap(31:50,1) = linspace(0,1,20); %001 -> 111 
    Cmap(31:50,2) = linspace(0,1,20);
    Cmap(31:50,3) = ones(1,20);
    Cmap(51:70,1) = ones(1,20); %001 -> 111 
    Cmap(51:70,2) = linspace(.9,0,20);
    Cmap(51:70,3) = linspace(.9,0,20);
    Cmap(71:100,1) = ones(1,30); %111 -> 100 
    Cmap(71:100,2) = zeros(1,30);
    Cmap(71:100,3) = zeros(1,30);
end


function C = GetColor(alpha,NormalizedSA,minval,maxval)

    C = zeros(1,3);
    
    if NormalizedSA < minval
        C = [0,0,1]; % blue
    elseif NormalizedSA > maxval
        C = [1,0,0]; % red
    else
        if NormalizedSA < alpha % blue to white
            C(1) = interp1([minval,alpha],[0.1,.9],NormalizedSA);
            C(2) = interp1([minval,alpha],[0.1,.9],NormalizedSA); 
            C(3) = 1; 
        else % white to red
            C(1) = 1; 
            C(2) = interp1([alpha,maxval],[.9,0.1],NormalizedSA);
            C(3) = interp1([alpha,maxval],[.9,0.1],NormalizedSA); 
        end
    end
end