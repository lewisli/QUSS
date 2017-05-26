%
% Distance-based generalized sensitivity analysis (dGSA)
% Function that plots the standardized sensitivity values for each interaction in the
% form of a table
 
% Author: Celine Scheidt
% Date: 2015, updated 2017



function DisplayInteractionsL1(SensitivityMainFactorperClass,SensitivityPerClassperBins,ParametersNames)

%% Input Parameters
%   - SensitivityMainFactor: Matrix (NbParameter x NbClass) containing the standardized L1norm for each parameter (row) and
%                               each cluster (column).
%   - SensitivityPerInteractions: 4D array (NbParams x NbParams-1 x NbClusters x max(NbBins)) containing the sensitivity
%                             values for each interaction, each class and each bin.
%   - ParametersNames: list containing the names of the parameters 


NbParams = size(SensitivityPerClassperBins,1);

MaxperInteractions = max(max(SensitivityPerClassperBins,[],4),[],3);
MaxInteractions = zeros(NbParams,NbParams);

for j = 1:NbParams
    c = 0;
    for i = 1:NbParams
        if i ~=j
            MaxInteractions(i,j) = MaxperInteractions(j,i-c);
        else
            MaxInteractions(i,j) = max(SensitivityMainFactorperClass(i,:)); 
            c = 1;
        end
    end
end

MeanInteractions = ComputeGlobalSA(SensitivityMainFactorperClass,SensitivityPerClassperBins,'L1norm');
 
if ~isnan(SensitivityMainFactorperClass)
    [~,SortedParameters] = sort(mean(SensitivityMainFactorperClass,2),'descend');  
    MeanInteractionsSorted = MeanInteractions(SortedParameters,SortedParameters);
    MaxInteractionsSorted = MaxInteractions(SortedParameters,SortedParameters);
else
    [~,SortedParameters] = sort(nansum(MeanInteractions,1) + nansum(MeanInteractions,2)' ,'descend');  
    MeanInteractionsSorted = MeanInteractions(SortedParameters,SortedParameters);
    MaxInteractionsSorted = MaxInteractions(SortedParameters,SortedParameters);
end
    

 % Work on the colors
SAcolor = zeros(size(MeanInteractionsSorted,1),size(MeanInteractionsSorted,2),3);
for i = 1:NbParams
    for j = 1:NbParams
        SAcolor(i,j,:)= GetColor(MaxInteractionsSorted(i,j),.80,1.2);
    end
end
figure; imagesc(SAcolor); axis equal; axis tight
set(gca,'XTick',1:NbParams)
set(gca,'YTick',1:NbParams)

if NbParams < 25
    if NbParams < 10, set(gca,'XTickLabel',ParametersNames(SortedParameters),'FontWeight','b');else set(gca,'XTickLabel','','FontWeight','b'); end
    set(gca,'YTickLabel',ParametersNames(SortedParameters),'FontWeight','b')

else 
    for i = 1:NbParams
        fprintf('%i: %s\n',i,ParametersNames{SortedParameters(i)})
    end
end

Cmap = CreateColorMap();
colormap(Cmap);
hCbar = colorbar('YTick',[.12 0.5 .88],'YTickLabel',{'Insensitive','Important','Critical'},'fontsize',10,'fontweight','b');  % to update

% Add line to colorbar
vMatlab= version('-release');
if str2double(vMatlab(1:end-1)) < 2014
    line('parent',hCbar,'xdata',[-5 5],'ydata',[.25 .25],'color','k','LineWidth',2)
    line('parent',hCbar,'xdata',[-5 5],'ydata',[.75 .75],'color','k','LineWidth',2)
end
 
% Write the max p-value of each interaction
if NbParams < 15
    for i = 1:NbParams
        for j = 1:NbParams
            if ~isnan(MeanInteractionsSorted(j,i))
                text(i-.2,j,num2str(MeanInteractionsSorted(j,i),'%.2f'),'Color','w','FontWeight','Bold','FontSize',8);
            end
        end
    end
end

spacing = 1;
for row = 0.5 : spacing : NbParams+.5
  line([.5, NbParams+.5], [row, row],'color','k','linewidth',2);
end
for column = 0.5 : spacing  : NbParams+.5
  line([column, column], [.5, NbParams+.5],'color','k','linewidth',2);
end


end


%% to define the colors


 function C = GetColor(MaxInteractions,minval,maxval)

    C = zeros(1,3);
   
    if MaxInteractions < minval
        C = [0,0,1]; % blue
    elseif MaxInteractions > maxval
        C = [1,0,0]; % red
    else
        if MaxInteractions <= 1 % blue to white
            C(1) = interp1([minval,1],[0.1,.9],MaxInteractions); 
            C(2) = interp1([minval,1],[0.1,.9],MaxInteractions); 
            C(3) = 1; 
        else % white to red
            C(1) = 1;
            C(2) = interp1([1,maxval],[.9,0.1],MaxInteractions); 
            C(3) = interp1([1,maxval],[.9,0.1],MaxInteractions); 
        end
    end
    if isnan(MaxInteractions); C = [.95,.95,.95]; end
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