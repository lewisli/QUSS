%
% Distance-based generalized sensitivity analysis (dGSA)
% Function that plots the ASL-based values for each interaction in the
% form of a table
 
% Author: Celine Scheidt
% Date: 2015, updated 2017


function DisplayInteractionsASL(SAinteractions,ParametersNames)

%% Input Parameters
%   - SAinteractions: Matrix (NbParams x NbParams) of the global
%                            ASL-based sensitivities for each interactions
%   - ParametersNames: list containing the names of the parameters 


NbParams = size(SAinteractions,1);
alpha = 0.95; 

if ~isnan(diag(SAinteractions))
    [~,SortedParameters] = sort(diag(SAinteractions),'descend');  
    SAinteractionsSorted = SAinteractions(SortedParameters,SortedParameters);
else
    [~,SortedParameters] = sort(nansum(SAinteractions,1) + nansum(SAinteractions,2)' ,'descend');  
    SAinteractionsSorted = SAinteractions(SortedParameters,SortedParameters);
end
     
% Work on the colors
SAinteractionsSortedRB = zeros(size(SAinteractionsSorted,1),size(SAinteractionsSorted,2),3);
for i = 1:NbParams
    for j = 1:NbParams
        SAinteractionsSortedRB(i,j,:)= GetColor(alpha,SAinteractionsSorted(i,j)./100,.90,.99);
    end
end
figure; imagesc(SAinteractionsSortedRB); axis equal; axis tight
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
            if ~isnan(SAinteractionsSorted(j,i))
                if SAinteractionsSorted(j,i) > 99.9
                    text(i-.2,j,'99.9+','Color','w','FontWeight','Bold','FontSize',8);
                else
                    text(i-.2,j,num2str(SAinteractionsSorted(j,i),'%.1f'),'Color','w','FontWeight','Bold','FontSize',8);
                end
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


 function C = GetColor(alpha,Pvalue,minval,maxval)

    C = zeros(1,3);
    alphaNorm = norminv(alpha); PvalueNorm = norminv(Pvalue);
    minvalNorm = norminv(minval); maxvalNorm = norminv(maxval);
    
    if PvalueNorm < minvalNorm
        C = [0,0,1]; % blue
    elseif PvalueNorm > maxvalNorm
        C = [1,0,0]; % red
    else
        if PvalueNorm < alphaNorm % blue to white
            C(1) = interp1([minvalNorm,alphaNorm],[0.1,.9],PvalueNorm); 
            C(2) = interp1([minvalNorm,alphaNorm],[0.1,.9],PvalueNorm); 
            C(3) = 1; 
        else % white to red
            C(1) = 1;
            C(2) = interp1([alphaNorm,maxvalNorm],[.9,0.1],PvalueNorm); 
            C(3) = interp1([alphaNorm,maxvalNorm],[.9,0.1],PvalueNorm); 
        end
    end
    if isnan(PvalueNorm); C = [.95,.95,.95]; end
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