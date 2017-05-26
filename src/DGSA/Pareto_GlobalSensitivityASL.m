%
% This function creates a Pareto plot showing the ASL-based sensitivities for each parameter.
% In red, the sensitive parameters, in blue the non-sensitive parameters 
 
% Author: Celine Scheidt
% Date: May 2016

function Pareto_GlobalSensitivityASL(SensitivityValues,ParamsNames,alpha)

%% Input Parameters
%   - SensitivityValues: vector (NbParams x 1) of the parameter sensitivities
%   - ParamsNames: Parameter names
%   - alpha: alpha value to return alpha-percentile from the resampling (default is 0.95)

NbParams = length(SensitivityValues);
if nargin <3; alpha = 0.95; end

if NbParams > 10
    TextSize = 10;
else
     TextSize = 8;
end  

% Sort from less sensitive to most sensitive parameter
[~, SortedSA] = sort(SensitivityValues(:),'ascend');  
SensitivityValues = SensitivityValues(SortedSA);
ParamsNames = ParamsNames(SortedSA);

SensitivityMainFactorsZscore = norminv(SensitivityValues./100,2,1); % to stretch the x axis for better vizualization
if any(~isfinite(SensitivityMainFactorsZscore)); SensitivityMainFactorsZscore(~isfinite(SensitivityMainFactorsZscore)) = max(SensitivityMainFactorsZscore(isfinite(SensitivityMainFactorsZscore)))+eps; end

figure
axes('FontSize',TextSize,'fontweight','b');  hold on;
for i = 1:NbParams
        C = GetColor(alpha,SensitivityValues(i)./100,.90,.99);
        barh(i,SensitivityMainFactorsZscore(i),'FaceColor',C,'BarWidth',0.8,'LineStyle','-','LineWidth',2); % blue
end

for i = 1:NbParams
    if (SensitivityValues(i) > 90) && (SensitivityValues(i) < 95)  
        text(norminv(.95,2,1)+.1,i,num2str(SensitivityValues(i),'%.1f'),'fontsize',8,'fontweight','b')
    else
        if SensitivityValues(i) > 99.9
            text(SensitivityMainFactorsZscore(i)+.1,i,'99.9+','fontsize',8,'fontweight','b')
        else
            text(SensitivityMainFactorsZscore(i)+.1,i,num2str(SensitivityValues(i),'%.1f'),'fontsize',8,'fontweight','b')
        end
    end
end
plot(norminv([alpha alpha],2,1),[0.5 NbParams+0.5],'--k','LineWidth',2) 
set(gca,'YTick',1:NbParams)
set(gca,'YTickLabel',ParamsNames)

box on; ylim([0 NbParams+1]);
xlim([0 max(SensitivityMainFactorsZscore(:))+1]);

ThickToDisplay = [5, 20, 50, 80, 95,99.3,99.95,99.999];
   
set(gca,'XTick',norminv(ThickToDisplay./100,2,1))
set(gca,'XTickLabel',ThickToDisplay)
Cmap = CreateColorMap();
colormap(Cmap);
colorbar('YTick',[1.15 1.5 1.85],'YTickLabel',{'Insensitive','Important','Critical'},'fontsize',TextSize,'fontweight','b'); 

end

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