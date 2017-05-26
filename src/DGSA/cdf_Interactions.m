%
% This function plots the class conditional cdfs of each parameter for bin
 
% Author: Celine Scheidt
% Date: August 2012


function cdf_Interactions(ParamsValues,ConditionalParams,Clustering,NbBins,WhichCluster,LabelName)

%% Input Parameters
%   - ParametersValues: vector (NbModels x 1) of the parameter values (e.g. x in x|y)
%   - ConditionalParams: vector (NbModels x 1) of the values of the conditional parameter values (e.g. y in x|y)
%   - Clustering: Clustering results
%   - NbBins: Number of bins for the conditional parameter
%   - WhichCluster: index of the cluster (class) for which cdfs should be plotted
%   - LabelNames: label for the interaction ('x|y')

type_line = {'--','-.',':','--','-.',':'};
C = definecolor(WhichCluster,length(Clustering.medoids),NbBins);


if length(unique(ConditionalParams)) == NbBins
    levels = sort(unique(ConditionalParams));
else 
    levels = quantile(ConditionalParams,(1:NbBins-1)/NbBins);
end

idx_c =  find(Clustering.T == WhichCluster); % points in the cluster

[fx_c,x_c] = ecdf(ParamsValues(idx_c));
figure;axes('FontSize',12,'FontWeight','b');hold on;box on;
stairs(x_c,fx_c,'LineWidth',4,'Color','k');

for l = 1:NbBins
    if l == 1
        idx_cl = idx_c(ConditionalParams(idx_c) <= levels(l));
    elseif l == NbBins
        idx_cl = idx_c(ConditionalParams(idx_c) > levels(l-1));
    else
        idx_cl = idx_c(all(horzcat(ConditionalParams(idx_c) <= levels(l),ConditionalParams(idx_c) > levels(l-1)),2));
    end
    if ~isempty(idx_cl)
        [f_xgy,xc_xgy] = ecdf(ParamsValues(idx_cl));
        stairs(xc_xgy,f_xgy,char(type_line(l)),'LineWidth',4,'Color',C(l,:));
    end
end
xlabel(LabelName,'fontsize',14,'fontweight','b');ylabel('cdf','fontsize',14,'fontweight','b'); box on;

end


%% Function defining the colors for each bin
function C = definecolor(WhichCluster,NbClusters,NbBins)

    Cs = jet(124);  
    
    % First and last cluster are blue or red
    if WhichCluster == 1
        Ds1 = floor(linspace(24,1,NbBins));
        C = Cs(Ds1,:);       
    elseif WhichCluster == max(NbClusters)
        Dsend = floor(linspace(100,124,NbBins));
        C = Cs(Dsend,:); 
    else
        % Definition of the color for the middle clusters
        Dint = [];
        Ds = floor(linspace(1,size(Cs,1),NbClusters));
        Ds = Ds(WhichCluster);
        
        if Ds+7 < 55 % i.e. in the blue color, need to reverse the color   
            Dint = floor([Dint,Ds+7:-14/(NbBins-1):Ds-7]);
        else
            Dint = floor([Dint,Ds-7:14/(NbBins-1):Ds+7]);
        end       
        C = Cs(Dint,:);       
    end
end
    


