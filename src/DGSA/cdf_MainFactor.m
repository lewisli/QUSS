%
% This function plots the cdf of each parameter for each class
 
% Author: Celine Scheidt
% Date: August 2012


function cdf_MainFactor(ParametersValues,Clustering,ParametersNames)

%% Input Parameters
%   - ParametersValues: matrix (NbModels x NbParams) of the parameter values
%   - ParametersNames: list containing the names of the parameters
%   - Clustering: Clustering results

nbparams = size(ParametersValues,2);
nbclusters = length(Clustering.medoids);
C = definecolor(nbclusters);

colmaxFig = 3;
rowFignb = nbparams/colmaxFig;
figure; 
axes('FontSize',12,'FontWeight','b');hold on;box on;

for i = 1:nbparams
    [f_prior,x_prior] = ecdf(ParametersValues(:,i));  % prior cdf

    subplot(rowFignb,colmaxFig,i);
    stairs(x_prior,f_prior,'LineWidth',4,'Color','k');
    hold on;
    
    for j = 1:nbclusters
        [f_cluster,x_cluster] = ecdf(ParametersValues(Clustering.T ==j,i));  % cdf per cluster
        stairs(x_cluster,f_cluster,'LineWidth',4,'Color',C(j,:));
    end
     xlabel(ParametersNames(i),'fontsize',14,'fontweight','b');ylabel('cdf','fontsize',14,'fontweight','b');
end


end

% Define the color for each class

function C = definecolor(nbclusters)
    Cs = jet(124);
    Ds = floor(linspace(1,size(Cs,1),nbclusters));
    C = Cs(Ds,:);
end

