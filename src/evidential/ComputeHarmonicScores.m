% function HarmonicScores = ComputeHarmonicScores(DataStruct, ObsStruct)
% Author: Lewis Li (lewisli@stanford.edu)
% Date:    Feburary 5th 2016
%
% Function   : ComputeHarmonicScores
%
% Description: Peforms functional principal component analysis on a time series.
%			   The first step is to perform functional data analysis with a
%			   spline basis, specified by SplineOrder and SplineKnots. The
%			   coefficients are used as in the input for Principal Component
%			   analysis.
%
% Parameters :
%   DataStruct: Struct containing data/prediction variables
%   TruthStruct: Struct containing observed data
%   PlotLevel: Level 0 -> Plot nothing
%   PlotLevel: Level 1 -> Plot basis function only
%   PlotLevel: Level 2 -> Plot eigenvalues variance
%   PlotLevel: Level 3 -> Plot reconstructions
%   PlotLevel: Level 4 -> Plot everything
%
%
% Return     :
%	dataFPCA:	The harmonic scores w
%
function [dataFPCA] = ComputeHarmonicScores(DataStruct,ObsStruct,PlotLevel)

% Default behaviour is to plot nothing
if (nargin < 3)
    PlotLevel = 0;
end

FontSize = 24;

Data = DataStruct.data;
Time = DataStruct.time;
Name = DataStruct.name;

StartTime = min(Time);
EndTime = max(Time);
norder = DataStruct.spline(1);
nknots = DataStruct.spline(2);
nbasis = nknots+norder-2;

% Builds the spline basis
emptyBasis = create_bspline_basis([StartTime EndTime], nbasis, norder);

% Number of responses that Data contains
NumResponses = size(Data,3);
dataFPCA = cell(NumResponses,1);

for r = 1:NumResponses
    
    if (PlotLevel == 1 || PlotLevel == 4 && r == 1)
        h = figure;
        hBasis = plot(emptyBasis);
        set(hBasis,'LineWidth',6);
        title([ 'Order: ' num2str(norder) ' With ' ...
            num2str(nknots) ' Knots Basis B-Splines'],'FontSize',34)
        xlabel('Time (Days)');
        axis tight;
        set(gcf,'color','w');
        
        set(h, 'Units', 'normalized', 'Position', [0,0,0.75,0.75]);
        set(gcf,'color','w');
        set(gca,'FontSize',FontSize);
    end
    
    % This is the array of all realizations over all time for response
    % indexed by r
    CurrentResponse = Data(:,:,r)';
    
    % If we are computing the harmonic scores of historical data, we
    % include the actual observed data
    if (~isempty(ObsStruct))
        CurrentResponse = [CurrentResponse ObsStruct.data(:,:,r)'];
    end
    
    predFd=data2fd(CurrentResponse,Time,emptyBasis,{'Time';'Model';Name});
    dataFPCA{r} = pca_fd(predFd, nbasis);
    
    if (PlotLevel == 2 || PlotLevel == 4 && r == 1)
        %subplot(1,2,2);
        hEigen=figure;
        plot(cumsum(dataFPCA{r}.varprop)*100,'LineWidth',8); axis tight;
        xlabel('Eigen Component');
        ylabel('% Variance');
        title(['Eigenvalue Variance for ' DataStruct.ObjNames{r}]);
        axis square; axis tight;
        set(gca,'FontSize',FontSize);
        set(gcf,'color','w');
        set(hEigen, 'Units', 'normalized', 'Position', [0,0,0.75,0.75]);
    end
    
    if (PlotLevel == 3 || PlotLevel == 4 && r == 1)
        
        % Plot two random fits
        h2=figure;
        i = randi([1 size(CurrentResponse,2)]);
        subplot(1,2,1);
        plotfit_fd(CurrentResponse(:,i),Time,predFd(i));
        title(['Realization ' num2str(i) ' ' DataStruct.ObjNames{r}]);
        axis square; axis tight;
        set(gca,'FontSize',FontSize);
       
        
        i = randi([1 size(CurrentResponse,2)]);
        subplot(1,2,2);
        plotfit_fd(CurrentResponse(:,i),Time,predFd(i));
        title(['Realization ' num2str(i) ' ' DataStruct.ObjNames{r}]);
        set(gcf,'color','w');
        set(gca,'FontSize',FontSize);
        axis square; axis tight;
        set(h2, 'Units', 'normalized', 'Position', [0,0,0.75,0.75]);
    end
end

end