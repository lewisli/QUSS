function [predPCA] = ComputeHarmonicScores(CFCAStruct,PlotLevel,SavePath)
%COMPUTEHARMONICSCORES Compute harmonic scores for time series using FDA/FPCA
%
% Peforms functional principal component analysis on a time series. The first 
% step is to perform functional data analysis with a spline basis, specified by 
% SplineOrder and SplineKnots. The coefficients are used as in the input for 
% Principal Component analysis.
%
% Parameters :
%   CFCAStruct: Struct of cfca object we wish to compress
%   PlotLevel: Level 0 -> Plot nothing
%   PlotLevel: Level 1 -> Plot basis function only
%   PlotLevel: Level 2 -> Plot eigenvalues variance
%   PlotLevel: Level 3 -> Plot reconstructions
%   PlotLevel: Level 4 -> Plot everything
%   SavePath: Directory to save figures
%
% Return :
%	predPCA: Struct containing the harmonic scores for object
%
% Author: Lewis Li (lewisli@stanford.edu)
% Date:    Feburary 5th 2016

% Default behaviour is to plot nothing
if (nargin < 2)
    PlotLevel = 0;
end
if (nargin < 3)
    SaveOn = false;
else
    SaveOn = true;
end

FontSize = 12;
Data = CFCAStruct.data;
Time = CFCAStruct.time;
Name = CFCAStruct.name;

StartTime = min(Time);
EndTime = max(Time);
norder = CFCAStruct.spline(1);
nknots = CFCAStruct.spline(2);
nbasis = nknots+norder-2;

% Builds the spline basis
% same spline basis
emptyBasis = create_bspline_basis([StartTime EndTime], nbasis, norder);

% Number of responses that data contains
NumResponses = size(Data,3);
predPCA = cell(NumResponses,1);

for r = 1:NumResponses
    
    if (PlotLevel == 1 && r == 1)
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
        
        if SaveOn == true
            export_fig([SavePath CFCAStruct.type '_Basis'], '-png','-m3');
        end
    end
    
    CurrentResponse = Data(:,:,r)';
    predFd=data2fd(CurrentResponse,Time,emptyBasis,{'Time';'Model';Name});
    predPCA{r} = pca_fd(predFd, nbasis);
    
    if (PlotLevel == 2 || PlotLevel == 4 && r == 1)
        hEigen=figure;
        plot(cumsum(predPCA{r}.varprop)*100,'LineWidth',8); axis tight;
        xlabel('Eigen Component');
        ylabel('% Variance');
        title(['Scree Plot ' CFCAStruct.ObjNames{r}]);
        axis square; axis tight;
        set(gca,'FontSize',FontSize);
        set(gcf,'color','w');
        set(hEigen, 'Units', 'normalized', 'Position', [0,0,0.75,0.75]);
        
        if SaveOn == true
            export_fig([SavePath CFCAStruct.type '_Eigenvalues'], '-png','-m3');
        end
    end
    
    if (PlotLevel == 3 || PlotLevel == 4 && r == 1)
        
        % Plot two random fits
        h2=figure;
        i = randi([1 size(CurrentResponse,2)]);
        subplot(1,2,1);
        plotfit_fd(CurrentResponse(:,i),Time,predFd(i));
        title(['Realization ' num2str(i) ' ' CFCAStruct.ObjNames{r}]);
        axis square; axis tight;
        set(gca,'FontSize',FontSize);
       
        
        i = randi([1 size(CurrentResponse,2)]);
        subplot(1,2,2);
        plotfit_fd(CurrentResponse(:,i),Time,predFd(i));
        title(['Realization ' num2str(i) ' ' CFCAStruct.ObjNames{r}]);
        set(gcf,'color','w');
        set(gca,'FontSize',FontSize);
        axis square; axis tight;
        set(h2, 'Units', 'normalized', 'Position', [0,0,0.75,0.75]);
        
        if SaveOn == true
            export_fig([SavePath CFCAStruct.type '_Reconstruction'], '-png','-m3');
        end
    end
end

end