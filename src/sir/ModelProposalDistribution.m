function [ SampleForParam, pdfParam] = ModelProposalDistribution(...
    Hist_Last, Hist_Obs, Param_Last,ParameterNames,TrueParameters,...
    EigenTolerance,PlotLevel)
%MODELPROPOSALDISTRIBUTION Models proposal distributions at subsequent time
%steps
%   Detailed explanation goes here
%
% Inputs:
%   Hist_Last: Struct containing responses of models from previous SIR
%   iteration, forward simulated to current time step
%   Hist_Obs: Struct containing observed data to current time step
%   Param_Last: Matrix containing parameter values of models from previous
%   time step of SIR
%   ParameterNames: Name of parameters
%   TrueParameters: True parameter value (used for plotting only)
%   EigenTolerance: Number of eigenvalues to be kept
%   PlotLevel: Controls plotting
%
% Outputs:
%   SampleForParam: 
%   pdfParam: 
%
% Author: Lewis Li (lewisli@stanford.edu)
% Date: May 30th 2016
%

% Compute harmonic scores of last iteration's data
histPCA = ComputeHarmonicScores(Hist_Last,Hist_Obs,0);

% Decide how many eigen values we are going to keep...
MinEigenValues = 2;
nHarm = max(MinEigenValues,...
    sum(cumsum(histPCA{1}.varprop)<EigenTolerance));
HarmonicScores = histPCA{1}.harmscr(:,1:nHarm);

if (PlotLevel > 0)
    scatter(histPCA{1}.harmscr(1:end-1,1),histPCA{1}.harmscr(1:end-1,2),...
        50,'filled')
    hold on;
    plot(histPCA{1}.harmscr(end,1),...
        histPCA{1}.harmscr(end,2),'^','MarkerSize',20,...
        'MarkerFaceColor','r');
end

addpath('../../thirdparty/likelihood_continuous/');
NumPDFPoints = 100;
NumParameters = size(Param_Last,2);

SampleForParam = zeros(NumPDFPoints,NumParameters);
pdfParam = zeros(NumPDFPoints,NumParameters);


for i = 1:NumParameters
    ParameterOfInterest = i;
    ParameterOfInterestName = ParameterNames{ParameterOfInterest};
    
    display(['Computing f(p_i|d_pbs) for ' ParameterOfInterestName]);
    
    [SampleForParam(:,i),pdfParam(:,i)] = UpdateProbabilityContinuous(...
        HarmonicScores,Param_Last(:,i));
    
    if (mod(i-1,4) == 0)
        figure(floor(i/4)+1);
    end
    
    [f,xi] = ksdensity(Param_Last(:,i));
    
    if (PlotLevel> 0)
        TruthParameterValue = TrueParameters.(ParameterOfInterestName);
        
        subplot(2,2,mod(i-1,4)+1);
        
        hold on;
        plot(SampleForParam(:,i),pdfParam(:,i),'-b','LineWidth',2);
        plot(xi,f,'-k','LineWidth',2);
        plot([TruthParameterValue,TruthParameterValue],...
            [0 max(pdfParam(:,i))],'r','LineWidth',3);
        ylabel('PDF','FontSize',10,'FontWeight','b');
        xlabel(ParameterNames{i},'FontSize',10,'FontWeight','b');
        
    end
end

end

