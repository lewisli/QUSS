function [ NextIterationModels ] = SampleNewProposalModels( ...
    SampleForParam, pdfParam,NumProposalModels,ParameterNames,PlotLevel)
%SampleNewProposalModels Summary of this function goes here
%   Detailed explanation goes here

NumParameters = size(SampleForParam,2);
NextIterationModels = zeros(NumProposalModels,NumParameters);


rng(8701);
for i = 1:NumParameters
    
    NextIterationModels(:,i) = ...
        SampleEmpiricalPDF( SampleForParam(:,i), ...
        pdfParam(:,i), NumProposalModels);
    
    if (PlotLevel > 0)
        figure(i);
        subplot(211);
        hold on;
        N = hist(NextIterationModels(:,i));
        hist(NextIterationModels(:,i));
        
        xlabel(ParameterNames{i},'FontSize',10,'FontWeight','b');
        title(['Histogram of Sampled Values for ' ParameterNames{i}]);
        subplot(212);
        hold on;
        plot(SampleForParam(:,i),pdfParam(:,i),'-b','LineWidth',2);
        
        ylabel('PDF','FontSize',10,'FontWeight','b');
        xlabel(ParameterNames{i},'FontSize',10,'FontWeight','b');
    end
end


end

