close all; clear all;
addpath('../src/evidential/');
addpath('../src/thirdparty/fda_matlab/');

load('../data/SIR/3dslresults/Time1000/Prior.mat');
load('../data/SIR/3dslresults/Time1000/Observed.mat')
load('../data/SIR/parameters/TimeStep0000/Sampled.mat');
load('../data/SIR/parameters/ParameterNames.mat');

%%
ObservedData = TrueData;
[ ProposalParameterVal,ProposalPDFs] = EstimateProposalDistribution(...
    PriorData,TrueData,CurrentIterationParameters);


addpath('../src/sir');
NumProposalModels = 500;
PlotLevel = 0;
ProposalModelParameters = SampleNewProposalModels(ProposalParameterVal, ...
    ProposalPDFs,NumProposalModels,ParameterNames,PlotLevel);

%%
clear all;
LoadCanonicalFromSave=1;

% Load last iteration proposal model parameters
load('../data/SIR/parameters/TimeStep1000/Sampled.mat');
LastIterationProposalParameters = CurrentIterationParameters;

load('../data/SIR/3dslresults/Time2000/Observed.mat')
load('../data/SIR/3dslresults/Time2000/Prior.mat')


% ObservedData = TrueData;
% [ ProposalParameterVal,ProposalPDFs] = EstimateProposalDistribution(...
%     PriorData,ObservedData,LastIterationProposalParameters);

load('../data/SIR/3dslresults/Time2000/Proposal.mat')

for CurrentTime = 3000:1000:8000
    
    data_dir = ['../data/SIR/3dslresults/Time' num2str(CurrentTime) '/'];
    load([data_dir '/Observed.mat']);
    load([data_dir '/Prior.mat']);
    load([data_dir '/Proposal.mat']);
    
    [ResampledModels, PosteriorQuantiles,PriorQuantiles] = UpdateSIRPosterior(...
        PriorData, PriorPrediction, ProposalData, ProposalPrediction, ...
        TrueData,CurrentTime,LoadCanonicalFromSave);
    figure;
    hold on;
    h1 = plot(PriorPrediction.time, PosteriorQuantiles,'b--',...
        'linewidth',3);
    h2 = plot(PriorPrediction.time, PriorQuantiles,'color',[0.5, 0.5, 0.5],...
        'linewidth',3);
    h3 = plot(PriorPrediction.time, TruePrediction.data,'color','r',...
        'linewidth',3);
    legend([h1(1),h2(1),h3(1)],'SIR','Prior','Reference');
    title(['SIR After ' num2str(CurrentTime) ' Days']);
    axis tight; xlabel('Time (Days)'); ylabel( PriorPrediction.name);
    set(gcf,'color','w');
end


