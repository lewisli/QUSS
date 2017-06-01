close all; clear all;
addpath('../src/evidential/');
addpath('../src/thirdparty/fda_matlab/');

load('../data/SIR/3dslresults/Prior.mat');
load('../data/SIR/3dslresults/Time1000/Observed.mat')

FontSize = 30;
PlotResponses(PriorData,TrueData,FontSize);

load('../data/SIR/parameters/Prior/Sampled.mat');
load('../data/SIR/parameters/ParameterNames.mat');


%% Dimension Reduction on Data Variable
addpath('../src/thirdparty/fda_matlab/');
% Compute FPCA on prior data + observed data
data_FPCA = ComputeHarmonicScores(PriorData,TrueData,0);

% ComputeHarmonicScore stacks the observed data as the last row of the
% resulting scores
obs_realization = size(PriorData.data,1) + 1;

% Perform Mixed PCA
eigentolerance = 0.9;
rmpath('../src/thirdparty/fda_matlab/');
[mpca_scores, mpca_obs] = MixedPCA(data_FPCA,obs_realization,...
    eigentolerance);

%%
addpath('../src/thirdparty/likelihood_continuous/');

NumParameters = length(ParameterNames);

NumPDFPoints = 100;
ProposalParameterVal = zeros(NumPDFPoints,NumParameters);
ProposalPDFs = zeros(NumPDFPoints,NumParameters);

% Stack mpca_scores and mpca_obs to compute probability
mpca_stacked = [mpca_scores; mpca_obs];

for i = 1:NumParameters
    ParameterName = ParameterNames{i};
    ParameterValues = ModelParameterValues.(ParameterName);
    
    % Kernel Density Estimation
    [ProposalParameterVal(:,i),ProposalPDFs(:,i)] = ...
        UpdateProbabilityContinuous(mpca_stacked,ParameterValues);
end

%%
FontSize = 12;
for i = 1:4
    subplot(2,2,i);
    
    ParameterName = ParameterNames{i};
    ParameterValues = ModelParameterValues.(ParameterName);
    
    [prior_f,prior_xi] = ksdensity(ParameterValues);
    
    hold on;
    plot(prior_xi,prior_f,'k','linewidth',2);
    plot(ProposalParameterVal(:,i),ProposalPDFs(:,i),'b','linewidth',2)
    legend('Prior','Proposal'); axis tight;
    set(gca,'fontsize',FontSize); xlabel(ParameterName); ylabel('PDF');
    
end


%%
addpath('../src/sir');
NumProposalModels = 500;
PlotLevel = 0;
ProposalModelParameters = SampleNewProposalModels(ProposalParameterVal, ProposalPDFs,...
    NumProposalModels,ParameterNames,PlotLevel);

%% Load proposal results
load('../data/SIR/3dslresults/Time1000/Observed.mat')
load('../data/SIR/3dslresults/Time1000/Proposal.mat');
load('../data/SIR/3dslresults/Time1000/Prior.mat');

%PlotResponses(ProposalData,TrueData,FontSize);
PlotResponses(ProposalPrediction,[],FontSize);
%%
EigenTolerance = 0.95;
C_D  = 0;
PlotLevel = 0;

% First step is to compute f(h|d_{obs}) using Evidential Learning
[mu_Target,C_Target,DcPrior,HcPrior, dobs_c]=ComputePosteriorPrediction(...
    PriorData, PriorPrediction, TrueData, EigenTolerance,C_D,PlotLevel);

%%
% We next project the proposal models into the canonical domain
[DcProposal,HcProposal] = ComputeCanonicalOnExternalBasis(PriorData,...
    PriorPrediction,ProposalData,ProposalPrediction,EigenTolerance,TrueData);

%%
hold on;
h1=scatter(DcPrior(:,1),HcPrior(:,1),50,[0.5,0.5,0.5],'filled');
h2=scatter(DcProposal(:,1),HcProposal(:,1),'filled');
h3=plot([dobs_c(1), dobs_c(1)], [min(HcPrior(:,1)), max(HcPrior(:,1))],...
    'r','linewidth',2.0);
legend([h1(1), h2(1), h3],'Prior','Proposal','Observed');
axis tight; grid on;
xlabel('d^c_1'); ylabel('h^c_1');
set(gcf,'color','w');

%%
% Compute densities
% Get density of f(h^\star)
[fPrior,xPrior] = ksdensity(HcPrior(:,1));

% Apply kernel density estimation on H_c2 to get density of
% \hat{f}(h^\star|d_obs)
[fProposal,xProposal] = ksdensity(HcProposal(:,1));

% Get density of f(h^\star|d_obs)
NumPointsForDensity = 1000;
xTarget = linspace(min(xPrior),max(xPrior),NumPointsForDensity);
fTarget = normpdf(xTarget,mu_Target(1),(C_Target(1,1)));

%%
CurrentTime=1000;
plot(xPrior,fPrior,'k','LineWidth',3);
plot(xProposal,fProposal,'r','LineWidth',3);
plot(xTarget,fTarget,'b','LineWidth',3);
title(['Distributions After ' ...
    num2str(CurrentTime) ' days']);
h = legend('Prior: $f(h^\star)$',...
    'Proposal: $\hat{f}(h^\star|d_{obs})$',...
    'Target: $f(h^\star|d_{obs})$');
set(h,'Interpreter','latex','FontSize',FontSize,'FontWeight','bold',...
    'Location','best')
xlab = xlabel('$h^\star$'); ylab = ylabel('PDF','FontSize',FontSize);
set(xlab,'Interpreter','latex','FontSize',FontSize);
axis tight; set(gcf,'color','w'); set(gca,'FontSize',FontSize); grid on;

%%

% Compute weights of models used to calculate \hat{f}(h*|d_obs)
[WNorm] = ComputeModelWeights( HcProposal, fProposal, xProposal,...
    fTarget, xTarget );


% Plot models by weights before resampling
PlotModelsByWeight( HcPrior,DcPrior,HcProposal,...
    DcProposal,dobs_c,WNorm,FontSize);

R = [DcProposal(:,1) HcProposal(:,1)];
[RNewSystematic, WNewSystematic, IndexAfterResampling] = ...
    SystematicResampling(R,WNorm);

%%
PlotModelsByWeight(HcPrior,DcPrior,RNewSystematic(:,2),...
    RNewSystematic(:,1),dobs_c,...
    WNewSystematic,FontSize);
colormap jet; caxis([0, 1/length(HcProposal)*1.5]);

title(['Resampled Proposal Models After ' num2str(CurrentTime)...
    ' days'],'FontSize',FontSize);


ResampledModels = unique(IndexAfterResampling);
Weights = WNorm;
%%
hQuantiles = PlotUpdatedQuantiles(PriorPrediction,...
    ProposalPrediction,TruePrediction,HcProposal,Weights,CurrentTime,FontSize);