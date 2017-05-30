function [IndexAfterResampling,WNorm,N_eff,hFigures] = ...
    CalculateModelWeightsAndResample(DcTarget, HcTarget, ...
    DcProposal, HcProposal,dobs_c,muTarget, CTarget, CurrentTime,PlotLevel)
%CALCULATEMODELWEIGHTS Summary of this function goes here
%   Detailed explanation goes here

FontSize = 40;
MarkerSize = 40;

%hFigures = [];

% Plot the proposal models on the same support as the prior models along
% with the actual observed data
if (PlotLevel >=1)
    
    MarkerSize = 150;
    FontSize=FontSize+4;
    hFigures(1) = figure('Position', [100, 100, 1000, 800]);
    hold on;
    h1 = scatter(DcTarget(:,1),HcTarget(:,1),MarkerSize,'k','filled');
    h2 = scatter(DcProposal(:,1),HcProposal(:,1),MarkerSize-50,'r','filled');
    set(gcf,'color','w'); set(gca,'FontSize',FontSize-4);
    xlabel('d_c^1','FontSize',FontSize);
    ylabel('h_c^1','FontSize',FontSize);
    grid on;
    plot([dobs_c(1) dobs_c(1)],[min(HcTarget(:,1)) ...
        max(HcTarget(:,1))],'b','LineWidth',2);
    tt1=text(dobs_c(1)+0.05,min(HcTarget(:,1))+0.5,'d_{obs}',...
        'Fontweight','b','FontSize',FontSize);
    title(['Prior and Proposal Models After ' ...
        num2str(CurrentTime) ' days'],'FontSize',FontSize);
    hlegend = legend([h1(1), h2(1)], 'Prior', 'Proposal');
    set(hlegend,'Location','southeast');
    set(tt1,'color','b','FontSize',FontSize);
    axis tight; axis square;
%     
%     hFigures(12) = figure('Position', [100, 100, 1000, 800]);
% 
%     hold on;
%     h1 = scatter(DcTarget(:,1),HcTarget(:,1),MarkerSize,[0.5, 0.5, 0.5],'filled');
%     set(gcf,'color','w'); set(gca,'FontSize',FontSize-4);
%     xlabel('d_c^1','FontSize',FontSize);
%     ylabel('h_c^1','FontSize',FontSize);
%     grid on;
%     plot([dobs_c(1) dobs_c(1)],[min(HcTarget(:,1)) ...
%         max(HcTarget(:,1))],'b','LineWidth',2);
%     tt1=text(dobs_c(1)+0.25,min(HcTarget(:,1))+1,'d_{obs}',...
%         'Fontweight','b','FontSize',FontSize);
%     hlegend = legend([h1(1)], 'Prior Models');
%     set(hlegend,'Location','southeast');
%     set(tt1,'color','b','FontSize',FontSize);
%     axis tight; axis square;
%     
%     yl = ylim;
%     
%     hFigure(13) = figure('Position', [100, 100, 200, 800]);
% 
%     x = [min(HcTarget(:,1)):.1:max(HcTarget(:,1))];
%     norm = normpdf(x,mean(HcTarget(:,1)),std(HcTarget(:,1)));
%     norm2 = normpdf(x,mean(HcProposal(:,1)),std(HcProposal(:,1)));
%     hold on;
%     plot(norm, x, 'color', [0.5, 0.5, 0.5], 'linewidth',4);
%     plot(norm2, x, 'color', 'b', 'linewidth',4);
%     set(gca,'XDir','reverse')
%     ylim(yl);
%     grid on;
%     legend('Prior', 'Posterior');
%     set(gcf,'color','w'); set(gca,'FontSize',FontSize-4);
end

% Compute densities
% Get density of f(h^\star)
[fPrior,xPrior] = ksdensity(HcTarget(:,1));

% Apply kernel density estimation on H_c2 to get density of
% \hat{f}(h^\star|d_obs)
[fProposal,xProposal] = ksdensity(HcProposal(:,1));

% Get density of f(h^\star|d_obs)
NumPointsForDensity = 1000;
xTarget = linspace(min(xPrior),max(xPrior),NumPointsForDensity);
fTarget = normpdf(xTarget,muTarget(1),(CTarget(1,1)));

if (PlotLevel >=2)
    hFigures(2) = figure('Position', [100, 100, 1000, 800]);
    hold on;
    plot(xPrior,fPrior,'k','LineWidth',3);
    plot(xProposal,fProposal,'r','LineWidth',3);
    plot(xTarget,fTarget,'b','LineWidth',3);
    title(['Distributions After ' ...
        num2str(CurrentTime) ' days']);
    h = legend('Prior: $f(h^\star)$',...
        'Proposal: $\hat{f}(h^\star|d_{obs})$',...
        'Target: $f(h^\star|d_{obs})$');
    set(h,'Interpreter','latex','FontSize',FontSize,...
        'FontWeight','bold','Location','best')
    xlab = xlabel('$h^\star$');
    ylab = ylabel('PDF','FontSize',FontSize);
    set(xlab,'Interpreter','latex','FontSize',FontSize);
    axis tight;
    set(gcf,'color','w');
    set(gca,'FontSize',FontSize-4);
    grid on;
end

% Compute weights of models used to calculate \hat{f}(h*|d_obs)
[WNorm] = ComputeModelWeights( HcProposal, fProposal, xProposal,...
    fTarget, xTarget );

if (PlotLevel >=3)
    % Plot models by weights before resampling
    hFigures(3) = PlotModelsByWeight( HcTarget,DcTarget,HcProposal,...
        DcProposal,dobs_c,WNorm,FontSize);
%     title(['Proposal Models After ' num2str(CurrentTime)...
%         ' days'],'FontSize',FontSize);
end

% Perform resampling only when effective model number drops beneath a
% certain threshold
N_eff = 1/sum(WNorm.^2);

R = [DcProposal(:,1) HcProposal(:,1)];
[RNewSystematic, WNewSystematic, IndexAfterResampling] = ...
    SystematicResampling(R,WNorm);

if (PlotLevel >=4)
    % Plot models by weights after resampling
    hFigures(4) = PlotModelsByWeight(HcTarget,DcTarget,RNewSystematic(:,2),...
        RNewSystematic(:,1),dobs_c,...
        WNewSystematic,FontSize);
    colormap jet; caxis([0, 1/length(HcProposal)*1.5]);
    
    title(['Resampled Proposal Models After ' num2str(CurrentTime)...
        ' days'],'FontSize',FontSize);
end





end

