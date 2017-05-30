function [SampleForParam,pdfParam] = UpdateProbabilityContinuous(Y,ParamValues)

%% Inputs:
    % - Y: coordinates in MDS space. Last line correspond to the observed data
    % - ParamValues: values of the continuous parameter W
    % - prior:  if provided, the prior distribution is provided as well

%% Outputs:
    % - SampleForParam: vector where the PDF is sampled
    % - pdfParam: value of the PDF
    
    SampleForParam = linspace(min(ParamValues),max(ParamValues),100);
    PDFestimation = horzcat(repmat(Y(end,:),length(SampleForParam),1),SampleForParam');

    % evaluate the bandwith for D and W 
    [sigDataProp,~] = GetBandwidthDataAndProp(Y(1:end-1,:),ParamValues);

    % f(w|dobs)
    pdfParam = EvaluateDensity(horzcat(Y(1:end-1,:),ParamValues),sigDataProp,PDFestimation); 
    pdfParam = pdfParam ./ trapz(SampleForParam,pdfParam);
    
%     figure; 
%     plot(SampleForParam,pdfParam,'-b','LineWidth',2);xlim([min(SampleForParam),max(SampleForParam)]); hold on
%     ylabel('PDF','FontSize',10,'FontWeight','b')

end