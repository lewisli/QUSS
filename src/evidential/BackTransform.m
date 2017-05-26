function [ Output ] = BackTransform(Input, Match)
%BACKTRANFORM Undo Normal Score Transform
%
% Computes a histogram transform on Gaussian variable Input such 
% that the resulting variable given in Output has a histogram that
% matches that of the variable given in Match.
%
% Parameters :
%   Input: (NReal x NDim) Gaussian variable to be transformed
%   Match: (NInputReal x NDim) Variable's who's histogram we are trying to match
%
% Return :
%   Output: (NReal x NDim) The transformed variable
%
% References
%   Isaaks, Edward H., and R. Mohan Srivastava. "Applied geostatistics."
%   (1989): 1-561.
%
% Author: Lewis Li (lewisli@stanford.edu)
% Date:    Feburary 7th 2016

Output = zeros(size(Input));

for i = 1:size(Input,1)

    OriginalScores = Match(:,i);
    TransformedSamples = Input(i,:);
    BackTransformedValue = zeros(size(TransformedSamples));
    
    [F,x] = ecdf(OriginalScores);
    
    % Get the cdf values at each value assuming a gaussian with mean/std
    % given by OriginalScores
    FStar = normcdf(TransformedSamples,mean(OriginalScores),...
        std(OriginalScores));
    
    % for each FStar, find closest F
    for j = 1:length(FStar)
        
        [~,index] = min(abs(F-FStar(j)));
        
        if (index == 1)
            BackTransformedValue(j) = x(index);
        elseif (index == length(x))
            BackTransformedValue(j) = x(end);
        elseif (F(index) < FStar(j))
            BackTransformedValue(j) = 0.5*(x(index)+x(index-1));
        else
            BackTransformedValue(j) = 0.5*(x(index)+x(index+1));
        end
    end

    
    Output(i,:) = BackTransformedValue;
end

end

