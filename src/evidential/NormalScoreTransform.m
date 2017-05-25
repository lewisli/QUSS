% function Output = NormalScoreTransform(Input)
% Author: Lewis Li (lewisli@stanford.edu)
% Date:    Feburary 7th 2016
%
% Function   : NormalScoreTransform
%
% Description:  Computes a histogram transform to take a variable into a
% Gaussian one
%
% Parameters :
% 	Input 	   		Variable to be transformed
%
% Return     :
%	Output:	The transformed variable

function [ Output ] = NormalScoreTransform( Input, PlotLevel )

if (nargin < 2)
    PlotLevel =1;
end
Output = Input;

for i=1:size(Input,2)
    h=Input(:,i);
    [f,x]=ecdf(h);
    y=norminv(f(2:end-1),mean(h),std(h));
    x=x(2:end-1);
    F = griddedInterpolant(x,y,'linear');
    Output(:,i) =F(h);
    
    
    
    [mu, sigma] = normfit(Output(:,i));
    
    if (PlotLevel == 1)
        subplot(ceil(size(Input,2) /2),2,i);
        histfit(Output(:,i));
        
        format short;
        title(['h^{cg}_' num2str(i)  ' ~ N(' ...
            num2str(mu) ',' num2str(sigma) ')' ] )
    end
end


end

