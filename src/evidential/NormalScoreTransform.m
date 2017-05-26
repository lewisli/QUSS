function [ Transformed ] = NormalScoreTransform( Untransformed, PlotLevel )
%NORMALSCORETRANFORM Perform normal score transform
%
% Computes a histogram transform to take a variable into a Gaussian one
%
% Parameters :
% 	Untransformed: (NReals x NDim) Variable to be transformed
%   PlotLevel: (int) Controls which plots to generate [optional]
%
% Return :
%	Transformed: The transformed variable (Gaussian)
%
% Author: Lewis Li (lewisli@stanford.edu)
% Date:    Feburary 7th 2016

if (nargin < 2)
    PlotLevel =1;
end
Transformed = Untransformed;

for i=1:size(Untransformed,2)
    h=Untransformed(:,i);
    [f,x]=ecdf(h);
    y=norminv(f(2:end-1),mean(h),std(h));
    x=x(2:end-1);
    F = griddedInterpolant(x,y,'linear');
    Transformed(:,i) =F(h);
    
    [mu, sigma] = normfit(Transformed(:,i));
    
    if (PlotLevel == 1)
        subplot(ceil(size(Untransformed,2) /2),2,i);
        histfit(Transformed(:,i));
        
        format short;
        title(['h^{cg}_' num2str(i)  ' ~ N(' ...
            num2str(mu) ',' num2str(sigma) ')' ] )
    end
end


end

