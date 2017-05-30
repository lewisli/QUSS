 %% This function performs kernel smoothing in nD and evaluates the pdf of 
 %% points X at the locations defined in loc

function pdf=EvaluateDensity(X,sig,loc)


%% Input Parameters
% - X: matrix (N x n) containing coordinates of the N prior models projected in nD space.  
%      Each row is a point, the columns represent the dimension n of the space. 
% - sig: vector of length n (1 value for each dimension) containing the bandwidth
% - loc: point(s) where to estimate the density.  It should have the same number of columns (n) than X

%% Output Parameters
% pdf:  pdf calculated at location loc

% reference: http://en.wikipedia.org/wiki/Kernel_density_estimation


% check dimensions
if size(X,2) ~= size(loc,2)
    error(' X and loc should have the same number of columns \n');
end

% evaluate the density: use of Gaussian kernels.
n = size(loc,1);
rho=0;
f = 1./(2*pi*prod(sig)*sqrt(1-rho.^2));
      
pdf = zeros(size(loc,1),1);
for i=1:size(X,1)
    pdf_tmp = f*exp(-1./2.*(1-rho.^2)*(sum(((bsxfun(@minus, loc, X(i,:))).^2)./(ones(n,1)*sig.^2),2)));
    pdf=pdf+pdf_tmp;
end
pdf = pdf./size(X,1);

end
