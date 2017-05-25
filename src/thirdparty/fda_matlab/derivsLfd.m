function dy = derivslfd(tnow, y) 
% DERIVS sets up the 1st order system corresponding to   
%   linear differential operator defined by wfd.

%  last modified 22 August 2003

global wfdcell;   

norder = length(wfdcell);
wmat   = zeros(norder, norder);
for j=1:(norder-1)
    wmat(j,j+1) = 1;
end

if norder == 1
    dy = eval_fd(tnow, -wfdcell{1})*y;
else
    for j=1:norder
        wj = eval_fd(tnow, wfdcell{j});
        wmat(norder,j) = -wj;
    end
end
    
dy = wmat * y;

