function breaks = quantbrk(x, rangeval, nbasis, norder)
%QUANTBRK sets knots for B-spline basis at quantiles of X
%Arguments:
%  X        ... Argument values 
%  RANGEVAL ... Range over which basis is defined
%  NBASIS   ... Number of basis functions
%  NORDER   ... Order of basis

%  Last modified:  8 December 2000

n       = length(x);
xsort   = sort(x);
indelta = floor(n/(nbasis-norder+1));
brkindx = (1:(nbasis-norder))*indelta;
intbrk  = xsort(brkindx);
breaks  = [rangeval(1),intbrk',rangeval(2)];
