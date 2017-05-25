function XtWY = trapzmat(X,Y,delta,wt)
%TRAPZMAT integrates the products of two matrices of values
%   using the trapezoidal rule, assuming equal spacing
%  X is the first  matrix of values
%  Y is the second matrix of values
%  DELTA is the spacing between argument values (one by default)
%  WT is a vector of weights (ones by default)
%
%  XtWY is a matrix of integral estimates, number of rows equal to
%  number of col of X, number of cols equal to number of cols of Y

n = size(X,1);

if size(Y,1) ~= n
    error('X and Y do not have same number of rows.');
end

%  set default arguments

if nargin < 4, wt = ones(n,1); end
if nargin < 3, delta = 1; end

if size(wt,1) ~= n
    error('X and WT do not have same number of rows.');
end

if size(wt,2) ~= 1
    error('WT is not a column vector');
end

if delta <= 0
    error('DELTA is not a positive value.');
end

wt([1,n]) = wt([1,n])/2;
wt = wt.*delta;

X = X.*(wt*ones(1,size(X,2)));
XtWY = X'*Y;


