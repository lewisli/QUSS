function  coef = project_basis(y, argvals, basismat, penmat)
% PROJECT_BASISMAT Project discrete values on to a set of basis fn. values
%  Arguments for this function:
%
%  Y        ... an array containing values of curves
%               If the array is a matrix, rows must correspond to argument
%               values and columns to replications, and it will be assumed
%               that there is only one variable per observation.
%               If Y is a three-dimensional array, the first dimension
%               corresponds to argument values, the second to replications,
%               and the third to variables within replications.
%               If Y is a vector, only one replicate and variable are assumed.
%  ARGVALS  ... A vector of argument values.  This must be of length
%               length(Y) if Y is a vector or size(Y)(1) otherwise.
%  BASISMAT  ... Matrix of values of basis functions at ARGVALS
%  PENMAT    ... A penalty matrix
%
%  Returns a coefficient vector or array. The first dimension is the number
%     of basis functions and the other dimensions (if any) match
%     the other dimensions of Y.
%

%  Last modified:  18 March 2002

  if nargin < 4
    penalize = 0;
  end

%  Calculate the basis and penalty matrices, using the default
%   for the number of derivatives = the penalty.
  Bmat     = basismat' * basismat + penmat;
  nbasis   = size(basismat,2);
  %  Do the fitting by a simple solution of the
  %    equations taking into account smoothing
  ydim = size(y);
  if(length(ydim) <= 2)
    Dmat = basismat' * y;
    coef = symsolve(Bmat,Dmat);
  else
    nvar = ydim(3);
    coef = zeros([nbasis, ydim(2), nvar]);
    for ivar = 1:nvar
      Dmat = basismat' * y(:,:,ivar);
      coef(:,:,ivar) = symsolve(Bmat,Dmat);
    end
  end

