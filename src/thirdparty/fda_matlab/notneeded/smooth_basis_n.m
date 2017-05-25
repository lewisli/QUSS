function smoothlist = smooth_basis_n(y, argvals, basisobj, wtvec, ...
                                     penspecs, fdnames)
%  SMOOTH_BASIS_N  Smooths discrete curve values using penalized basis functions
%  This is a version of SMOOTH_BASIS that inputs a cell whose entries contain
%  a penalty structure.
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
%  ARGVALS  ... A set of argument values, set by default to equally spaced on
%               the unit interval (0,1).
%  BASISOBJ ... A basis.fd object created by function create_basis.fd.
%  WTVEC    ... A vector of N weights, set to one by default, that can
%               be used to differentially weight observations = the
%               smoothing phase
%  PENSPECS ... A cell whose elements are each penalty structures.  
%               A penalty structure has elements
%  PEN.LFD    ... The order of derivative or a linear differential
%                 operator to be penalized = the smoothing phase.
%                 By default Lfd is set = function GETBASISPENALTY
%  PEN.LAMBDA ... The smoothing parameter determining the weight to be
%                 placed on the size of the derivative = smoothing.  This
%                 is 0 by default.
%  PEN.RNG    ... A vector containing 1 or 2 elements defining a value or
%                 a range of values over which the penalization is defined
%  PEN.WTFD   ... A functional data object specifying a weight function
%                 The default is 1, meaning a constant = 1 basis.                  
%  FDNAMES  ... A cell of length 3 with names for
%               1. argument domain, such as 'Time'
%               2. replications or cases
%               3. the function.
%  Returns a cell object containing:
%  FD    ...  an object of class fd containing coefficients
%  DF    ...  a degrees of freedom measure
%  GCV   ...  a measure of lack of fit discounted for df.

%  last modified 14 January 2003

  n = length(argvals);

  if nargin < 6
    fdnames{1} = 'time';
    fdnames{2} = 'reps';
    fdnames{3} = 'values';
  end

 if nargin < 5
    penspecs = {0};
  end

  if nargin < 4
    wtvec = ones(n,1);
  end

  if ~iscell(penspecs)
    error('Argument PENSPECS is not a cell.')
  end
  
  sizew = size(wtvec);
  if (length(sizew) > 1 & sizew(1) > 1 & sizew(2) > 1) | ...
      length(sizew) > 2
    error ('WTVEC must be a vector.');
  end
  if length(sizew) == 2 & sizew(1) == 1
    wtvec = wtvec';
  end
  if length(wtvec) ~= n
    error('WTVEC of wrong length');
  end
  if min(wtvec) <= 0
    error('All values of WTVEC must be positive.');
  end

  ydim = size(y);
  ndim  = length(ydim);

  if ydim(1) ~= n
    error('Number of arguments differs from first dimension of matrix of values');
  end

  switch ndim
    case 1
      ncurves = 1;
      nvar    = 1;
    case 2
      ncurves = ydim(2);
      nvar    = 1;
    case 3
      ncurves = ydim(2);
      nvar    = ydim(3);
    otherwise
      error('Second argument must not have more than 3 dimensions');
  end

  nbasis   = getnbasis(basisobj);
  onebasis = ones(1,nbasis);

  basismat = getbasismatrix(argvals, basisobj);

  nobasis = isnumeric(penspecs{1});
  if nobasis
    npenspec = 0;
    lambda = 0;
  else
    penspecd = size(penspecs);
    if length(penspecd) > 2 | (penspecd(1) > 1 & penspecd(2) > 1)
      error('Argument PENSPECS must be a vector of cells');
    end
    npenspec = penspecd(1);
    penspec1 = penspecs{1};
  end
  
    %  The following code is for the coefficients completely determined

    basisw = basismat .* (wtvec * ones(1,nbasis));
    Bmat = basisw' * basismat;

    Cmat = Bmat;
    %  smoothing required, set up coefficient matrix for normal equations
    for i=1:npenspec
      penspeci = penspecs{i};
      penmat = eval_penalty_n(basisobj, penspeci);
      Cmat = Cmat + penspeci.lambda.*penmat;
    end
    Cmat = (Cmat + Cmat')./2;
    
   %  compute inverse of Cmat

    if is_diag(Cmat)
      Cmatinv = diag(1/diag(Cmat));
    else
      Lmat    = chol(Cmat);
      Lmatinv = inv(Lmat);
      Cmatinv = Lmatinv * Lmatinv';
    end

    %  compute degrees of freedom of smooth

    df = sum(diag(Cmatinv * Bmat));

    %  solve normal equations for each observation

    if ndim < 3
      Dmat = basisw' * y;
      coef = Cmatinv * Dmat;
    else
      for ivar = 1:nvar
        Dmat = basisw' * y(:,:,ivar);
        coef(:,:,ivar) = Cmatinv * Dmat;
      end
    end

  %  compute  GCV index

  if df < n
    if ndim < 3
      yhat = basismat * coef;
      SSE = sum((y - yhat).^2);
    else
      SSE = 0;
      for ivar = 1:nvar
        yhat = basismat * coef(:,:,ivar);
        SSE = SSE + sum((y(:,:,ivar) - yhat).^2);
      end
    end
    gcv = (SSE/n)/(nvar*(n - df)/n)^2;
  else
    gcv = NaN;
  end

  fdobj = fd(coef, basisobj, fdnames);

  smoothlist.fdobj = fdobj;
  smoothlist.df    = df;
  smoothlist.gcv   = gcv;



