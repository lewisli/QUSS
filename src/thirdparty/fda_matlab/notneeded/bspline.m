function basismat = bspline (x, breaks, norder, nderiv)
%  BSPLINE  Computes values or derivative values of B-spline basis functions
%  Arguments:
%  X      ... Argument values for which function values are computed
%  BREAKS ... Strictly increasing knot sequence spanning argument range
%  NORDER ... Order of B-spline (one greater than degree) max = 19
%             Default 4.
%  NDERIV ... Order of derivative required, default 0.
%  Calls spline toolbox function SPCOL

%  last modified 15 June 1999

  x = squeeze(x);
  sizex = size(x);
  ndim  = length(sizex);
  switch ndim
    case 2
      if sizex(1) > 1 & sizex(2) > 1
        error('First argument must be a vector');
      else
        n = length(x);
      end
    case 1
      n = length(x);
      x = x';
    otherwise
      error('First argument must be a vector');
  end

  if nargin < 4
    nderiv = 0;
  end

  if nargin < 3
    norder = 4;
  end

  if nargin < 2
    error('Knots must be supplied as the second argument');
  end

  nbreaks = length(breaks);
  if nbreaks < 2, error('Number of knots less than 2.'); end
  if any(diff(breaks) <= 0)
    error ('Knots are not strictly increasing.');
  end
  if norder  < 1 | norder > 19
    error ('Order of basis out of admissible range 1 ... 19');
  end
  if nderiv < 0, error('NDERIV is negative');  end
  if nderiv >= norder
     error ('NDERIV cannot be as large as order of B-spline.');
  end

  nderivp1 = nderiv + 1;

  %  compute basis matrix

  knots  = [breaks(1)      *ones(1,norder-1), breaks, ...
            breaks(nbreaks)*ones(1,norder-1)];
  if nderiv == 0
    tau    = reshape(x,[1,n]);
  else
    onevec = ones(nderivp1,1);
    xmat   = reshape(x,[1,n]);
    tau    = reshape(onevec * xmat,[1,n*nderivp1]);
  end

  basismat = spcol(knots,norder,tau);
  if nderiv > 0
    index = nderivp1:nderivp1:n*nderivp1;
    basismat = basismat(index,:);
  end

