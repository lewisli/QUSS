function basismat = bsplineFDA(x, breaks, norder, nderiv)
%  BSPLINEFDA  Computes values or derivative values of B-spline basis functions
%  Arguments:
%  X      ... Argument values for which function values are computed
%  BREAKS ... Strictly increasing knot sequence spanning argument range
%  NORDER ... Order of B-spline (one greater than degree) max = 19
%             Default 4.
%  NDERIV ... Order of derivative required, default 0.
%  This version calls function BsplineM.  
%  Earlier version called spline toolbox function SPCOL

%  last modified 14 March 2001

%  Note:  the earlier version was called bspline, but this conflicted with
%    a function already available in the spline toolbox.

  x     = squeeze(x);
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
  
  if min(diff(x)) < 0 
      [x,isrt] = sort(x); 
      sortwrd = 1;
  else
      sortwrd = 0;
  end
  if x(1) - breaks(1) < -1e-10 | x(n) - breaks(nbreaks) > 1e-10
    disp([x(1), x(n)])
    error ('Argument values out of range.')
  end
  tau    = reshape(x,[1,n]);

  basismat = BsplineM(tau, breaks, norder, nderiv, 1);
  
  if sortwrd, basismat(isrt,:) = basismat;  end
   
