function basisobj = bspline_basis(rangeval, nbasis, norder, breaks)
%SPLINE_BASIS Creates a bspline functional data basis.
%  This function is identical to CREATE_BSPLINE_BASIS.
%  Arguments ...
%  RANGEVAL ... an array of length 2 containing the lower and upper
%               boundaries for the rangeval of argument values
%  NBASIS   ... the number of basis functions
%  NORDER   ... order of b-splines (one higher than their degree).  The
%                 default of 4 gives cubic splines.
%  BREAKS   ... also called knots, these are a strictly increasing sequence
%               of junction points between piecewise polynomial segments.
%               They must satisfy BREAKS(1) = RANGEVAL(1) and
%               BREAKS(NBREAKS) = RANGEVAL(2), where NBREAKS is the total
%               number of BREAKS.  There must be at least 3 BREAKS.
%  There is a potential for inconsistency among arguments NBASIS, NORDER, and
%  BREAKS.  It is resolved as follows:
%     If BREAKS is supplied, NBREAKS = length(BREAKS), and
%     NBASIS = NBREAKS + NORDER - 2, no matter what value for NBASIS is
%     supplied.
%     If BREAKS is not supplied but NBASIS is, NBREAKS = NBASIS - NORDER + 2,
%        and if this turns out to be less than 3, an error message results.
%     If neither BREAKS nor NBASIS is supplied, NBREAKS is set to 21.
%  Returns
%  BASISOBJ  ... a functional data basis object

%  A B-spline basis may also be constructed using CREATE_EASY_BASIS 
%    CREATE_BASIS or MAKE_BASIS.

%  last modified 29 November 2000

  if nargin < 4
    breaks = 0;
  end
  if nargin < 3
    norder = 4;
  end
  if nargin < 2
    nbasis = 0;
  end
  type = 'bspline';

  if nbasis == 0 & breaks == 0
    nbreaks = 21;
    nbasis  = 19 + norder;
    breaks  = linspace(rangeval(1), rangeval(2), nbreaks);
  end
  if nbasis == 0 & breaks ~= 0
    nbreaks = length(breaks);
    nbasis  = nbreaks + norder - 2;
  end
  if nbasis ~= 0 & breaks == 0
    nbreaks = nbasis - norder + 2;
    breaks  = linspace(rangeval(1), rangeval(2), nbreaks);
  end
  nbreaks = length(breaks);
  
  if (nbreaks < 3)
    error ('Number of values in BREAKS less than 3.');
  end
  if (nbasis < nbreaks)
    error ('NBASIS is less than number of values=BREAKS.');
  end
  if (breaks(1) ~= rangeval(1))
    error('Smallest value in BREAKS not equal to RANGEVAL(1).');
  end
  if (breaks(nbreaks) ~= rangeval(2))
    error('Largest  value in BREAKS not equal to RANGEVAL(2).');
  end
  if (min(diff(breaks)) <= 0)
    error('Values in BREAKS not strictly increasing');
  end

  %  The PARAMS field contains only the interior knots

  params = breaks(2:(nbreaks-1));
  basisobj = basis(type, rangeval, nbasis, params);
