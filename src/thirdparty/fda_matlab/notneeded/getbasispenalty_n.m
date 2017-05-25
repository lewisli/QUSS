function penaltymat = getbasispenalty_n(basisfd, penspec)
%  GETBASISPENALTY   Computes penalty matrix for basis.fd object BASISSTR.
%  Arguments:
%  BASISOBJ ... a basis functional data object
%  PENSPEC  ... A penalty structure. This has elements
%  PEN.LFD    ... The order of derivative or a linear differential
%                 operator to be penalized = the smoothing phase.
%                 By default penspec is set = function GETBASISPENALTY
%  PEN.LAMBDA ... The smoothing parameter determining the weight to be
%                 placed on the size of the derivative = smoothing.  This
%                 is 0 by default.
%  PEN.RNG    ... A vector containing 1 or 2 elements defining a value or
%                 a range of values over which the penalization is defined
%  PEN.WTFD   ... A functional data object specifying a weight function
%                 The default is 1, meaning a constant = 1 basis.                  

%    last modified 26 April 1999

  if nargin < 2
    penspec = 0;
  end

  if ~isa_basis(basisfd)
    error('Argument BASISSTR is not a functional basis object.');
  end

  type   = getbasistype(basisfd);

  switch type
    case 'fourier'
      penaltymat = fourierpen(basisfd, penspec.Lfd);
    case 'bspline'
      penaltymat = bsplinepen_n(basisfd, penspec);
    case 'monom'
      penaltymat = monompen(basisfd, penspec.Lfd);
    case 'polyg'
      penaltymat = polygpen(basisfd, penspec.Lfd);
    case 'expon'
      penaltymat = exponpen(basisfd, penspec.Lfd);
    case 'const'
      penaltymat = basisfd.rangeval(2) - basisfd.rangeval(1);
    otherwise
      error('Basis type not recognizable');
  end

