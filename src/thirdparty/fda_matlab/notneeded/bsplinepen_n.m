function penaltymat = bsplinepen_n(basisobj, penspec)
%  BSPLINEPEN    Computes the bspline penalty matrix.
%  Arguments:
%  BASISFD  ... a basis.fd object
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
%  Returns PENALTYMAT the penalty matrix.

%  Last modified:  28 March 2000

  if nargin < 2
    penspec = 0;
  end

  if ~isa_basis(basisobj)
    error('First argument is not a basis object.');
  end

  type = getbasistype(basisobj);
  if ~strcmp(type, 'bspline')
    error('BASISOBJ not of type bspline');
  end

  Lfd = penspec.Lfd;
  [nderiv, derivcoef] = Lset(Lfd);
  
  rng = penspec.rng;
  if all(rng == 0)
    rng = getbasisrange(basisobj);
  end
  
  if nderiv < 0
    error('NDERIV is negative');
  end

  nbasis = getnbasis(basisobj);
  params = getbasispar(basisobj);
  norder = nbasis - length(params);

  if nderiv >= norder
      disp(['Derivative of order', num2str(nderiv), ...
            'cannot be taken for B-spline of order', num2str(norder)]);
      disp('Probable cause is a value of the nbasis argument');
      disp(' in function create.basis.fd that is too small.');
  end

  if length(rng) == 1
    Lbasismat = eval_fd(basisfd, rng, Lfd);
    penaltymat = Lbasismat' * Lbasismat;
  else
    if strcmp(class(Lfd),'double') & round(Lfd) == norder - 1
      rangeval = getbasisrange(basisobj);
      breakvals  = [rangeval(1), params, rangeval(2)];
      nbreakvals = length(breakvals);
      if nderiv >= norder
        error ('NDERIV cannot be as large as order of B-spline.');
      end
      halfseq = (breakvals(2:nbreakvals) + ...
                 breakvals(1:(nbreakvals-1)))./2;
      halfmat = bsplineFDA(halfseq, breakvals, norder, nderiv);
      brwidth = diff(breakvals);
      penaltymat = halfmat' * diag(brwidth) * halfmat;
    else
      penaltymat = inprod_n(basisobj, basisobj, Lfd, Lfd, penspec.wtfd);
    end
  end

