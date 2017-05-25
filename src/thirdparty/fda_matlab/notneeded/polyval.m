function polymat = polyval(x, norder, nderiv)
%  POLYVAL  Values of the monomials, or their derivatives.
%  NORDER is the number of terms in the highest order polynomial.
%  The degree of the highest order polynomial is one less than NORDER.
%  The default is the linear function.
%  Arguments are as follows:
%  X      ... array of values at which the polynomials are to
%             evaluated
%  NORDER ... order of spline (1 more than degree), so that 1 gives a
%             step function, 2 gives triangle functions,
%             and 4 gives cubic splines
%  NDERIV ... highest order derivative.  0 means only function values
%             are returned.
%  Return is:
%              a matrix with length(X) rows and NORDER columns containing
%  the values of the monomials

%  last modified 15 March 1999

  % set default arguments

  if nargin < 3
    nderiv = 0;
  end

  if nargin < 2
    norder = 2;
  end

  n = length(x);

  if nderiv >= norder
    polymat = zeros(n,norder);
    return
  end

  n = length(x);
  ndegree = norder - 1;

  if nderiv == 0
    %  use the recursion formula to compute polynomial values
    polymat = ones(n,norder);
    if ndegree > 0
      polymat(:,2) = x;
      if ndegree > 1
        for j = 2:ndegree
          polymat(:,j+1) = x.*polymat(:,j);
        end
      end
    end
  else
    polymat = zeros(n,norder);
    fac = 1;
    for j=1:nderiv, fac=fac*j; end
    polymat(:,nderiv+1) = fac;
    for j=(nderiv+2):norder
      fac = (j-nderiv);
      for k=2:nderiv, fac = fac*(j-nderiv+k-1); end
      polymat(:,j) = fac.*x.^(j-nderiv-1);
    end
  end

