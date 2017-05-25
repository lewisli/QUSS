function powermmat = powerbasis(x, exponents, nderiv)
%POWERBASIS computes values of monomials, or their derivatives.
%  The powers of X are the NBASIS nonnegative integers in EXPONENTS.
%  The default is 1, meaning X itself.
%  Arguments are as follows:
%  X         ... vector of values at which the polynomials are to
%                evaluated
%  EXPONENTS ... vector of exponents
%  NDERIV    ... order of derivative to be returned.
%  Return is:
%  A matrix with length(X) rows and NBASIS columns containing
%    the values of the monomials or their derivatives

%  last modified 22 January 2002

n = length(x);
xdim = size(x);
if length(xdim) == 2
    if xdim(2) == n
        x = x';
        xdim = size(x);
    end
    if xdim(1) ~= n | xdim(2) ~= 1
        error('X is not a column vector.');
    end
end

% set default arguments

if nargin < 3, nderiv = 0; end

nbasis = length(exponents);

powermat = zeros(n,nbasis);
if nderiv == 0
    for ibasis=1:nbasis
        powermat(:,ibasis) = x.^exponents(ibasis);
    end
else
    if any(exponents - neriv < 0) & any(x == 0)
        return;
    else
        for ibasis=1:nbasis
            degree = exponents(ibasis);
            if nderiv <= degree
                fac = degree;
                for ideriv=2:nderiv
                    fac = fac*(degree-ideriv+1);
                end
                powermat(:,ibasis) = fac.*x.^(degree-nderiv);
            end
        end
    end
end

