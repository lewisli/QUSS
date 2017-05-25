function penaltymat = polynompen(basis, Lfd)
%POLYNOMPEN computes the monomial penalty matrix for penalty LFD.
%  Arguments:
%  BASISFD ... a monomial basis object
%  Lfd     ... a linear differential operator object.
%  Returns the penalty matrix.

%  Last modified:  2 January 2003

if ~isa_basis(basis)
    error('First argument is not a basis object.');
end

if ~isa_Lfd(Lfd)
    error (['Argument Lfd is neither a functional data object', ...
             ' nor an integer.']);
end
  
type = getbasistype(basis);
if ~strcmp(type, 'poly')
    error('basis not of type POLY');
end

%  set up default linear differential operator
if nargin < 2, 
    range = getbasisrange(basis);
    Lfd = Lfd(2, fd(zeros(1,2), create_constant_basis(range))); 
end

if ~isa_Lfd(Lfd)
    error ('Argument Lfd is not a linear differential operator object');
end

nderiv = getnderiv(Lfd);

if nderiv < 0, error('NDERIV is negative'); end

ctr    = getbasispar(basis);
nbasis = getnbasis(basis);
ndegre = nbasis - 1;
exponents = 0:ndegre;

if isinteger(Lfd)
    nbasis = getnbasis(basis);
    penaltymat = zeros(nbasis);
    xrange = getbasisrange(basis);
    for ibasis=1:nbasis
        ideg = exponents(ibasis);
        if nderiv == 0
            ifac = 1;
        else
            ifac = ideg;
            for k=2:nderiv
              ifac = ifac*(ideg - k + 1);
            end
        end
        for jbasis=1:ibasis
            jdeg = exponents(jbasis);
            if nderiv == 0
                jfac = 1;
            else
                jfac = jdeg;
                for k=2:nderiv
                    jfac = jfac*(jdeg - k + 1);
                end
            end
            if ideg >= nderiv & jdeg >= nderiv
                penaltymat(ibasis,jbasis) = ifac*jfac* ...
                    ((xrange(2)-ctr)^(ideg+jdeg-2*nderiv+1) -  ...
                    (xrange(1)-ctr)^(ideg+jdeg-2*nderiv+1));
                penaltymat(jbasis,ibasis) = penaltymat(ibasis,jbasis);
            end
        end
    end
else
    penaltymat = inprod(basis, basis, Lfd, Lfd);
end
