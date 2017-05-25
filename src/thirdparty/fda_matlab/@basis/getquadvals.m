function quadvals = getquadvals(basisobj)
%  GETQUADVALS   Extracts the quadrature points and weights
%     from basis object BASISOBJ.

%  last modified 19 May 2004

if ~isa_basis(basisobj)
    error('Argument is not a functional basis object.');
end

quadvals = basisobj.quadvals;
