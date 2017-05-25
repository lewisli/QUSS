function temp = eval_Lphi(jvar, bwtcell, basisobj)
%  Evaluates the linear combination of basis function derivatives 
%  for variable JVAR in differential operator L
%  defined by weight cell array BWTCELL.
%  The evaluation is at the quadrature points.

%  Last modified 10 August 2004

Dorder    = size(bwtcell, 2);
nbasis    = getnbasis(basisobj);
onebas    = ones(1,nbasis);

%  the function term in the operator

bfd       = getfd(bwtcell{jvar,1});
bbasis    = getbasis(bfd);
bcoef     = getcoef(bfd);
bbasismat = getvalues(bbasis);
bmat      = (bbasismat*bcoef)*onebas;
basismat  = getvalues(basisobj);
temp      = basismat.*bmat;

%  the terms for the derivatives

for jderiv = 2:Dorder
    Dbfd       = getfd(bwtcell{jvar, jderiv});
    Dbbasis    = getbasis(Dbfd);
    Dbcoef     = getcoef(Dbfd);
    Dbbasismat = getvalues(Dbbasis);
    Dbmat      = (Dbbasismat*Dbcoef)*onebas;
    Dbasismat  = getvalues(basisobj,jderiv-1);
    temp       = temp + Dbasismat.*Dbmat;
end

