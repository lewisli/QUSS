function fdcell = fd2cell(fdobj)
%FD2CELL converts a univariate functional data object to a cell
%  object, mainly for purposes of defining a linear differential
%  operator object where the arguments are required to be cells.

%  Last modified 5 November 2003

%  check FDOBJ

if ~isa_fd(fdobj)
    error('FDOBJ is not a functional data object.');
end

%  get the coefficient matrix and the basis

coef     = getcoef(fdobj);
coefsize = size(coef);

%  check wether FDOBJ is univariate

if length(coefsize) > 2
    error('FDOBJ is not univariate.');
end

for i=1:coefsize(2)
    fdcell{i} = fdPar(fdobj(i));
end
