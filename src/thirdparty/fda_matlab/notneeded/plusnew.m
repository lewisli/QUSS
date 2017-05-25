function plusfd = plus(fd1, fd2, basisobj)
%  PLUS: Pointwise sum of two functional data objects,
%    the sum of a scalar and a functional data object,
%    or the sum of a vector and a functional data obect
%       where the length of the vector is the same as the
%       number of replications of the object.
%  When both arguments are functional data objects, 
%  they need not have the same bases, 
%  but they must either (1)  have the same number of replicates, or
%  (2) one function must have a single replicate and other multiple 
%  replicates.  In the second case, each function in the multiple
%  replicate object is multiplied by the singleton function in the 
%  other objects.
%  In either case, they must have the same number of functions.
%  When both arguments are functional data objects, the
%  basis used for the sum is constructed to be of higher 
%  dimension than the basis for either factor according to rules
%  described in function TIMES for two basis objects.  

%  last modified 22 October 2003

if ~(isa_fd(fd1) | isa_fd(fd2))
      error('Neither argument for * is a functional data object.');
end

if isa_fd(fd1) & isa_fd(fd2)
    %  both arguments are functional data objects
    coef1  = getcoef(fd1);
    coef2  = getcoef(fd2);
    coefd1 = size(coef1);
    coefd2 = size(coef2);
    ndim1  = length(coefd1);
    ndim2  = length(coefd2);
    if ndim1 ~= ndim2
        error('Dimensions of coefficient matrices not compatible.');
    end
    %  allow for one function being a single replicate,
    %  and if so, copy it as many times as there are replicates
    %  in the other function.
    if coefd1(2) == 1 & coefd2(2) >  1
        if     ndim1 == 2
            coef1 = coef1*ones(1,coefd2(2));
        elseif ndim1 == 3
            temp = zeros(coefd2);
            for j=1:coefd1(3)
                temp(:,:,j) = squeeze(coef1(:,1,j))*ones(1,coefd2(2));
            end
            coef1 = temp;
        else
            error('Dimensions of coefficient matrices not compatible.');
        end
        coefd1 = size(coef1);
        fd1    = putcoef(fd1, coef1);
    end
    if coefd1(2) >  1 & coefd2(2) == 1
        if     ndim2 == 2
            coef2 = coef2*ones(1,coefd1(2));
        elseif ndim1 == 3
            temp = zeros(coefd1);
            for j=1:coefd2(3)
                temp(:,:,j) = squeeze(coef2(:,1,j))*ones(1,coefd1(2));
            end
            coef2 = temp;
        else
            error('Dimensions of coefficient matrices not compatible.');
        end
        coefd2 = size(coef2);
        fd2    = putcoef(fd2, coef2);
    end
    if coefd1(2) ~= coefd2(2) 
        error('Number of replications are not equal.');
    end
    if ndim1 > 2 & ndim2 > 2 & ndim1 ~= ndim2
        error(['Both arguments multivariate, ',  ...
               'but involve different numbers ', ...
               'of functions.']);
    end
    basisobj1 = getbasis(fd1);
    basisobj2 = getbasis(fd2);
    nbasis1   = getnbasis(basisobj1);
    nbasis2   = getnbasis(basisobj2);
    rangeval1 = getbasisrange(basisobj1);
    rangeval2 = getbasisrange(basisobj2);
    if (any(rangeval1 ~= rangeval2))
        error('The ranges of the arguments are not equal.');
    end
    neval     = max(10*max(nbasis1+nbasis2) + 1, 101);
    evalarg   = linspace(rangeval1(1), rangeval2(2), neval)';
    fdarray1  = eval_fd(fd1, evalarg);
    fdarray2  = eval_fd(fd2, evalarg);
    if (ndim1 <= 2 & ndim2 <= 2) | (ndim1 > 2 & ndim2 > 2)
        fdarray = fdarray1 + fdarray2;
    end
    if ndim1 == 2 & ndim2 > 2
        fdarray = zeros(coefd2);
        for ivar = 1:coefd2(3)
          fdarray(:,:,ivar) = fdarray1 + fdarray2(:,:,ivar);
        end
    end
    if ndim1 > 2 & ndim2 == 2
        fdarray = zeros(coefd1);
        for ivar = 1:coefd1(3)
            fdarray(:,:,ivar) = fdarray1(:,:,ivar) + fdarray2;
        end
    end
    %  set up basis for sum
    if nargin < 3, basisobj = basisobj1.*basisobj2;  end
    coefprod = project_basis(fdarray, evalarg, basisobj, 1);
    fdnames1 = fd1.fdnames;
    fdnames2 = fd2.fdnames;
    fdnames  = fdnames1;
    fdnames{3} = [fdnames1{3},'*',fdnames2{3}];
 else
    %  one argument is numeric and the other is functional
    if ~(isnumeric(fd1) | isnumeric(fd2))
        error('Neither argument for * is numeric.');
    end
    if isnumeric(fd1) & isa_fd(fd2)
        fac = fd1;
        fd  = fd2;
    elseif isa_fd(fd1) & isnumeric(fd2)
        fac = fd2;
        fd  = fd1;
    else
        error('One of the arguments for * is of the wrong class.');
    end
    coef     = getcoef(fd);
    coefd    = size(coef);
    N        = coefd(2);
    if ~(length(fac) == 1 | all(size(fac) == [1,N]))
        error('Numerical factor not of correct dimensions');
    end
    basisobj = getbasis(fd);
    nbasis   = getnbasis(basisobj);
    rangeval = getbasisrange(basisobj);
    neval    = max(10*nbasis + 1,101);
    neval    = min(neval,201);
    if all(size(fac) == [1,N])
        fac = ones(neval,1)*fac;
    end
    evalarg  = linspace(rangeval(1),rangeval(2), neval)';
    fdarray  = fac + eval_fd(evalarg, fd);
    coefprod = project_basis(fdarray, evalarg, basisobj);
    fdnames  = fd.fdnames;
    if length(fac) == 1
        fdnames{3} = [num2str(fac),'*',fdnames{3}];
    end
end

plusfd.coef     = coefprod;
plusfd.basisobj = basisobj;
plusfd.fdnames  = fdnames;

plusfd = class(plusfd, 'fd');