function [ffd, wfd, resfd] = pdasystem(fdobj, difeorder, ...
                          fbasis, flambda, ffd0, festimate, ...
                          wbasis, wlambda, wfd0, westimate, n)
%PDA computes the basis function expansion of the
%  estimate of the coefficient functions w_j(t) and forcing function f(t)
%  in the nonhomogeneous linear differential operator
%
%    Lx(t) = f(t) +
%            w_0(t)x(t) + w_1(t)Dx(t) + ... + w_{m-1}D^{m-1}x(t) + 
%            D^m x(t)
%
%  of order m = DIFEORDER that minimizes in a least squares sense the residual
%  functions Lx(t).  
%  The functions x(t) are in functional data object FDOBJ.
%  The coefficient functions are expanded in terms of the
%  basis functions specified in FBASIS and WBASIS.
%
%  Unlike PDASCALAR, the functions x can be vector valued. 
%  That is, the differential equation can be a system of K equations rather
%    than a single equation.   
%  Each coefficient function w_j(t) is matrix valued, with a column
%    for each scalar function in the system.  
%  The forcing function f(t) is vector valued.

%  Arguments:
%  FDOBJ     ...  functional data object
%  DIFEORDER ...  order of the linear differential operator, that is, 
%                 the order of the highest derivative.
%  FBASIS    ...  basis object for the forcing functions.
%  FLAMBDA   ...  the parameter for penalizing the departure of the
%                 estimated forcing function from those defined in FFD0.
%  FFD0      ...  A specification of a functional data object that is used 
%                 for the forcing function if estimated,
%                 or as a target function toward which the estimated 
%                 forcing function is smoothed, or, if the forcing 
%                 functions are not estimated, the fixed forcing
%                 functions.
%                 FFD0 can be:
%                 1.  a constant, in which case the forcing function
%                     functional data object has all coefficients
%                     equal to this constant uses FBASIS as its basis.
%                 2.  a coefficient matrix, used to define the forcing
%                     functional data object.  The coefficient matrix
%                     must be compatible with the structure of FBASIS.
%                 3.  a functional data object with the same structure as 
%                     FFD that is returned by this function.
%  FESTIMATE ...  logical vector or length DIFEORDER, 
%                 if a value is T, the corresponding forcing function  
%                 is estimated, otherwise the target value is used.  
%  WBASIS    ...  basis object for the weight functions.
%  WLAMBDA   ...  a numerical array of length difeorder containing
%                 penalty parameters for penalizing the departure of the
%                 estimated weight functions from those defined in WFD0.
%  WFD0      ...  A specification of a functional data object that is used 
%                 for the weight functions if estimated,
%                 or as a target function toward which the estimated 
%                 weight functions are smoothed, or, if the weight 
%                 functions are not estimated, the fixed weight
%                 functions.
%                 WFD0 can be:
%                 1.  a constant, in which case the weight function
%                     functional data object has all coefficients
%                     equal to this constant uses WBASIS as its basis.
%                 2.  a coefficient matrix, used to define the weight
%                     functional data object.  The coefficient matrix
%                     must be compatible with the structure of WBASIS.
%                 3.  a functional data object with the same structure as 
%                     WFD that is returned by this function.
%  WESTIMATE ...  matrix with NVAR*DIFEORDER rows and  DIFEORDER columns, 
%                 if a value is T, the corresponding coefficient function 
%                 is estimated, otherwise the target value is used.  
%  N         ...  number of sampling points for numerical integration

%  Returns:
%  FFD       ...  FD object for estimated forcing function.  
%  WFD       ...  FD object for estimated weight  functions.  
%                 It has DIFEORDER*NVAR replicates and NVAR functions.
%  RESFD     ...  FD object for residual functions.

%  last modified 14 January 2003

if nargin < 3
    error('There are less than three arguments.');
end

if ~isa_fd(fdobj)
    error('FDOBJ is not a functional data object.');
end

%  get the dimensions of the data in FDOBJ

coef   = getcoef(fdobj);
coefd  = size(coef);
ndim   = length(coefd);
ncurve = coefd(2);
nvar   = 1;
if ndim == 3, nvar = coefd(3);  end
if nvar == 1
    error(['PDAN is designed for a system of more than one equation.', ...
           'Use PDA instead.']);
end
nwfun = nvar*difeorder;

difeorder = floor(difeorder);
if difeorder < 1, error('difeorder is less than one.');  end
nordp1 = difeorder + 1;

%  get information about data basis

basis    = getbasis(fdobj);
nbasis   = getnbasis(basis);
rangeval = getbasisrange(basis);

%  set default arguments

if nargin < 11, n = 5*nbasis; end

if nargin < 10, westimate =  ones(nwfun,nvar); end
if nargin <  9, wfd0      = 0; end
if nargin <  8, wlambda   = zeros(nwfun,nvar); end
if nargin <  7, wbasis    = create_constant_basis(rangeval); end
if nargin <  6, festimate = zeros(nvar,1); end
if nargin <  5, ffd0      = 0; end
if nargin <  4, flambda   = 0; end
if nargin <  3, fbasis    = create_constant_basis(rangeval); end

%  get the characteristics of FBASIS and WBASIS

if ~isa_basis(fbasis)
    error('FBASIS is not a basis object.');
end
ftype   = getbasistype(fbasis);
fnbasis = getnbasis(fbasis);
frange  = getbasisrange(fbasis);

if ~isa_basis(wbasis)
    error('WBASIS is not a basis object.');
end
wtype   = getbasistype(wbasis);
wnbasis = getnbasis(wbasis);
wrange  = getbasisrange(wbasis);

%  check array FLAMBDA

if size(flambda,1) ~= nvar | size(flambda,2) ~= 1
    error('FLAMBDA is not a vector of length NVAR.');
end
if ~all(flambda >= 0)
    error('FLAMBDA does not contain only nonegative values.');
end

%  check array FESTIMATE

if size(festimate,1) ~= nvar | size(festimate,2) ~= 1
    error('FESTIMATE is not a vector of length NVAR.');
end
if ~all(festimate == 0 | festimate == 1)
    error('FESTIMATE does not contain either zero or one.');
end

%  check FFD0 and make into a FD object if numeric

if ~strcmp(class(ffd0), 'fd') & ~isnumeric(ffd0)
    error('FFD0 is neither a vector nor a FD object');
else
    if strcmp(class(ffd0), 'fd')
        ffd0size = size(getcoef(ffd0));
        if length(ffd0size) ~= 2
            error(['Coefficient matrix for FFD0 is not ', ...
                'a two-dimensional array.']);
        end
        if ffd0size(1) ~= fnbasis | ffd0size(2) ~= nvar  
            error(['Coefficient matrix for FFD0 has ', ...
                   'incorrect dimension size.']);
        end
    end
end
if ~strcmp(class(ffd0), 'fd') 
    if length(ffd0) == 1
        ffd0 = ffd0.*ones(fnbasis,nvar);
    end
    if size(ffd0,1) ~= fnbasis | size(ffd0,2) ~= nvar  
        error('FFD0 has incorrect dimensions.');
    end
    ffd0 = fd(ffd0, fbasis);
end

%  check array WLAMBDA

if size(wlambda,1) ~= nwfun | size(wlambda,2) ~= nvar
    error('WLAMBDA does not have the correct dimensions.');
end
if ~all(all(wlambda >= 0))
    error('WLAMBDA does not contain only nonegative values.');
end

%  check array WESTIMATE

if size(westimate,1) ~= nwfun | size(westimate,2) ~= nvar
    error('WESTIMATE does not have the correct dimensions.');
end
if ~all(all(westimate == 0 | westimate == 1))
    error('WESTIMATE does not contain only zeros or ones.');
end

%  check WFD0 and make into a FD object if numeric

if ~strcmp(class(wfd0), 'fd') & ~isnumeric(wfd0)
    error('WFD0 is neither a vector nor a FD object');
else
    if strcmp(class(wfd0), 'fd')
        wfd0size = size(getcoef(wfd0));
        if length(wfd0size) ~= 3
            error(['Coefficient matrix for WFD0 is not ', ...
                   'a three-dimensional array.']);
        end
        if wfd0size(1) ~= wnbasis | ...
           wfd0size(2) ~= nwfun   | ...
           wfd0size(3) ~= nvar
            error(['Coefficient matrix for WFD0 has ',...
                   'incorrect dimension size.']);
        end
    end
end
if ~strcmp(class(wfd0), 'fd') 
    if length(wfd0) == 1
        wfd0 = wfd0.*ones(wnbasis,nwfun,nvar);
    end
    wfd0size = size(wfd0);
    if length(wfd0size) ~= 3
        error(['A numerical WFD0 is neither a constant ',...
               'nor a three-dimensional array.']);
    end
    if wfd0size(1) ~= wnbasis | ...
       wfd0size(2) ~= nwfun   | ...
       wfd0size(3) ~= nvar
        error('WFD0 has incorrect dimensions.');
    end
    wfd0 = fd(wfd0, wbasis);
end

nfcoef = sum(festimate);
nwcoef = sum(sum(westimate));
ncoef  = nfcoef + nwcoef;

if     frange ~= rangeval
    error('Forcing function range not equal to FDOBJ range');
end
if any(wrange ~= rangeval)
    error('Weight function range not equal to FDOBJ range');
end

%  set up sampling values to be used in numerical integration
%    and set up matrix of basis values

delta     = (wrange(2)-wrange(1))/(n-1);
x         = wrange(1):delta:wrange(2);
fbasismat = getbasismatrix(x, fbasis);
wbasismat = getbasismatrix(x, wbasis);

%  set up array to hold values of functions and their derivatives
%   for each variable

yarray = zeros([n,ncurve,nvar,nordp1]);
for ivar=1:nvar
    for j=1:nordp1
        yarray(:,:,ivar,j) = eval_fd(fdobj(:,ivar), x, j-1);
    end
end

%  set up array YPROD to hold mean of products of values in YARRAY

mmat  = m2ij(nordp1,nvar);

yprod = zeros([n,nvar,nordp1,nvar,nordp1]);
for m1=1:nvar*nordp1
    i1 = mmat(m1,2);
    j1 = mmat(m1,1);
    for m2=1:m1;
        i2 = mmat(m2,2);
        j2 = mmat(m2,1);
        if ncurve == 1
            yprodval = yarray(:,1,i1,j1).*yarray(:,1,i2,j2);
        else
            yprodval = mean(squeeze(yarray(:,:,i1,j1)).* ...
                            squeeze(yarray(:,:,i2,j2)),2);
        end
        yprod(:,i1,j1,i2,j2) = yprodval;
        yprod(:,i2,j2,i1,j1) = yprodval;
    end
end

%  set up array FPROD to hold mean of functions and their derivatives
%    these are used in computing forcing function coefficients

fprod = zeros(n, nvar, nordp1);
for m1=1:nvar*nordp1
    i1 = mmat(m1,2);
    j1 = mmat(m1,1);
    if ncurve == 1
        fprodval = yarray(:,1,i1,j1);
    else
        fprodval = mean(squeeze(yarray(:,:,i1,j1)),2);
    end
    fprod(:,i1,j1) = fprodval;
end

clear yarray

indexf  = 1:fnbasis;
indexw  = 1:wnbasis;
mmat    = m2ij(difeorder,nvar);
onesn   = ones(n,1);

if any(flambda > 0)
    Hmatf = eval_penalty(fbasis);
end
if any(any(wlambda > 0))
    Hmatw = eval_penalty(wbasis);
end

%  set up array to hold coefficients for basis expansions

alpha = getcoef(ffd0);
beta  = getcoef(wfd0);

onen  = ones(n,1);
onef  = ones(1,fnbasis);
onew  = ones(1,wnbasis);
ind1n = [1,n];

%  --------------  beginning of loop through variables  -------------------

fpdacoef = zeros(fnbasis,nvar);
wpdacoef = zeros(wnbasis,nvar*difeorder,nvar);

mlimit = nvar*difeorder;
for ivar=1:nvar
    %  get number of coefficients to be estimated
    neqns  = festimate(ivar)*fnbasis + sum(westimate(:,ivar))*wnbasis;
    cmat   = zeros(neqns, neqns);
    dmat   = zeros(neqns, 1);
    indexj = indexf;
    Cprod  = trapzmat(fbasismat,fbasismat,delta,onesn);
    if festimate(ivar)
        %  first equations: forcing function entries
        %  right side entry for forcing function
        dmat(indexj) = ...
            trapzmat(fbasismat,onesn,delta,fprod(:,ivar,nordp1));
        %  upper right corner
        cmat(indexj,indexj) = Cprod;
        %  remaining entries: loop through variable-derivative pairs
        mkbasis = fnbasis;
        for m2=1:mlimit;
            i2 = mmat(m2,2);
            j2 = mmat(m2,1);
            Cprod = trapzmat(fbasismat,wbasismat,delta,fprod(:,ivar,j2));
            if westimate(m2,ivar)
                indexk  = indexw  + mkbasis;
                mkbasis = mkbasis + wnbasis;
                cmat(indexj,indexk) = Cprod;
            else
                dmat(indexj) = dmat(indexj) + Cprod*beta(:,m2,ivar);
            end
        end
        %  add penalty if required
        if flambda(ivar) > 0
            cmat(indexj,indexj) = cmat(indexj,indexj) - ...
                        flambda(ivar).*Hmatf;
            if any(alpha(:,ivar))
                dmat(indexj) = dmat(indexj) + ...
                        flambda(ivar).*inprod(fbasisfd,ffd0(ivar));
            end
        end
        mjbasis = fnbasis;
    else
        mjbasis = 0;
    end
    %  loop through variable-derivative pairs
    for m1=1:mlimit
        i1 = mmat(m1,2);
        j1 = mmat(m1,1);
        if westimate(m1,ivar)
            indexj  = indexw  + mjbasis;
            mjbasis = mjbasis + wnbasis;
            %  right side entry for variable-derivative pair
            dmat(indexj) = ...
                trapzmat(wbasismat,onesn,  delta,yprod(:,i1,j1,ivar,nordp1));
            %  first entry: forcing function entry
            Cprod = trapzmat(wbasismat,fbasismat,delta,fprod(:,i1,j1));
            if festimate(ivar)
                %  forcing function estimated
                indexk = indexf;
                cmat(indexj,indexk) = Cprod;
                mkbasis = fnbasis;
            else
                % forcing function fixed
                dmat(indexj) = dmat(indexj) + Cprod*alpha(:,ivar);
                mkbasis = 0;
            end
            %  remaining entries: loop through variable-derivative pairs
            for m2=1:mlimit;
                i2 = mmat(m2,2);
                j2 = mmat(m2,1);
                Cprod = ...
                    trapzmat(wbasismat,wbasismat,delta,yprod(:,i1,j1,i2,j2));
                if westimate(m2,ivar)
                    indexk  = indexw  + mkbasis;
                    mkbasis = mkbasis + wnbasis;
                    cmat(indexj,indexk) = Cprod;
                else
                    dmat(indexj) = dmat(indexj) + Cprod*beta(:,m2,ivar);
                end
            end
            %  add penalty if required
            if wlambda(m1,ivar) > 0
                cmat(indexj,indexj) = cmat(indexj,indexj) - ...
                        wlambda(m1,ivar).*Hmatw;
                if any(getcoef(wfd0(1,ivar)))
                    dmat(indexj) = dmat(indexj) + ...
                        wlambda(m1,ivar).*inprod(wbasisfd,wfd0(m1,ivar));
                end
            end
        end
    end
      
    dvec = -cmat\dmat;
  
    if festimate(ivar)
        fpdacoef(:,ivar) = dvec(indexf);
        mjbasis = fnbasis;
    else
        mjbasis = 0;
    end
    for m1=1:mlimit;
        i1 = mmat(m1,2);
        j1 = mmat(m1,1);
        if westimate(m1,ivar) == 1;
            indexj  = indexw  + mjbasis;
            mjbasis = mjbasis + wnbasis;
            wpdacoef(:,m1,ivar) = dvec(indexj);
        end
    end
end

%  --------------  end of loop through variables  -------------------

%  set up the functional data object FFD
ffdnames = getnames(fdobj);
ffdnames{2} = 'Forcing function';
ffdnames{3} = 'Forcing function value';
ffd = fd(fpdacoef, fbasis, ffdnames);

%  set up the functional data object WFD
wfdnames = getnames(fdobj);
wfdnames{2} = 'Weight functions';
wfdnames{3} = 'Weight value';
wfd = fd(wpdacoef, wbasis, wfdnames);

%  set up residual functional data object RESFD
wzero   = zeros(wnbasis,1);
onesnc  = ones(1,ncurve);
ind     = 1:difeorder;
resmat  = zeros(n, ncurve, nvar);
fdarray = zeros(n, ncurve, nvar, nordp1);
for j = 0:difeorder
    fdarray(:,:,:,j+1)  = eval_fd(fdobj, x, j);
end
wmat = eval_fd(wfd, x);
fmat = eval_fd(ffd, x);
for ivar=1:nvar
    m = 0;
    for iord=1:difeorder
        for jvar=1:nvar
           m = m + 1;
            resmat(:,:,ivar) = resmat(:,:,ivar) + ...
                (wmat(:,m,ivar)*onesnc).*fdarray(:,:,jvar,iord);
        end
    end
    resmat(:,:,ivar) = resmat(:,:,ivar) + ...
        fmat(:,ivar)*onesnc + fdarray(:,:,ivar,nordp1);
end
resfd = data2fd(resmat, x, basis);

resfdnames = getnames(resfd);
resfdnames{2} = 'Residual function';
resfdnames{3} = 'Residual function value';
resfd = putnames(resfd, resfdnames);



