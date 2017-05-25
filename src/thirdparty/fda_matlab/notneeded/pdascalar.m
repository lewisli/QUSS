function [ffd, wfd, resfd] = pdascalar(fdobj, norder, ...
                          fbasis, flambda, ffd0, festimate, ...
                          wbasis, wlambda, wfd0, westimate, n)
%PDASCALAR computes the basis function expansion of the
%  estimate of the coefficient functions w_j(t) and forcing function f(t)
%  in the nonhomogeneous linear differential operator
%
%    Lx(t) = f(t) +
%            w_0(t)x(t) + w_1(t)Dx(t) + ... + w_{m-1}D^{m-1}x(t) + 
%            D^m x(t)
%
%  of order m = norder that minimizes in a least squares sense the residual
%  functions Lx(t).  
%  The functions x(t) are in functional data object FDOBJ.
%  The coefficient functions are expanded in terms of the
%  basis functions specified in FBASIS and WBASIS.

%  Arguments:
%  FDOBJ     ...  functional data object
%  NORDER    ...  order of the linear differential operator, that is, 
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
%  FESTIMATE ...  logical constant, if a value is T, the forcing function  
%                 is estimated, otherwise the target value is used.  
%  WBASIS    ...  basis object for the weight functions.
%  WLAMBDA   ...  a numerical array of length NORDER containing
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
%  WESTIMATE ...  logical array of length NORDER, if a value is T, the
%                 corresponding coefficient function is estimated, 
%                 otherwise the target value is used.  
%  N         ...  number of sampling points for numerical integration

%  Returns:
%  FFD       ...  FD object for estimated forcing function.  
%  WFD       ...  FD object for estimated weight  functions.  
%                 It has NORDER replicates.
%  RESFD     ...  FD object for residual functions

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

norder = floor(norder);
if norder < 1, error('NORDER is less than one.');  end
nordp1 = norder + 1;

%  get information about data basis

basis    = getbasis(fdobj);
nbasis   = getnbasis(basis);
rangeval = getbasisrange(basis);

%  set default arguments

if nargin < 11, n = 5*nbasis; end

if nargin < 10, westimate =  ones(norder,1); end
if nargin <  9, wfd0      = 0; end
if nargin <  8, wlambda   = zeros(norder,1); end
if nargin <  7, wbasis    = create_constant_basis(rangeval); end
if nargin <  6, festimate = 0; end
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

if size(flambda,1) ~= 1 | size(flambda,2) ~= 1
    error('FLAMBDA is not a scalar constant.');
end
if ~(flambda >= 0)
    error('FLAMBDA does not contain a nonegative value.');
end

%  check array FESTIMATE

if size(festimate,1) ~= 1 | size(festimate,2) ~= 1
    error('FESTIMATE is not a scalar constant.');
end
if ~(festimate == 0 | festimate == 1)
    error('FESTIMATE does not contain either zero or one.');
end

%  check FFD0 and make into a FD object if numeric

if ~strcmp(class(ffd0), 'fd') & ~isnumeric(ffd0)
    error('FFD0 is neither a vector nor a FD object');
end
if ~strcmp(class(ffd0), 'fd') 
    if length(ffd0) == 1
        ffd0 = ffd0.*ones(fnbasis,1);
    end
    if size(ffd0,2) ~= 1 | size(ffd0,1) ~= fnbasis
        error('FFD0 has incorrect dimensions.');
    end
    ffd0 = fd(ffd0, fbasis);
end

%  check array WLAMBDA

if size(wlambda,1) ~= norder | size(wlambda,2) ~= 1
    error('WLAMBDA is not a vector of length NORDER.');
end
if ~all(wlambda >= 0)
    error('WLAMBDA does not contain only nonegative values.');
end

%  check array WESTIMATE

if size(westimate,1) ~= norder | size(westimate,2) ~= 1
    error('WESTIMATE is not a vector of length NORDER.');
end
if ~all(westimate == 0 | westimate == 1)
    error('WESTIMATE does not contain only zeros or ones.');
end

%  check WFD0 and make into a FD object if numeric

if ~strcmp(class(wfd0), 'fd') & ~isnumeric(wfd0)
    error('WFD0 is neither a vector nor a FD object');
end
if ~strcmp(class(wfd0), 'fd') 
    if length(wfd0) == 1
        wfd0 = wfd0.*ones(wnbasis,norder);
    end
    wfd0size = size(wfd0);
    if length(wfd0size) ~= 2
        error(['A numerical WFD0 is neither a constant ',...
               'nor a two-dimensional array.']);
    end
    if wfd0size(1) ~= wnbasis | wfd0size(2) ~= norder  
        error('WFD0 has incorrect dimensions.');
    end
    wfd0 = fd(wfd0, wbasis);
end

nfcoef = festimate;
nwcoef = sum(westimate);
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
if flambda > 0
    Hmatf = eval_penalty(fbasis, 1);
end
if any(wlambda > 0)
    Hmatw = eval_penalty(wbasis, 1);
end

%  set up array to hold values of functions and their derivatives

yarray = zeros([n,ncurve,nordp1]);

%  set up array to hold coefficients for basis expansions

alpha = getcoef(ffd0);
beta  = getcoef(wfd0);

indexf  = 1:fnbasis;
indexw  = 1:wnbasis;

neqns = nfcoef*fnbasis + nwcoef*wnbasis;
onen  = ones(n,1);
onef  = ones(1,fnbasis);
onew  = ones(1,wnbasis);
ind1n = [1,n];
Swgt  = zeros(n,nordp1);
Rwgt  = zeros(n,nordp1,nordp1);

%  --------------  beginning of loop through variables  -------------------

for ivar=1:nvar
    
    %  fill yarray with values of functions and their derivatives
    
    for j = 0:norder 
        if nvar == 1
            yarray(:,:,j+1) = eval_fd(fdobj, x, j);
        else
            yarray(:,:,j+1) = squeeze(eval_fd(fdobj(:,ivar), x, j));
        end
    end
    
    %  now compute all the mean derivative products
    %  SWGT has forcing and weight functions paired with mth derivative
    %  RWGT has all pairs from forcing + weight functions
    
    %  forcing function
    %  1 is paired with mth derivative
    Swgt(:,1) = mean((yarray(:,:,nordp1)),2);
    %  forcing function paired with itself:  multipying function is 1
    Rwgt(:,1,1) = onen;
    %  each weight function
    for j = 1:norder
        %  jth derivative paired with mth derivative
        Swgt(:,j+1) = mean((yarray(:,:,j).*yarray(:,:,nordp1)),2);
        %  jth derivative paired with forcing function
        Rwgt(:,1,j+1) = mean((yarray(:,:,j)),2);
        Rwgt(:,j+1,1) = Rwgt(:,1,j);
        %  jth function paired with each derivative up to j
        for k = 1:j
            Rwgt(:,j+1,k+1) = mean((yarray(:,:,j).*yarray(:,:,k)),2);
            Rwgt(:,k+1,j+1) = Rwgt(:,j+1,k+1);
        end
    end

    %  set up left and right sides of linear equation
    %  if function is to be estimated, set up a block of lines  
    %  in the linear equation, one line per basis function, 
    %  else add that function to the left side
    
    Dmat = zeros(neqns, 1);
    Cmat = zeros(neqns, neqns);
    
    %  forcing function equation
    Cprod  = trapzmat(fbasismat,fbasismat,delta,Rwgt(:,1,1));
    indexj = indexf;
    if festimate
        %  ffd to be estimated
        Dprod = trapzmat(fbasismat, onen, delta, Swgt(:,1));
        Dmat(indexj)        = Dprod;
        Cmat(indexj,indexj) = Cprod;
        mkbasis = fnbasis;
        for k = 1:norder
            Cprod  = trapzmat(fbasismat,wbasismat,delta,Rwgt(:,1,k+1));
            if westimate(k)
                indexk = indexw + mkbasis;
                mkbasis = mkbasis + wnbasis;
                Cmat(indexj,indexk) = Cprod;
            else
                Dmat(indexj) = Dmat(indexj) + Cprod*beta(:,j,ivar);
            end
        end
       %  add penalty if required
        if flambda > 0
            Cmat(indexj,indexj) = Cmat(indexj,indexj) - ...
                        flambda.*Hmatf;
            if any(alpha ~= 0)
                Dmat(indexj) = Dmat(indexj) + ...
                        flambda.*inprod(fbasis,ffd0);
            end
        end
        mjbasis = fnbasis;
        mkbasis = fnbasis;
    else
        mjbasis = 0;
        mkbasis = 0;
    end
    %  equation for jth derivative
    for j = 1:norder
        if westimate(j)
            %  beta_j to be estimated, set up equation
            indexj = indexw + mjbasis;
            mjbasis = mjbasis + wnbasis;
            %  right hand values
            Dprod = trapzmat(wbasismat,onen,delta,Swgt(:,j+1));
            Dmat(indexj) = Dprod;
            %  jth derivative paired with forcing function
            Cprod  = trapzmat(wbasismat,fbasismat,delta,Rwgt(:,j+1,1));
            if festimate
                %  forcing function estimated
                indexk = indexf;
                Cmat(indexj,indexk) = Cprod;
                mkbasis = fnbasis;
            else
                %  forcing function fixed
                Dmat(indexj) = Dmat(indexj) + Cprod*alpha(:,ivar);
                mkbasis = 0;
            end
            %  jth derivative paired with kth derivative
            for k = 1:norder
                Cprod  = trapzmat(wbasismat,wbasismat,delta,Rwgt(:,j+1,k+1));
                if westimate(k)
                    indexk = indexw + mkbasis;
                    mkbasis = mkbasis + wnbasis;
                    Cmat(indexj,indexk) = Cprod;
                else
                    Dmat(indexj) = Dmat(indexj) + Cprod*beta(:,j,ivar);
                end
            end
            %  add penalty of required
            if wlambda(j) > 0
                Cmat(indexj,indexj) = Cmat(indexj,indexj) - ...
                        wlambda(j).*Hmatw;
                if nvar == 1
                    wfdj = wfd0(j);
                else
                    wfdj = wfd0(j,ivar);
                end
                if any(beta(:,j) ~= 0)
                    Dmat(indexj) = Dmat(indexj) + ...
                            wlambda(j).*inprod(wbasis,wfdj);
                end
            end
        end
    end
    
    %  solve the equation 
 
    dvec = -Cmat\Dmat;

    %  set up the coefficient matrices

    fDmat = zeros(fnbasis,1);
    wDmat = zeros(wnbasis,norder);
    if festimate
        fDmat = dvec(indexf);
        mjbasis = fnbasis;
    else
        fDmat = getcoef(ffd0);
        mjbasis = 0;
    end
    for i = 1:norder
        if westimate(i)
            indexj = indexw + mjbasis;
            mjbasis = mjbasis + wnbasis;
            wDmat(:,i) = dvec(indexj);
        else
            if nvar == 1
                wDmat(:,i) = getcoef(wfd0(i));
            else
                wDmat(:,i) = getcoef(wfd0(i,ivar));
            end
        end
    end
    if nvar == 1 
        fpdacoef = fDmat;
        wpdacoef = wDmat;
    else 
        fpdacoef(:,  ivar) = fDmat;
        wpdacoef(:,:,ivar) = wDmat;
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
onesn = ones(1,ncurve);
if nvar == 1
    wcoef  = [zeros(wnbasis,1), wpdacoef];
    Lfd    = fd(wcoef,wbasis);
    resmat = eval_fd(fdobj, x, Lfd) + eval_fd(ffd, x)*onesn;
    resfd  = data2fd(resmat, x, basis);
else
    resmat = zeros(n, ncurve, nvar);
    for ivar=1:nvar
        wcoef = [zeros(wnbasis,1), squeeze(wpdacoef(:,:,ivar))];
        Lfd   = fd(wcoef,wbasis);
        resmat(:,:,ivar) = eval_fd(fdobj(:,ivar), x, Lfd)  + ...
                           eval_fd(ffd(ivar), x)*onesn;
    end
    resfd  = data2fd(resmat, x, basis);
end
resfdnames = getnames(resfd);
resfdnames{2} = 'Residual function';
resfdnames{3} = 'Residual function value';
resfd = putnames(resfd, resfdnames);



