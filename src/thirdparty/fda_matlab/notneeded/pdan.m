function wfd = pdan(fdobj, difeorder, wbasis, estimate, lambda, wfd0, n)
%PDAN computes the basis function expansion of the
%    estimate of the coefficient functions w_j(t) and forcing function f(t)
%    in the nonhomogeneous linear differential operator
%
%    Lx(t) = f(t) +
%       w_0(t)x(t) + w_1(t)Dx(t) + ... + w_{m-1}D^{m-1}x(t) + D^m x(t)
%
%    of order m = NORDER that minimizes in a least squares sense the residual
%    functions Lx(t).  The functions x(t) are in functional data object FDOBJ.
%
%  Unlike PDA, the functions x can be vector valued. That is, the differential 
%    equation can be a system of K equations rather than a single equation.   
%  The coefficient function w_j(t) is matrix valued, with a column
%    for each scalar function in the system.  
%  The forcing function f(t) is vector valued.
%  The forcing functions and weight functions are returned as a multivariate
%    FD object, with the forcing function and coefficients stacked associated
%    with a particular equation or variable.  If there are K equations, then
%    the FD object has K variables, and the number of curves is 1 + m*K.  
%    The first curve is the forcing function for that equation,
%    the next K the weight functions w_0(t), the next K the weight functions 
%    w_1(t), and so on.  
%
%  For example, suppose we have a system of two second order equations, so that
%    Lx(t) = f(t) + w_0(t)x(t) + w_1(t)Dx(t) + D^2 x(t)
%  Then both w_0 and w_1 are 2 by 2 matrix-valued.  The FD object WFD that
%    is returned will have 2 variables and 4 replications.  For each variable
%    the replications are:  first the two weight functions w_0 and then the two
%    weight functions w_1 associated with that variable or equation.
%  The coefficient functions are expanded in terms of the
%    basis functions specified in wbasis.
%  Arguments:
%  FDOBJ     ...  functional data object, or a matrix of function values
%  DIFEORDER ...  order of the linear differential operator, that is, the order
%                 of the highest derivative.
%  WBASIS    ...  basis object for weight functions
%  ESTIMATE  ...  A matrix of 0's and 1's with the number of rows being
%                 1+m*K (see discription of WFD above) and number cols K.
%                 If a value is 1, the corresponding coefficient function is
%                 estimated, otherwise the target value is used.
%  LAMBDA    ...  matrix with the same dimensions as ESTIMATE containing
%                 penalty parameters for penalizing the departure of the
%                 estimated weight functions from those defined = WFD0
%  WFD0      ...  A specification of a functional data object that is used for
%                 those weight functions not estimated, or as target functions
%                 toward which the estimated weight functions are smoothed. WFD0
%                 can either be a vector of DIFEORDER constants, or a functional
%                 data object with the same structure as WFN that is returned
%                 by this function.
%  N         ...  number of sampling points for numerical integration

%  Returns:
%  WFN       ...  estimated weight functional data object.  It has 1+DIFEORDER*K
%                 replicates, and K variables.

%  last modified 14 January 2003

if nargin < 3
    error('There are less than three arguments.');
end

if ~isa_fd(fdobj)
    error('FDOBJ is not a functional data object.');
end

norder = floor(difeorder);
if norder < 1
    error('NORDER is less than one.');
end
nordp1 = norder + 1;
nordp2 = norder + 2;

if ~isa_basis(wbasis)
    error('wbasis is not a basis object.');
end

basisfd  = getbasis(fdobj);
nbasis   = getnbasis(basisfd);
rangeval = getbasisrange(basisfd);

%  get the dimensions of the data in FD

coef   = getcoef(fdobj);
coefd  = size(coef);
ndim   = length(coefd);
ncurve = coefd(2);
if ndim == 3
    nvar = coefd(3);
else
    nvar = 1;
end
if nvar == 1
    error(['PDAN is designed for a system of more than one equation.', ...
           'Use PDA instead.']);
end
nfun = 1 + nvar*norder;

%  set default arguments

if nargin < 7, n        = 5*nbasis;           end
if nargin < 6, estimate = ones (nfun,nvar);  end
if nargin < 5, wfd0     = zeros(nfun,nvar);  end
if nargin < 4, lambda   = zeros(nfun,nvar);  end

% check array ESTIMATE

if size(estimate,1) ~= nfun
    error('ESTIMATE has wrong number of rows.');
end
if size(estimate,2) ~= nvar
    error('ESTIMATE has wrong number of columns.');
end
if ~all(estimate == 0 | estimate == 1)
    error('ESTIMATE does not contain only zeros or ones.');
end

% check array LAMBDA

if size(lambda,1) ~= nfun
    error('LAMBDA has wrong number of rows.');
end
if size(lambda,2) ~= nvar
    error('LAMBDA has wrong number of columns.');
end
if any(lambda < 0)
    error('LAMBDA contains negative numbers.');
end

%  check WFD0 and make into a FD object if numeric

if ~strcmp(class(wfd0), 'fd') & ~isnumeric(wfd0)
    error('WFD0 is neither a vector nor a FD object');
end

if ~strcmp(class(wfd0), 'fd') 
    if length(wfd0) ~= nordp1 
        error('WFD0 is a vector of incorrect length');
    end
    wbasis0 = create_constant_basis(rangeval);
    wfd0    = fd(wfd0, wbasis0);
end

%  check and get the characteristics of the basis to be used

typew   = getbasistype(wbasis);
nbasisw = getnbasis(wbasis);
rangew  = getbasisrange(wbasis);

if any(rangew ~= rangeval)
    error('Weight function range not equal to FDOBJ range');
end

%  why do we need this?
if strcmp(typew, 'bspline')
    params   = getbasispar(wbasis);
    nbreaksw = length(params);
    norderw  = nbasisw - nbreaksw;
end

%  set up sampling values to be used in numerical integration
%    and set up matrix of basis values

delta    = (rangew(2)-rangew(1))/(n-1);
x        = rangew(1):delta:rangew(2);
basismat = getbasismatrix(x, wbasis);

%  set up array to hold values of functions and their derivatives
%   for each variable

yarray = zeros([n,ncurve,nvar,nordp1]);
for icurve=1:ncurve
    for ivar=1:nvar
        for j=1:nordp1
            yarray(:,icurve,ivar,j) = eval_fd(fdobj(:,ivar), x, j-1);
        end
    end
end

%  set up array YPROD to hold mean of products of values in YARRAY

mmat  = m2ij(nvar,nordp1);
yprod = zeros([n,nvar,nordp1,nvar,nordp1]);
for m1=1:nvar*nordp1
    i1 = mmat(m1,1);
    j1 = mmat(m1,2);
    for m2=1:m1;
        i2 = mmat(m2,1);
        j2 = mmat(m2,2);
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
    i1 = mmat(m1,1);
    j1 = mmat(m1,2);
    if ncurve == 1
        fprodval = yarray(:,1,i1,j1);
    else
        fprodval = mean(squeeze(yarray(:,:,i1,j1)),2);
    end
    fprod(:,i1,j1) = fprodval;
end
      
indexw  = 1:nbasisw;
mmat    = m2ij(nvar,norder);
onesn   = ones(n,1);
coefmat = zeros(nbasisw,nfun,nvar);

if any(any(lambda > 0))
    Hmat = eval_penalty(wbasis);
end

%  loop through equations

for ivar=1:nvar
    %  get number of coefficients to be estimated
    ncoef  = sum(estimate(:,ivar));
    cmat   = zeros(ncoef*nbasisw, ncoef*nbasisw);
    dmat   = zeros(ncoef*nbasisw, 1);
    index1 = indexw;
    if estimate(1,ivar) == 1
        %  first row: forcing function entries
        %  right side entry for forcing function
        dmat(index1) = ...
            trapzmat(basismat,onesn,   delta,fprod(:,ivar,nordp1));
        %  upper right corner
        cmat(index1,index1) = ...
            trapzmat(basismat,basismat,delta,onesn);
        if lambda(1,ivar) > 0
            cmat(index1,index1) = cmat(index1,index1) - lambda(1,ivar).*Hmat;
            if any(getcoef(wfd0(1,ivar)))
                dmat(index1) = dmat(index1) + ...
                            lambda(1,ivar).*inprod(basisfd,wfd0(1,ivar));
            end
        end
        index1 = index1 + nbasisw;
    end
    %  loop through variable-derivative pairs
    for m1=1:nvar*norder
        i1 = mmat(m1,1);
        j1 = mmat(m1,2);
        if estimate(m1+1,ivar) == 1
            index2 = indexw;
            %  right side entry for variable-derivative pair
            dmat(index1) = ...
                trapzmat(basismat,onesn,  delta,yprod(:,i1,j1,ivar,nordp1));
            %  first entry: forcing function entry
            cmat(index1,index2) = ...
                trapzmat(basismat,basismat,delta,fprod(:,i1,j1));
            if lambda(m1+1,ivar) > 0
                cmat(index1,index2) = cmat(index1,index2) - ...
                        lambda(m1+1,ivar).*Hmat;
                if any(getcoef(wfd0(1,ivar)))
                    dmat(index1) = dmat(index1) + ...
                        lambda(m1+1,ivar).*inprod(basisfd,wfd0(m1+1,ivar));
                end
            end
            cmat(index2,index1) = cmat(index1,index2);
            index2 = index2 + nbasisw;
            %  remaining entries: loop through variable-derivative pairs
            for m2=1:m1;
                i2 = mmat(m2,1);
                j2 = mmat(m2,2);
                if estimate(m2+1,ivar) == 1
                    cmat(index1,index2) = ...
                        trapzmat(basismat,basismat,delta,yprod(:,i1,j1,i2,j2));
                    if lambda(m2+1,ivar) > 0
                        cmat(index1,index2) = cmat(index1,index2) - ...
                                lambda(m2+1,ivar).*Hmat;
                        if any(getcoef(wfd0(1,ivar)))
                            dmat(index1) = dmat(index1) + ...
                                lambda(m2+1,ivar).*inprod(basisfd,wfd0(m2+1,ivar));
                        end
                    end
                    cmat(index2,index1) = cmat(index1,index2);
                    index2 = index2 + nbasisw;
                end
            end
            index1 = index1 + nbasisw;
        end
    end
    
  
    dvec = -symsolve(cmat, dmat);
  
    index = indexw;
    if estimate(1,ivar) == 1
        coefmat(:,1,ivar) = dvec(index);
        index = index + nbasisw;
    end
    for m1=1:nvar*norder;
        i1 = mmat(m1,1);
        j1 = mmat(m1,2);
        if estimate(m1+1,ivar) == 1;
            coefmat(:,m1+1,ivar) = dvec(index);
            index = index + nbasisw;
        end
    end
end

wfdnames{1} = 'Time';
wfdnames{2} = 'Weight functions';
wfdnames{3} = 'Weight value';
wfd = fd(coefmat, wbasis, wfdnames);

