function wfd = ptwiseLM(xfd, yfd, wbasis, estimate, ...                        lambda, wfd0, constant, n)
%PTWISELM computes the basis function expansion of the estimate of the 
%    coefficient functions for a pointwise linear model predicting the 
%    function(S) in functional data object YFD from the p functions in
%    functional data object XFD.
%  The coefficient functions are expanded in terms of the basis functions 
%    specified in WBASIS.

%  Arguments:
%  XFD       ...  Functional data object for P independent variable
%                   functions.  These must be univariate.
%  YFD       ...  Functional data object for dependent variable
%                   functions
%  WBASIS    ...  Basis object for regression functions.
%  ESTIMATE  ...  Logical array of length P, if a value is T, the
%                   corresponding coefficient function is estimated, otherwise
%                   the target value is used.
%  LAMBDA    ...  Penalty parameter for penalizing the departure of the
%                   estimated weight functions from those defined in WFD0
%  WFD0      ...  A specification of a functional data object that is used for
%                   those weight functions not estimated, or as target functions
%                   toward which the estimated weight functions are smoothed. WFD0
%                   can either be a vector of NCOEF constants, or a functional
%                   data object with the same structure as WFN that is returned
%                   by this function.
%  CONSTANT  ...  If nonzero, a constant function is added to the fit.
%  N         ...  Number of sampling points for numerical integration

%  Returns:
%  WFD       ...  Estimated weight functional data object.  It has P + 1
%                   replicates if CONSTANT is T, and P otherwise

%  Last modified 14 January 2003

if ~(isa(xfd, 'fd'))
    error('Argument XFD not a functional data object.');
end
if ~(isa(yfd, 'fd'))
    error('Argument YFD not a functional data object.');
end

if nargin < 7, constant = 1; end
coefx  = getcoef(xfd);
coefdx = size(coefx);
ndimx  = length(coefdx);
ncurve = coefdx(2);

coefy  = getcoef(yfd);
coefdy = size(coefy);
ndimy  = length(coefdy);
if ndimy > 2 
    nvar = coefdy(3); 
else 
    nvar = 1;
end

if coefdy(2) ~= ncurve
    error('Number of replications for XFD and YFD are not the same.');
end

xbasis  = getbasis(xfd);
nbasisx = getnbasis(xbasis);
if nargin < 8, n = max([10*nbasisx; 101]);  end

ybasis  = getbasis(yfd);
nbasisy = getnbasis(ybasis);
if ndimx == 2
    coefx = reshape(coefx,[nbasisx,ncurve,1]);end
if ndimx > 2, p  = coefdx(3); else p = 1; end
if nargin < 4
    if constant
        estimate = ones(p+1,1);
    else
        estimate = ones(p,1);
    end
end

if nargin < 3, wbasis = xbasis;  end
typew   = getbasistype(wbasis);
nbasisw = getnbasis(wbasis);
rangew  = getbasisrange(wbasis);

if any(rangew ~= getbasisrange(xbasis))  
    error('Weight function range not equal to range in XFD');
end

if typew == 'bspline'
    nbreaksw = length(getbasispar(wbasis));
    norderw  = nbasisw - nbreaksw;
end

delta = (rangew(2)-rangew(1))/(n-1);
tfine = (rangew(1):delta:rangew(2))';

yarray = eval_fd(yfd, tfine);
if constant
    ncoef  = p + 1 - length(find(estimate==0));
    xarray = ones(n,ncurve,p+1);
    xarray(:,:,2:(p+1)) = eval_fd(tfine, xfd);
else
    ncoef  = p - length(find(estimate==0));
    xarray = eval_fd(xfd, tfine);
end

if nargin < 6, wfd0   = zeros(ncoef,1); end
if nargin < 5, lambda = zeros(ncoef,1); end

basismat = getbasismatrix(tfine, wbasis);

if ncurve == 1
    DV = -delta.*yarray;
    IV = zeros(n,ncoef*nbasisw);
    mi = 0;
    for i = 1:ncoef
        if estimate(i)
            mi = mi + 1;
            index = (1 + (mi-1)*nbasisw):(mi*nbasisw);
            IV(:,index) = delta.*(xarray(:,i)*ones(1,nbasisw)).*basismat;
        end
    end
    dvec = IV\DV;
else
    mi  = 0;
    mij = 0;
    Swgt = zeros(n,ncoef);
    Rwgt = zeros(n,ncoef*(ncoef+1)/2);
    for i = 1:ncoef
        if estimate(i)
            mi = mi + 1;
            index = (1 + (mi-1)*nbasisw):(mi*nbasisw);
            Swgt(:,mi) = mean((squeeze(xarray(:,:,i)).*yarray)')';
            mj = 0;
            for j = 1:i
                if estimate(j)
                    mj  = mj  + 1;
                    mij = mij + 1;
                    temp = squeeze(xarray(:,:,mi)).*squeeze(xarray(:,:,mj));
                    Rwgt(:,mij) = mean(temp')';
                end
            end
        end
    end
    
    [Dmat,Cmat] = SRsetup(ncoef, nbasisw, Swgt, Rwgt, basismat);
    
    if any(lambda > 0)
        if ~isa(wfd0, 'fd') & isa(wfd0, 'double')
            if length(wfd0) ~= ncoef
                error('WFD0 is a vector of incorrect length');
                wfd0 = fd(reshape(wfd0,1,ncoef), wbasis0);
            end
        else
            error('WFD0 is neither a vector nor a FD object');
        end
        Hmat = eval_penalty(wbasis0);
        for i = 1:ncoef
            index = (1 + (i-1)*nbasisw):(i*nbasisw);
            if lambda(i) > 0
                Cmat(index,index) = Cmat(index,index) - lambda(i)*Hmat;
                Dmat(index,1)     = Dmat(index,1) + ...
                    lambda(i).*inprod(wbasis0,wfn0(i));
            end
        end
    end
    dvec   = Cmat\Dmat;
end

dmat = zeros(nbasisw,ncoef);
mi = 0;
for i = 1:ncoef
    if estimate(i)
        mi = mi + 1;
        index = (1 + (mi-1)*nbasisw):(mi*nbasisw);
        dmat(:,i) = dvec(index);
    end
end

wfd = fd(dmat, wbasis);

%  ---------------------------------------------------------------------

function [Dmat, Cmat] = SRsetup(ncoef, nbasis, Swgt, Rwgt, basismat)
%  sets up coefficient matrices for basis expansion of weight functions
[nx,m1] = size(basismat);
nrow = ncoef*nbasis;
Dmat = zeros(nrow, 1);
Cmat = zeros(nrow, nrow);
m = 0;
onem1 = ones(1,m1);
for i = 1:ncoef
    indexi = (1:nbasis) + (i-1)*nbasis;
    temp     = basismat .* (Swgt(:,i)*onem1);
    temp([1,nrow],:) = temp([1,nrow],:)./2;
    Dmat(indexi) = sum(temp)';
    for j = 1:i
        m = m + 1;
        indexj = (1:nbasis) + (j-1)*nbasis;
        temp = basismat .* (Rwgt(:,m)*onem1);
        temp([1,nrow],:) = temp([1,nrow],:)./2;
        Cmat(indexi,indexj) = temp'*basismat;
        if i ~= j, Cmat(indexj,indexi) = Cmat(indexi,indexj); end
    end
end


