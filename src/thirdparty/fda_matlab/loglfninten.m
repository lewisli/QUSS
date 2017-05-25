function [logl, Dlogl] = loglfninten(x, basisfd, cvec)
nobs    = length(x);
Cval    = normalize_phi(x, basisfd, cvec);
phimat  = full(eval_basis(x, basisfd));
logl    = sum(phimat * cvec) - Cval;
EDW     = expect_phi(x, basisfd, cvec);
Dlogl   = sum(phimat)' - EDW;

function Cval = normalize_phi(eventtimes, basisfd, cvec, JMAX, EPS)

%  Computes integral from RNG(1) to EVENTTIMES(nobs) of
%      \exp [phi'(x) * cvec]
%  by numerical integration using Romberg integration

%  Arguments:
%  EVENTTIMES ...  Vector of event times
%  BASISFD    ...  Basis function object with basis functions phi.
%  CVEC       ... coefficient vector defining density, of length NBASIS
%  JMAX       ... maximum number of allowable iterations
%  EPS        ... convergence criterion for relative error

%  Return:
%  The integral of the function.

%  check arguments, and convert basis objects to functional data objects

if ~strcmp(class(basisfd), 'basis')
    error('First argument must be a basis function object.');
end

nbasis = getnbasis(basisfd);
oneb   = ones(1,nbasis);
rng    = getbasisrange(basisfd);
nobs   = length(eventtimes);
rng(2) = eventtimes(nobs);

%  set default arguments

if nargin < 5, EPS  = 1e-7; end
if nargin < 4, JMAX = 15;   end

%  set up first iteration

width = rng(2) - rng(1);
JMAXP = JMAX + 1;
h = ones(JMAXP,1);
h(2) = 0.25;
%  matrix SMAT contains the history of discrete approximations to the integral
smat = zeros(JMAXP,1);
%  the first iteration uses just the endpoints
x  = rng';
nx = length(x);
ox = ones(nx,1);
fx = full(eval_basis(x, basisfd));
wx = fx * cvec;
wx(wx < -50) = -50;
px = exp(wx);
smat(1)  = width.*sum(px)./2;
tnm = 0.5;
j   = 1;

%  now iterate to convergence

for j = 2:JMAX
    tnm  = tnm*2;
    del  = width/tnm;
    if j == 2
        x = (rng(1) + rng(2))/2;
    else
        x = rng(1)+del/2 : del : rng(2);
    end
    fx = full(eval_basis(x, basisfd));
    wx = fx * cvec;
    wx(wx < -50) = -50;
    px = exp(wx);
    smat(j) = (smat(j-1) + width.*sum(px)./tnm)./2;
    if j >= 5
        ind = (j-4):j;
        [ss, dss] = polintarray(h(ind),smat(ind),0);
        if ~any(abs(dss) >= EPS.*max(abs(ss)))
            %  successful convergence
            Cval = ss;
            return;
        end
    end
    smat(j+1) = smat(j);
    h(j+1)   = 0.25*h(j);
end
warning(['No convergence after ',num2str(JMAX),' steps in NORMALIZE.PHI'])

%  ---------------------------------------------------------------

function ss = expect_phi(eventtimes, basisfd, cvec, JMAX, EPS)
%  Computes inner products of basis functions with 
%      p(x) = exp [c'phi(x)],
%  where the upper limit of integration is 
%  EVENTTIMES(NOBS),
%  by numerical integration using Romberg integration

%  Arguments:
%  EVENTTIMES ...  Vector of event times
%  BASISFD  ...  A basis function
%           object.  In the latter case, a functional data object
%           is created from a basis function object by using the
%           identity matrix as the coefficient matrix.
%           The functional data objects must be univariate.
%  CVEC ... coefficient vector defining density, of length NBASIS
%  JMAX ... maximum number of allowable iterations
%  EPS  ... convergence criterion for relative error

%  Return:
%  A vector SS of length NBASIS of inner products.

%  check arguments, and convert basis objects to functional data objects

if ~strcmp(class(basisfd),'basis')
    error('First argument must be a basis function object.');
end

nbasis = getnbasis(basisfd);
oneb   = ones(1,nbasis);
rng    = getbasisrange(basisfd);
nobs   = length(eventtimes);
rng(2) = eventtimes(nobs);

%  set default arguments

if nargin < 5, EPS  = 1e-7; end
if nargin < 4, JMAX = 15;   end

%  set up first iteration

width = rng(2) - rng(1);
JMAXP = JMAX + 1;
h = ones(JMAXP,1);
h(2) = 0.25;
%  matrix SMAT contains the history of discrete approximations to the integral
smat = zeros(JMAXP,nbasis);
sumj = zeros(1,nbasis);
%  the first iteration uses just the endpoints
x  = rng';
nx = length(x);
ox = ones(nx);
fx = full(eval_basis(x, basisfd));
wx = fx * cvec;
wx(wx < -50) = -50;
px   = exp(wx);
sumj = fx'*px;
smat(1,:) = width.*sumj'./2;
tnm = 0.5;
j   = 1;

%  now iterate to convergence

for j = 2:JMAX
    tnm  = tnm*2;
    del  = width/tnm;
    if j == 2
        x = (rng(1) + rng(2))/2;
    else
        x = (rng(1)+del/2 : del : rng(2))';
    end
    fx = full(eval_basis(x, basisfd));
    wx = fx * cvec;
    wx(wx < -50) = -50;
    px   = exp(wx);
    sumj = fx' * px;
    smat(j,:) = (smat(j-1,:) + width.*sumj'./tnm)./2;
    if j >= 5
        ind = (j-4):j;
        temp = squeeze(smat(ind,:));
        [ss, dss] = polintarray(h(ind),temp,0);
        if ~any(abs(dss) > EPS*max(abs(ss)))
            %  successful convergence
            return;
        end
    end
    smat(j+1,:) = smat(j,:);
    h(j+1) = 0.25*h(j);
end
warning(['No convergence after ',num2str(JMAX),' steps in EXPECT.PHI'])

%  ---------------------------------------------------------------

function [y,dy] = polintarray(xa, ya, x)
%  YA is an array with up to 4 dimensions
%     with 1st dim the same length same as the vector XA
n     = length(xa);
yadim = size(ya);
nydim = length(yadim);
if yadim(2) == 1, nydim = 1; end
if yadim(1) ~= n, error('First dimension of YA must match XA'); end
difx = xa - x;
absxmxa = abs(difx);
tmp = 1:n;
ns = min(tmp(absxmxa == min(absxmxa)));
cs = ya;
ds = ya;
if nydim == 1, y = ya(ns);  end
if nydim == 2, y = ya(ns,:);  end
if nydim == 3, y = ya(ns,:,:);  end
if nydim == 4, y = ya(ns,:,:,:);  end
ns = ns - 1;
for m = 1:(n-1)
    if nydim == 1
        for i = 1:(n-m)
            ho      = difx(i);
            hp      = difx(i+m);
            w       = (cs(i+1) - ds(i))./(ho - hp);
            ds(i) = hp.*w;
            cs(i) = ho.*w;
        end
        if 2*ns < n-m
            dy = cs(ns+1);
        else
            dy = ds(ns);
            ns = ns - 1;
        end
    end
    if nydim == 2
        for i = 1:(n-m)
            ho      = difx(i);
            hp      = difx(i+m);
            w       = (cs(i+1,:) - ds(i,:))./(ho - hp);
            ds(i,:) = hp.*w;
            cs(i,:) = ho.*w;
        end
        if 2*ns < n-m
            dy = cs(ns+1,:);
        else
            dy = ds(ns,:);
            ns = ns - 1;
        end
    end
    if nydim == 3
        for i = 1:(n-m)
            ho      = difx(i);
            hp      = difx(i+m);
            w       = (cs(i+1,:,:) - ds(i,:,:))./(ho - hp);
            ds(i,:,:) = hp.*w;
            cs(i,:,:) = ho.*w;
        end
        if 2*ns < n-m
            dy = cs(ns+1,:,:);
        else
            dy = ds(ns,:,:);
            ns = ns - 1;
        end
    end
    if nydim == 4
        for i = 1:(n-m)
            ho      = difx(i);
            hp      = difx(i+m);
            w       = (cs(i+1,:,:,:) - ds(i,:,:,:))./(ho - hp);
            ds(i,:,:,:) = hp.*w;
            cs(i,:,:,:) = ho.*w;
        end
        if 2*ns < n-m
            dy = cs(ns+1,:,:,:);
        else
            dy = ds(ns,:,:,:);
            ns = ns - 1;
        end
    end
    y = y + dy;
end
