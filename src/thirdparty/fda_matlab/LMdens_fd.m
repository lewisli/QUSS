function [Wfdobj, beta, Fstr, iternum, iterhist] = ...
    LMdens_fd(y, WfdPar, zmat, beta0, sigma0, ...
              conv, iterlim, active, dbglev)
%LMDENS_FD estimates a regression and the density of the residuals.
%  If betaO is [], or empty, then only the density is estimated. 

%  Arguments are:
%  Y       ... vector of dependent variable values
%  WFDPAR  ... functional parameter object specifying the initial log
%              density, the linear differential operator used to smooth
%              smooth it, and the smoothing parameter.
%  ZMAT    ... matrix of covariates
%  BETA0   ... initial model   coefficients
%  SIGMA0  ... initial standard error
%  CONV    ... convergence criterion
%  ITERLIM ... iteration limit for scoring iterations
%  ACTIVE  ... indices among 1:NBASIS of parameters to optimize 
%              Normally the first index is set inactive.
%  DBGLEV  ... level of output of computation history

%  Returns:
%  WFD      ... functional data basis object defining final density
%  BETA     ... final model coefficients
%  FSTR     ... Struct object containing
%               FSTR.F    ... final log likelihood
%               FSTR.NORM ... final norm of gradient
%  ITERNUM  ... Number of iterations
%  ITERHIST ... History of iterations

%  last modified 1 April 2004

if nargin < 2
    error('WFDPAR is not supplied.');
end

%  check WfdPar

if ~isa_fdPar(WfdPar) 
    if isa_fd(WfdPar) | isa_basis(WfdPar)
        WfdPar = fdPar(WfdPar);
    else
        error(['WFDPAR is not a functional parameter object, ', ...
                'not a functional data object, and ', ...
                'not a basis object.']);
    end
end

%  set up WFDOBJ

Wfdobj = getfd(WfdPar);

%  set up LFDOBJ

Lfdobj = getLfd(WfdPar);
Lfdobj = int2Lfd(Lfdobj);
nderiv = getnderiv(Lfdobj);

%  set up LAMBDA

lambda = getlambda(WfdPar);

%  set up BASIS

basis  = getbasis(Wfdobj);
nbasis = getnbasis(basis);
rangex = getbasisrange(basis);

nobs   = length(y);

%  set some default arguments

if nargin < 9, dbglev  = 1;         end
if nargin < 8, active  = 2:nbasis;  end
if nargin < 7, iterlim = 20;        end
if nargin < 6, conv    = 1e-2;      end

%  initialize some arrays

climit    = [-50,0;0,50]*ones(2,nbasis);
cvec0     = getcoef(Wfdobj);
hmat      = zeros(nbasis,nbasis);
inactive  = ones(1,nbasis);
inactive(active) = 0;
inactive  = find(inactive);
ninactive = length(inactive);
dbgwrd    = dbglev > 1;

ind1      = 1:nbasis;

%  Set up linear model

covwrd = ~isempty(beta0);
if covwrd
    %  Independent variables are present ... estimate linear model
    if size(zmat,1) ~= nobs 
        error('ZMAT must have as many rows as length(X)');
    end
    ncov = size(zmat,2);
    Zsum = sum(zmat)';
    res0 = (y - zmat * beta0);
    ind2 = (nbasis+1):(nbasis+ncov);
else
    %  No independent variables are present
    ncov = 0;
    zmat = [];
    Zsum = [];
    res0 = y;
end

%  bring residuals out of range to range and set 
%    U0 = res./sigma0;

[res0, U0, indlo, indhi] = reschk(res0, rangex, sigma0);

%  initialize matrix Kmat defining penalty term

if lambda > 0 
    Kmat = lambda.*eval_penalty(basis, Lfdobj);
end

%  evaluate log likelihood
%    and its derivatives with respect to these coefficients

[logl, Dlogl] = loglfnLM(basis, cvec0, U0, zmat, sigma0);

%  compute initial badness of fit measures

Foldstr.f  =  -logl;
gvec       = -Dlogl;
if lambda > 0 
    gvec(ind1) = gvec(ind1) + 2.*(Kmat * cvec0);
    Foldstr.f  = Foldstr.f  + cvec0' * Kmat * cvec0;
end
if covwrd
    gvec(ind2)  = gvec(ind2) - 2.*sum(res0).*Zsum;
end
if ninactive > 0, gvec(inactive) = 0;  end
Foldstr.norm = sqrt(mean(gvec.^2));

%  compute the initial expected Hessian

hmat = VarfnLM(basis, cvec0, U0, zmat, sigma0);
if lambda > 0 
    hmat(ind1,ind1) = hmat(ind1,ind1) + 2.*Kmat;
end
if covwrd
    hmat(ind2,ind2) = hmat(ind2,ind2) + 2.*Zsum*Zsum';
end
if ninactive > 0 
    hmat(inactive,:) = 0;
    hmat(:,inactive) = 0;
    hmat(inactive,inactive) = eye(ninactive);
end

%  evaluate the initial update vector for correcting the initial bmat

deltac    = -hmat\gvec;
cosangle  = -gvec'*deltac/sqrt(sum(gvec.^2)*sum(deltac.^2));

%  initialize iteration status arrays

iternum = 0;
status = [iternum, Foldstr.f, -logl, Foldstr.norm];
if dbglev > 0
    fprintf('\nIteration  Criterion  Neg. Log L  Grad. Norm\n')  
    fprintf('\n%5.f     %10.4f %10.4f %10.4f\n', status);
end
iterhist = zeros(iterlim+1,length(status));
iterhist(1,:) = status;

%  quit if ITERLIM == 0

if iterlim == 0
    Fstr = Foldstr;
    iterhist = iterhist(1,:);
    beta = beta0;
    if length(indlo) > 0
        warning([num2str(length(indlo)),' lower residuals trimmed']);
    end
    indhi = find(res > rangex(2)*sigma0);
    if length(indhi) > 0
        warning([num2str(length(indhi)),' upper residuals trimmed']);
    end
    return;
end

%  -------  Begin iterations  -----------

STEPMAX = 5;
MAXSTEP = 100;
trial   = 1;
cvec    = cvec0;
beta    = beta0;
linemat = zeros(3,5);

for iter = 1:iterlim
    iternum = iternum + 1;
    Fstr    = Foldstr;
    %  set initial switches
    dblwrd = [0,0]; limwrd = [0,0]; stpwrd = 0; ind = 0; ips = 0;
    %  normalize search direction vector
    sdg     = sqrt(sum(deltac.^2));
    deltac  = deltac./sdg;
    %  compute initial slope
    linemat(2,1) = sum(deltac.*gvec);
    %  return with error condition if initial slope is nonnegative
    if linemat(2,1) >= 0
        fprintf('Initial slope nonnegative.\n');
        ind = 3;
        iterhist = iterhist(1:(iternum+1),:);
        break;
    end
    %  return successfully if initial slope is very small
    if linemat(2,1) >= -1e-5;
        if dbglev > 1, fprintf('Initial slope too small\n'); end
        iterhist = iterhist(1:(iternum+1),:);
        break;
    end
    %  load up initial search matrix 
    linemat(1,1:4) = 0;
    linemat(2,1:4) = linemat(2,1);
    linemat(3,1:4) = Foldstr.f;
    %  output initial results for stepsize 0
    stepiter  = 0;
    if dbglev > 1
        fprintf('      %3.f %10.4f %10.4f %10.4f\n', ...
            [stepiter, linemat(:,1)']); 
    end
    %  first step set to trial
    linemat(1,5)  = trial;
    %  Main iteration loop for linesrch
    for stepiter = 1:STEPMAX
        %  ensure that step does not go beyond limits on parameters
        %  check the step size
        [linemat(1,5),ind,limwrd] = ...
            stepchk(linemat(1,5), cvec, deltac, limwrd, ind, ...
            climit, active, dbglev);
        if linemat(1,5) <= 1e-9 
            %  Current step size too small ... terminate
            Fstr    = Foldstr;
            cvecnew = cvec;
            betanew = beta;
            gvecnew = gvec;
            if dbglev > 1
                fprintf('Stepsize too small: %10.4f', linemat(1,5));
            end
            break;
        end
        cvecnew = cvec + linemat(1,5).*deltac(ind1);
        if covwrd
            betanew = beta + linemat(1,5).*deltac(ind2);
            resnew  = y - zmat * betanew;
        else
            resnew  = y;
        end
        %  compute new function value and gradient
        [resnew, Unew, indlo, indhi] = reschk(resnew, rangex, sigma0);
        [logl, Dlogl]  = loglfnLM(basis, cvecnew, Unew, zmat, sigma0);
        Fstr.f  =  -logl;
        gvecnew = -Dlogl;
        if lambda > 0 
            gvecnew(ind1) = gvecnew(ind1) + 2.*Kmat * cvecnew;
            Fstr.f = Fstr.f + cvecnew' * Kmat * cvecnew;
        end
        if covwrd
            gvecnew(ind2) = gvecnew(ind2) - 2.*sum(resnew).*Zsum;
        end
        if ninactive > 0, gvecnew(inactive) = 0;  end
        Fstr.norm  = sqrt(mean(gvecnew.^2));
        %  update search matrix
        linemat(2,5) = sum(deltac.*gvecnew);
        linemat(3,5) = Fstr.f;
        %  output current results
        if dbglev > 1 
            fprintf('      %3.f %10.4f %10.4f %10.4f\n', ...
                [stepiter, linemat(:,5)']); 
        end
        %  compute next step
        [linemat,ips,ind,dblwrd] = ...
            stepit(linemat, ips, ind, dblwrd, MAXSTEP, dbglev);
        trial  = linemat(1,5);
        %  ind == 0 implies convergence
        if ind == 0 | ind == 5, break; end
        %  end iteration loop
    end
    
    %  update current parameter vectors
    
    cvec   = cvecnew;
    gvec   = gvecnew;
    Wfdobj = putcoef(Wfdobj, cvec);
    if covwrd
        beta = betanew;
        res  = y - zmat * beta;
    else
        res = y;
    end
    %  check residuals and truncate if needed
    [res, U, indlo, indhi] = reschk(res, rangex, sigma0);
    %  update and output iteration status
    status = [iternum, Fstr.f, -logl, Fstr.norm];
    iterhist(iter+1,:) = status;
    fprintf('%5.f     %10.4f %10.4f %10.4f\n', status);
    %  test for convergence
    if abs(Fstr.f-Foldstr.f) < conv
        iterhist = iterhist(1:(iternum+1),:);
        if length(indlo) > 0
            warning([num2str(length(indlo)),' lower residuals trimmed']);
        end
        indhi = find(res > rangex(2)*sigma0);
        if length(indhi) > 0
            warning([num2str(length(indhi)),' upper residuals trimmed']);
        end
        break;
    end
    %  exit loop if convergence
    if Fstr.f >= Foldstr.f,  
        break;  
    end
    %  compute the new Hessian
    hmat = VarfnLM(basis, cvec, U, zmat, sigma0);
    if lambda > 0
        hmat(ind1,ind1) = hmat(ind1,ind1) + 2.*Kmat;
    end
    if covwrd
        hmat(ind2,ind2) = hmat(ind2,ind2) + 2.*Zsum*Zsum';
    end
    if ninactive > 0 
        hmat(inactive,:) = 0;
        hmat(:,inactive) = 0;
        hmat(inactive,inactive) = eye(ninactive);
    end
    %  evaluate the update vector
    deltac    = -hmat\gvec;
    cosangle  = -gvec'*deltac/sqrt(sum(gvec.^2)*sum(deltac.^2));
    if cosangle < 0
        if dbglev > 1, disp('cos(angle) negative');  end
        deltac = -gvec;
    end
    Foldstr = Fstr;    
end

%  ------------------------------------------------------

function [logl, Dlogl] = loglfnLM(basis, cvec, res, zmat, sigma)
nbasis  = getnbasis(basis);
covwrd = ~isempty(zmat);
if covwrd
    [nobs, ncov] = size(zmat);
else
    nobs = length(res);
    ncov = 0;
end
npar    = nbasis + ncov;
ind1    = 1:nbasis;
oneobs  = ones(nobs,1);
if covwrd
    ind2   = nbasis+(1:ncov);
    onecov = ones(ncov,1);
end

phimat = getbasismatrix(res, basis);
cval   = normalize_phi(basis, cvec);
logl   = sum(phimat * cvec  - log(cval*sigma));

Dlogl = zeros(npar,1);
DcW   = phimat;
EDcW  =  oneobs*expect_phi(basis, cvec, cval);
Dlogl(ind1) = sum(DcW - EDcW);
if covwrd
    Dphimat = getbasismatrix(res, basis, 1);
    Dphivec = Dphimat * cvec;
    U2      = expect_phi(basis, cvec, cval, 1);
    DdW     = -(Dphivec*onecov').*zmat./sigma;
    EDdW    = -zmat.*(U2*cvec)./sigma;
    Dlogl(ind2) = sum(DdW - EDdW);
end

%  ------------------------------------------------------

function  Varphi = VarfnLM(basis, cvec, res, zmat, sigma) 
nbasis = getnbasis(basis);
ind1   = 1:nbasis;
covwrd = ~isempty(zmat);
if covwrd 
    [nobs, ncov] = size(zmat);
    ind2   = nbasis+(1:ncov);
    onecov = ones(ncov,1);
else
    nobs = length(res);
    ncov = 0;
end
npar   = nbasis + ncov;
oneobs = ones(nobs,1);

cval = normalize_phi(basis, cvec);
EDcW = oneobs*expect_phi(basis, cvec, cval);
if covwrd
    U2   = expect_phi(basis, cvec, cval, 1);
    EDdW = -zmat.*(U2*cvec)./sigma;
    EDw  = [EDcW,EDdW];
else
    EDw = EDcW;
end

EDwDwt = zeros(npar,npar);
U11 = squeeze(expect_phiphit(basis, cvec, cval));
EDwDwt(ind1,ind1) = nobs.*U11;  
if covwrd
    U21 = squeeze(expect_phiphit(basis, cvec, cval, 1, 0));
    U22 = squeeze(expect_phiphit(basis, cvec, cval, 1, 1));
    EDwDwt(ind2,ind1) = -sum(zmat)' * (cvec'*U21)./sigma;
    EDwDwt(ind1,ind2) = EDwDwt(ind2,ind1)';
    EDwDwt(ind2,ind2) = (cvec'*U22*cvec) .* (zmat'*zmat)./sigma.^2;
end

Varphi = EDwDwt - EDw'*EDw;

%  ------------------------------------------------------

function  [res, U, indlo, indhi] = reschk(res, rangex, sigma0) 
%RESCHK brings residuals outside of limits to limits
indlo = find(res < rangex(1)*sigma0);
if length(indlo) > 0
    res(indlo) = rangex(1)*sigma0;
end
indhi = find(res > rangex(2)*sigma0);
if length(indhi) > 0
    res(indhi) = rangex(2)*sigma0;
end
U   = res./sigma0;

%  ------------------------------------------------------

function ss = expect_phi(basisfd, cvec, Cval, nderiv, ...
                         rng, JMAX, EPS)  
%  Computes expectations of basis functions with respect to density
%      p(x) = Cval^{-1} exp c'phi(x)
%  by numerical integration using Romberg integration

%  Arguments:
%  BASISFD  ...  A basis function
%           object.  In the latter case, a functional data object
%           is created from a basis function object by using the
%           identity matrix as the coefficient matrix.
%           The functional data objects must be univariate.
%  CVEC ... coefficient vector defining density, of length NBASIS
%  CVAL ... normalizing constant defining density
%  MU   ... mean value to be subtracted from variates
%  SIGMA .. standard deviation to define u = (x - mu)/sigma
%  RNG ...  vector of length 2 giving the interval over which the
%           integration is to take place
%  NDERIV . order of derivative required for basis function expectation
%  UWRD ... If T, expectation is of (D PHI)*U
%  JMAX ... maximum number of allowable iterations
%  EPS  ... convergence criterion for relative error

%  Return:
%  A vector SS of length NBASIS of integrals of functions.

  %  check arguments, and convert basis objects to functional data objects

  if ~strcmp(class(basisfd),'basis')
    error('First argument must be a basis function object.');
  end

  nbasis = getnbasis(basisfd);
  oneb   = ones(1,nbasis);
  rangeval = getbasisrange(basisfd);

  %  set default arguments

  if nargin < 7, EPS = 1e-7;     end
  if nargin < 6, JMAX = 15;      end
  if nargin < 5, rng = rangeval; end
  if nargin < 4, nderiv = 0;     end
  if nargin < 3, Cval = 1;       end

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
  fx = getbasismatrix(x, basisfd);
  wx = fx * cvec;
  wx(wx < -50) = -50;
  px = exp(wx)./Cval;
  if nderiv == 0 
    Dfx = fx;
  else 
    Dfx = getbasismatrix(x, basisfd, 1);
  end
  sumj = Dfx' * px;
  smat(1,:)  = width.*sumj'./2;
  tnm = 0.5;
  j   = 1;

  %  now iterate to convergence

  for j = 2:JMAX
    tnm  = tnm*2;
    del  = width/tnm;
    if j == 2 
      x = (rng(1) + rng(2))/2; 
    else 
      x = (rng(1)+del/2 : del : rng(2)-del/2)';
    end
    nx = length(x);
    fx = getbasismatrix(x, basisfd);
    wx = fx * cvec;
    wx(wx < -50) = -50;
    px = exp(wx)./Cval;
    if nderiv == 0 
      Dfx = fx; 
    else 
      Dfx = getbasismatrix(x, basisfd, 1);
    end
    sumj = Dfx' * px;
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

  %  ------------------------------------------------------

function ss = expect_phiphit(basisfd, cvec, Cval, nderiv1, nderiv2, ...
                             rng, JMAX, EPS) 

%  Computes expectations of cross product of basis functions with
%  respect to density
%      p(x) = Cval^{-1} exp c'phi(x)
%  by numerical integration using Romberg integration

%  Arguments:
%  BASISFD  ...  A basis function
%           object.  In the latter case, a functional data object
%           is created from a basis function object by using the
%           identity matrix as the coefficient matrix.
%           The functional data objects must be univariate.
%  CVEC ... coefficient vector defining density
%  CVAL ... normalizing constant defining density
%  RNG ...  vector of length 2 giving the interval over which the
%           integration is to take place
%  JMAX ... maximum number of allowable iterations
%  EPS  ... convergence criterion for relative error

%  Return:
%  A matrix of order NBASIS of integrals of functions.

  %  check arguments, and convert basis objects to functional data objects

  if ~strcmp(class(basisfd),'basis') 
    error('First argument must be a basis function object.');
  end

  nbasis = getnbasis(basisfd);
  oneb   = ones(1,nbasis);
  rangeval = getbasisrange(basisfd);

  %  set default arguments

  if nargin < 8, EPS = 1e-7;       end
  if nargin < 7, JMAX = 9;         end
  if nargin < 6, rng = rangeval;   end
  if nargin < 5, nderiv2 = 0;      end
  if nargin < 4, nderiv1 = 0;      end
  if nargin < 3, Cval = ones(n,1); end

  %  set up first iteration

  width = rng(2) - rng(1);
  JMAXP = JMAX + 1;
  h = ones(JMAXP,1);
  h(2) = 0.25;
  %  matrix SMAT contains the history of discrete approximations to the integral
  smat = zeros([JMAXP,nbasis,nbasis]);
  %  the first iteration uses just the endpoints
  x  = rng';
  nx = length(x);
  fx = getbasismatrix(x, basisfd);
  wx = fx * cvec;
  wx(wx < -50) = -50;
  px = exp(wx)./Cval;
  if nderiv1 == 0 
    Dfx1 = fx; 
  else 
    Dfx1 = getbasismatrix(x, basisfd, 1);
  end
  if nderiv2 == 0 
    Dfx2 = fx; 
  else 
    Dfx2 = getbasismatrix(x, basisfd, 1);
  end
  sumj = Dfx1' * ((px * oneb) .* Dfx2);
  smat(1,:,:)  = width.*sumj./2;
  tnm = 0.5;
  j   = 1;

  %  now iterate to convergence
  for j = 2:JMAX
    tnm  = tnm*2;
    del  = width/tnm;
    if j == 2 
      x = (rng(1) + rng(2))/2; 
    else 
      x = (rng(1)+del/2 : del : rng(2)-del/2)';
    end
    nx = length(x);
    fx = getbasismatrix(x, basisfd);
    wx = fx * cvec;
    wx(wx < -50) = -50;
    px = exp(wx)./Cval;
    if nderiv1 == 0 
      Dfx1 = fx;
    else 
      Dfx1 = getbasismatrix(x, basisfd, 1);
    end
    if nderiv2 == 0 
      Dfx2 = fx; 
    else 
      Dfx2 = getbasismatrix(x, basisfd, 1);
    end
    sumj = Dfx1' * ((px * oneb) .* Dfx2);
    smat(j,:,:) = (squeeze(smat(j-1,:,:)) + width.*sumj./tnm)./2;
    if j >= 5
      ind = (j-4):j;
      temp = squeeze(smat(ind,:,:));
      [ss, dss] = polintarray(h(ind),temp,0);
      if ~any(abs(dss) > EPS.*max(max(abs(ss))))
        %  successful convergence
        return;
      end
    end
    smat(j+1,:,:) = smat(j,:,:);
    h(j+1) = 0.25*h(j);
  end

