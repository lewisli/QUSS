function ss = normalize_phi(basisfd, cvec, mu, sigma, rng, JMAX, EPS) 

%  Computes integrals of
%      p(x) = exp phi'(x) * cvec
%  by numerical integration using Romberg integration

%  Arguments:
%  BASISFD  ...  Basis function object with basis functions phi.
%  CVEC ... coefficient vector defining density, of length NBASIS
%  MU   ... mean values to be subtracted from variates
%  SIGMA .. standard deviation to define u = (x - mu)/sigma
%  RNG ...  vector of length 2 giving the interval over which the
%           integration is to take place.  Multiply a standard interval
%           like (-5,5) by sigma to make it scale free
%  JMAX ... maximum number of allowable iterations
%  EPS  ... convergence criterion for relative error

%  Return:
%  The integral of the function.

  %  check arguments, and convert basis objects to functional data objects

  if ~strcmp(class(basisfd), 'basis') 
    error('First argument must be a basis function object.');
  end

  nbasis = getnbasis(basisfd);
  oneb   = ones(1,nbasis);
  rangeval = getbasisrange(basisfd);

  %  set default arguments

  if nargin < 7, EPS = 1e-7;     end
  if nargin < 6, JMAX = 15;      end
  if nargin < 5, rng = rangeval; end
  if nargin < 4, sigma = 1;      end
  if nargin < 3, mu = 0;         end

  %  set up first iteration

  width = rng(2) - rng(1);
  JMAXP = JMAX + 1;
  h = ones(JMAXP,1);
  h(2) = 0.25;
  %  matrix SMAT contains the history of discrete approximations to the integral
  smat = zeros(JMAXP,1);
  %  the first iteration uses just the endpoints
  x  = rng;
  nx = length(x);
  ox = ones(nx,1);
  u  = (x - mu)./sigma;
  fx = getbasismatrix(u, basisfd);
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
    u  = (x - mu)./sigma;
    fx = getbasismatrix(u, basisfd);
    wx = fx * cvec;
    wx(wx < -50) = -50;
    px = exp(wx);
    smat(j) = (smat(j-1) + width.*sum(px)./tnm)./2;
    if j >= 5 
      ind = (j-4):j;
      [ss, dss] = polintarray(h(ind),smat(ind),0);
      if ~any(abs(dss) >= EPS.*max(abs(ss))) 
        %  successful convergence
        return;
      end
    end
    smat(j+1) = smat(j);
    h(j+1)   = 0.25*h(j);
  end
  warning(['No convergence after ',num2str(JMAX),' steps in NORMALIZE.PHI'])

