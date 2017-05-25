function  ss = inprod_n(fdstr1, fdstr2, Lfd1, Lfd2, wtfd, JMAX, EPS)
%  INPROD   Computes matrix of inner products of functions
%    by numerical integration
%    using Romberg integration with the trapezoidal rule.
%  Arguments:
%  FD1STR and FDSTR2 ...  these may be either functional data or basis function
%                    objects.  In the latter case, a functional data object
%                    is created from a basis function object by using the
%                    identity matrix as the coefficient matrix.
%                    Both functional data objects must be univariate.
%                    If  inner products for multivariate objects are needed,
%                    use a loop and call inprod(fdstr1(i),fdstr2(i)).
%  Lfd1 and Lfd2 ...  dif ferential operators for inner product for
%                    FD1 and FD2, respectively
%  WTFD ...  A functional data object defining a weight
%  JMAX ...  maximum number of allowable iterations
%  EPS  ...  convergence criterion for relative error
%  Return:
%  A NREP1 by NREP2 matrix SS of inner products for each possible pair
%  of functions.

%  last modified 26 April 1999

  %  set up default values of arguments

  if  nargin < 7
    EPS = 1e-4;
  end
  if  nargin < 6
    JMAX = 15;
  end
  if nargin < 5
    wtfd = 0;
  end
  if  nargin < 4
    Lfd2 = 0;
  end
  if  nargin < 3
    Lfd1 = 0;
  end

 JMIN = 11;
 
  %  check arguments, and convert basis objects to functional data objects

  if ~isa_Lfd(Lfd1)
    error (['Argument Lfd1 is neither a functional data object', ...
             ' nor an integer.']);
  end

  if ~isa_Lfd(Lfd2)
    error (['Argument Lfd2 is neither a functional data object', ...
             ' nor an integer.']);
  end
  
  if  isa_fd(wtfd)
    coefw = getcoef(wtfd);
    coefd = size(coef);
    if coefd(2) > 1
      error('Argument WTFD is not a single function');
    end
  end
  
  fdclass = 1;
  if  isa_fd(fdstr1)
    coef1 = getcoef(fdstr1);
  elseif isa_basis(fdstr1)
    coef1  = eye(getnbasis(fdstr1));
    temp1  = fd(coef1, fdstr1);
    fdstr1 = temp1;
  else
    fdclass = 0;
  end

  if  isa_fd(fdstr2)
    coef2 = getcoef(fdstr2);
  elseif isa_basis(fdstr2)
    coef2  = eye(getnbasis(fdstr2));
    temp2  = fd(coef2, fdstr2);
    fdstr2 = temp2;
  else
    fdclass = 0;
  end

  if ~fdclass
    error (['The two first arguments are not', ...
            ' functional data or basis objects.']);
  end

  %  determine NREP1 and NREP2, and check for common range

  coefd1 = size(coef1);
  coefd2 = size(coef2);
  if length(coefd1) > 2 | length(coefd2) > 2
    error('Functional data objects must be univariate');
  end
  nrep1 = coefd1(2);
  nrep2 = coefd2(2);
  basisobj1 = getbasis(fdstr1);
  basisobj2 = getbasis(fdstr2);
  range1 = getbasisrange(basisobj1);
  range2 = getbasisrange(basisobj2);
  if  any(range1-range2) ~= 0
    error('Ranges are not equal');
  end

  %  set up first iteration

  width = range1(2) - range1(1);
  JMAXP = JMAX + 1;
  h = ones(JMAXP,1);
  h(2) = 0.25;
  s = reshape(zeros(JMAXP*nrep1*nrep2,1),[JMAXP,nrep1,nrep2]);
  %  the first iteration uses just the endpoints
  fx1 = eval_fd(fdstr1, range1, Lfd1);
  fx2 = eval_fd(fdstr2, range1, Lfd2);
  if ~isnumeric(wtfd)
    wtd = eval_fd(wtfd, range1);
    fx2 = (wtd * ones(1,nrep2)) .* fx2
  end
  s(1,:,:)  = width.*(fx1' * fx2)./2;
  tnm = 0.5;

  %  now iterate to convergence

  for j = 2:JMAX
    tnm = tnm.*2;
    del = width./tnm;
    x   = range1(1)+del/2:del:range1(2)-del/2;
    fx1 = eval_fd(fdstr1, x, Lfd1);
    fx2 = eval_fd(fdstr2, x, Lfd2);
    if ~isnumeric(wtfd)
      wtd = eval_fd(wtfd, x);
      fx2 = (wtd * ones(1,nrep2)) .* fx2
    end
    chs = reshape(width.*(fx1' * fx2)./tnm,[1,nrep1,nrep2]);
    s(j,:,:) = (s(j-1,:,:) + chs)./2;
    if j >= 5
      ind = (j-4):j;
      ya = s(ind,:,:);
      xa = h(ind);
      absxa = abs(xa);
      [absxamin, ns] = min(absxa);
      cs = ya;
      ds = ya;
      y  = squeeze(ya(ns,:,:));
      ns = ns - 1;
      for m = 1:4
        for i = 1:(5-m)
          ho      = xa(i);
          hp      = xa(i+m);
          w       = (cs(i+1,:,:) - ds(i,:,:))./(ho - hp);
          ds(i,:,:) = hp.*w;
          cs(i,:,:) = ho.*w;
        end
        if 2*ns < 5-m
          dy = squeeze(cs(ns+1,:,:));
        else
          dy = squeeze(ds(ns,:,:));
          ns = ns - 1;
        end
        y = y + dy;
      end
      ss  = y;
      crit = max(max(abs(dy)))./max(max(abs(ss)));
      if crit < EPS& j >= JMIN
        return
      end
    end
    s(j+1,:,:) = s(j,:,:);
    h(j+1)   = 0.25.*h(j);
  end
  disp(['No convergence after',num2str(JMAX),' steps in INPROD']);

