function [gval, tval, fval] = mongrad(x, Wfd)
%MONGRAD evaluates the gradient of a monotone function of the form
%             h(x) = (D^{-1} exp Wfd)(x)
%  where  D^{-1} means taking the indefinite integral.
%  The interval over which the integration takes places is defined in
%  the basisobj object = Wfd.
%  Wfd may contain several curves, but must only one variable.
%  If Wfd contains only one curve, use MONGRAD1 instead to avoid
%     using multidimensional arrays.  It will be faster.
%  Arguments:
%  X         ...  Vector of argument values at which gradient is evaluated
%  WFD       ...  Functional data object defining monotone function(s)
%  Returns:
%  GVAL  ... values of derivatives in NOBS by NBASIS by NCURVE array
%  TVAL  ... Arguments used for trapezoidal approximation to integral
%  FVAL  ... Values of exp Wfd corresponding to TVAL

%  Last modified 10 December 2000

  if nargin < 2,  error('There are less than three arguments');  end
 
  %  set some constants
  
  EPS  = 1e-5;  JMIN = 11;  JMAX   = 15;
  
  %  get coefficient matrix and dimensions of problem 
  
  coef  = getcoef(Wfd);
  if ndims(coef) > 2  
    error('WFD is not a univariate function');
  end
  ncurve = size(coef,2);

  basisobj = getbasis(Wfd);
  rangeval = getbasisrange(basisobj);
  nbasis   = getnbasis(basisobj);
  onebas   = ones(1,nbasis);

  %  set up first iteration

  width = rangeval(2) - rangeval(1);
  JMAXP = JMAX + 1;
  h     = ones(JMAXP,1);
  h(2)  = 0.25;
  %  matrix SMAT contains the history of discrete approximations to the
  %    integrals
  smat = zeros(JMAXP,nbasis,ncurve);
  %  array TVAL contains the argument values used = the approximation
  %  array FVAL contains the integral values at these argument values,
  %     rows corresponding to argument values
  %  the first iteration uses just the endpoints
  tval = rangeval;
  bmat = eval_basis(basisobj, rangeval);
  Dfx0 = exp(bmat*coef);
  grad = zeros(2,nbasis,ncurve);
  for icurve=1:ncurve
     grad(:,:,icurve) = (Dfx0(:,icurve)*onebas).*bmat;
  end
  smat(1,:,:) = (width/2).*sum(grad);
  fval = grad;
  tnm = 0.5;
  j   = 1;
  %  now iterate to convergence
  for j = 2:JMAX
    tnm  = tnm*2;
    del  = width/tnm;
    hdel = del/2;
    tj   = (rangeval(1)+hdel):del:rangeval(2);
    ntj  = length(tj);
    tval = [tval, tj];
    bmat = eval_basis(basisobj, tj);
    Dfxj = exp(bmat*coef);
    grad = zeros(ntj,nbasis,ncurve);
    for icurve=1:ncurve
        grad(:,:,icurve) = (Dfxj(:,icurve)*onebas).*bmat;
    end
    smat(j,:,:) = (smat(j-1,:,:) + del.*sum(grad))./2;
    fval = cat(1,fval,grad);
    if j >= max([JMIN,5])
      ind = (j-4):j;
      [ss,dss] = polintmat1(h(ind),smat(ind,:,:),0);
      if all(abs(dss) < EPS.*max(max(abs(ss))))
        % successful convergence
        % sort argument values and corresponding function values
        [tval,ordind] = sort(tval);
        fval  = fval(ordind,:,:);
        % set up partial integral values
        integfval = (tval(2) - tval(1)).*cumtrapz(fval);
        gval = zeros(length(x),nbasis,ncurve);
        for icurve=1:ncurve
           gval(:,:,icurve) = ...
              interp1(tval, integfval(:,:,icurve), x, 'cubic');
        end
        return;
      end
    end
    smat(j+1,:,:) = smat(j,:,:);
    h(j+1)      = 0.25*h(j);
  end
  error(['No convergence after ',num2str(JMAX),' steps in MONGRAD']);
