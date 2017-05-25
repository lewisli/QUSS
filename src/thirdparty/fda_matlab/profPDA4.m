function [SSE, PENSSE, penmat, fdobj, df, gcv] = ...
    profPDA4(bvec, y, fitcell, bwtcell, awtcell, ufdcell)
% profPDA4 estimates a homogeneous linear differential equation
%  directly from discrete noisy observations of a process.  
%  In this version forcing functions can be accommodated. 
%
%profPDA4 works with the basis function expansions of the
%  estimates of the coefficient functions a_k(t) and b_j(t) 
%  in the possibly nonhomogeneous linear differential operator
%
%    Lx(t) = a_1(t)u_1(t) + ... + a_k(t)u_K(t) + 
%       b_0(t)x(t) + ... + b_{m-1}(t)D^{m-1}x(t) + 
%       \exp(b_m(t))D^m x(t)
%
%  of order m = NORDER that minimizes in a least squares sense the residual
%  functions f(t) = Lx(t).  
%
%  Arguments:
%  BVEC       ... A super-vector containing vectors of coefficients defining 
%                 the weight functions to be estimated.  
%                 Coefficients are stacked on top of each other as follows:
%                 coefficients for b_0 first, b_1 next, and so on, and
%                 then continues on to the forcing function coefficients
%                 to be estimated.  See BVEC2LFD for more details.
%  Y          ... Matrix of curve values to be fit.
%  FITCELL    ... Cell object with one cell per variable.
%                 The contents of each cell are a struct containing
%    BASISOBJ ... A basis object for representing the curve(s)
%                   that are the solutions of the differential equations
%                   that fit the data
%    BMAT0    ... Weighted cross-product matrix for basis matrix values.
%    DMAT0    ... Weighted product of transposed basis matrix and Y.
%    LAMBDA   ... The smoothing parameter controlling the penalty on
%                   the roughness penalty.  In order to estimate the 
%                   DIFE, this should be large but not too large.  
%                   If it is too small, the roughness is ignored, and
%                   the differential equation is badly estimated.
%                   If it is too large, though, the data are ignored and
%                   the DIFE is also badly estimated.  
%  BWTCELL    ... Cell object for the weight functions for the 
%                 homogeneous part of the equation.
%  AWTCELL    ... Cell object for the weight functions for the 
%                 forcing functions
%  UFDCELL    ... Cell object containing functional data objects for
%                 the forcing functions
%  Returns:
%  PENSSE  ...  The penalized error sum of squares.  This is what is
%               required for numerical optimization software since it
%               is the objective function being minimized.
%               It is PENSSE = SSE + lambda.*C'KC, where C is the 
%               coefficient vector or matrix defining the curve(s) 
%               fitting the data and K is the penalty matrix corresponding
%               to the estimated DIFE.                 
%  SSE     ...  The error sum of squares.
%  PENALTY ...  The quantity C'KC'.
%  FDOBJ   ...  Functional data object fitting the data

%  Last modified 7 January 2004

%  check number of variables

nvar = length(fitcell);

if nvar > 1 & length(bwtcell) > 1
    error('This version cannot handle multiple variables.');
end

%  check basis

if ~isa_basis(basisobj)
    error('basisobj is not a basis object.');
end

n        = size(y,1);
nbasis   = getnbasis(basisobj);
onebasis = ones(1,nbasis);
ncurves  = size(Dmat,2);

%  check LAMBDA

if lambda < 0
    warning ('Value of LAMBDA was negative, and 0 used instead.');
    lambda = 0;
end

%  calculate Lfd object

Lfdobj = bvec2Lfd(bvec, bwtcell, awtcell, ufdcell);
nderiv = getnderiv(Lfdobj);
    
%  smooth the data

afdcell = getafd(Lfdobj);  % multiplier(s) of forcing function(s)
ufdcell = getufd(Lfdobj);  % forcing function(s)
if ~isempty(afdcell) & ~isempty(ufdcell)
    %  here the linear differential operator is not homogeneous.  
    %  first set up the homogeneous counterpart
    Lfdhom = Lfd(nderiv, getwfd(Lfdobj));
    %  evaluate the penalty matrix for the homogeneous operator
    penmat = eval_penalty(basisobj, Lfdhom);
    %  set up the part of the roughness penalty affected by the
    %  presence of forcing function(s)
    penvec = zeros(nbasis,1);
    nforce = length(ufdcell);
    for k=1:nforce
        afdk = afdcell{k};
        ufdk = ufdcell{k};
        ffdk = times(afdk,ufdk,basisobj);
%         ffdk = times(afdk,ufdk,getbasis(ufdk));
        penvec = penvec - inprod(basisobj, ffdk, Lfdhom, int2Lfd(0));
    end
else
    %  here the linear differential operator is homogeneous
    %  only the penalty matrix is needed.
    nforce = 0;
    penmat = eval_penalty(basisobj, Lfdobj);
    penvec = zeros(nbasis,1);
end

%  check value of LAMBDA for overflow problems

Bnorm   = sqrt(sum(sum(Bmat.^2)));
pennorm = sqrt(sum(sum(penmat.^2)));
condno  = pennorm/Bnorm;
if lambda*condno > 1e12
    lambda = 1e12/condno;
    warning(['LAMBDA reduced to ',num2str(lambda),...
            ' to prevent overflow']);
end

%  set up coefficient matrix BMAT

Bmat = Bmat0;
if lambda > 0
    Bmat = Bmat + lambda .* penmat;
end
%  set up right side vector DMAT
Dmat = Dmat0;
if ~all(penvec == 0) & lambda > 0
    %  if the linear differential operator is nonhomogeneous
    %  use PENVEC to alter the right side of the equation.
    Dmat = Dmat + lambda.*(penvec*ones(1,ncurves));
end

%  compute inverse of Bmat

if is_diag(Bmat)
    Bmatinv = diag(1./diag(Bmat));
else
    Bmatinv = inv(Bmat);
end

%  compute degrees of freedom of smooth

df = sum(diag(Bmatinv*Bmat0));

%  solve normal equations for each observation

coef = Bmatinv * Dmat;

%  compute error sum of squares

SSE = 0;
for ivar=1:nvar
    coefi = coef(:,ivar);
    yi    = y(:,ivar);
    SSE = SSE + sum(yi.^2) - 2*.Dmat0(:,ivar)'*coefi + ...
        coefi'*Bmat0(:,:,ivar)*coefi;
end
yhat = basismat * coef;
SSE  = sum(sum((y - yhat).^2));

%  compute  GCV index

if df < n
    gcv = (SSE/n)/((n - df)/n)^2;
else
    gcv = NaN;
end

fdobj = fd(coef, basisobj);

%  penalized least squares

PENSSE = SSE + lambda.*sum(diag(coef'*penmat*coef));

%  update penalized least squares for terms for 
%  smoothing derivative weight coefficients

m2 = 0;
for j=1:nderiv
    bwtstructj = bwtcell{j};
    if getestimate(bwtstructj) 
        basisj  = getbasis(getfd(bwtstructj));
        nbasisj = getnbasis(basisj);
        m1 = m2 + 1;
        m2 = m2 + nbasisj;
        lambdaj = getlambda(bwtstructj);
        if lambdaj > 0
            Lfdobjj = getLfd(bwtstructj);
            penmatj = eval_penalty(basisj, Lfdobjj);
            bcoefj  = bvec(m1:m2);
            termj   = lambdaj.*bcoefj'*penmatj*bcoefj;
            SSE     = SSE + termj;
        end
    end
end

%  update penalized least squares for terms for 
%  smoothing forcing function coefficients

for k=1:nforce
    awtstructk = awtcell{k};
    if getestimate(awtstructk) 
        basisk  = getbasis(getfd(awtstructk));
        nbasisk = getnbasis(basisk);
        m1 = m2 + 1;
        m2 = m2 + nbasisk;
        lambdak = getlambda(awtstructk);
        if lambdak > 0
            Lfdobjk = getLfd(awtstructk);
            penmatk = eval_penalty(basisk, Lfdobjk);
            bcoefk  = bvec(m1:m2);
            SSE     = SSE + lambdak.*bcoefk'*penmatk*bcoefk;
        end
    end
end

%  ----------------------------------------------------------------

function [Lfdobj, bfdcell, afdcell] = ...
    bvec2Lfd(bvec, bwtcell, awtcell, ufdcell)
%BVEC2LFD sets up linear differential operator object from a
%  super-vector of coefficients for expansions of coefficient functions
%  to be estimated.  The first coefficient vector in BVEC is for the
%  first homogeneous coefficient to be estimated, and so on, 
%  progressing through remaining homogeneous coefficients, and then
%  moving to the first forcing function to be estimated, and so on.

%  Last modified 18 August 2003

%  calculate order

norder = size(bwtcell,length(size(bwtcell)));

%  calculate number of forcing functions

if isempty(awtcell) | isempty(ufdcell)
    nforce = 0;
    afdcell = {};
else
    nforce = size(ufdcell,2);
end

%  loop through derivatives

m2 = 0;
for ideriv=1:norder
    fdobji = getfd(bwtcell{ideriv});
    if getestimate(bwtcell{ideriv})
        psi_basis  = getbasis(fdobji);
        npsi_basis = getnbasis(psi_basis);
        m1 = m2 + 1;
        m2 = m2 + npsi_basis;
        fdobji = putcoef(fdobji, bvec(m1:m2));
    end
    bfdcell{ideriv} = fdobji;
end

%  loop through forcing functions

for iforce=1:nforce
    fdobji = getfd(awtcell{iforce});
    if getestimate(awtcell{iforce})
        basisa  = getbasis(fdobji);
        nbasisa = getnbasis(basisa);
        m1 = m2 + 1;
        m2 = m2 + nbasisa;
        fdobji = putcoef(fdobji, bvec(m1:m2));
    end
    afdcell{iforce} = fdobji;
end
    
%  compute differential operator that defines the penalty

if isempty(awtcell) | isempty(ufdcell)
    Lfdobj = Lfd(norder, bfdcell);
else
    Lfdobj = Lfd(norder, bfdcell, afdcell, ufdcell);
end

