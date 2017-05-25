function [SSE, DSSE, PENSSE, fdobj, df, gcv] = ...
    profPDA(bvec, y, basisobj, basismat, Bmat, Dmat, ...
            bwtcell, awtcell, ufdcell, lambda, derivs)
% profPDA estimates a homogeneous linear differential equation
%  directly from discrete noisy observations of a process.  
%  In this version forcing functions can be accommodated, but not
%  multiple variables. 
%
%profPDA works with the basis function expansions of the
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
%  BASISOBJ   ... A basis object for representing the curve(s)
%                 that are the solutions of the differential equations
%                 that fit the data
%  BASISMAT   ... Matrix of basis matrix values.
%  BMAT       ... Weighted cross-product matrix for basis matrix values.
%  DMAT       ... Weighted product of basis matrix and Y.
%  BWTCELL    ... Cell object for the weight functions for the 
%                 homogeneous part of the equation.
%  AWTCELL    ... Cell object for the weight functions for the 
%                 forcing functions
%  UFDCELL    ... Cell object containing functional data objects for
%                 the forcing functions
%  LAMBDA     ... The smoothing parameter controlling the penalty on
%                 the roughness penalty.  In order to estimate the 
%                 DIFE, this should be large but not too large.  
%                 If it is too small, the roughness is ignored, and
%                 the differential equation is badly estimated.
%                 If it is too large, though, the data are ignored and
%                 the DIFE is also badly estimated.  
%  DERIVS     ... Returns the gradient in DSSE if nonzero, and []
%                 otherwise.
%  Returns:
%  SSE     ...  The error sum of squares.  This is what is
%               required for numerical optimization software since it
%               is the objective function being minimized.
%  DSSE    ...  The gradient of SSE.  Computed only if DERIVS is nonzero,
%               otherwise is returned as an empty matrix.
%  PENSSE  ...  The penalized error sum of squares.  
%               It is PENSSE = SSE + lambda.*C'KC, where C is the 
%               coefficient vector or matrix defining the curve(s) 
%               fitting the data and K is the penalty matrix corresponding
%               to the estimated DIFE.                 
%  FDOBJ   ...  Functional data object fitting the data
%  DF      ...  A measure of the equivalent degrees of freedom in the
%               smoothing matrix.
%  GCV     ...  The generalized cross-validation criterion

%  Last modified 16 June 2004

if nargin < 11
    derivs = 0;
end

%  check number of variables

J = size(bwtcell,1);

if J > 1 & length(bwtcell) > 1
    error('This version cannot handle multiple variables.');
end

%  get number of data points and number of curves

[n, ncurves] = size(y);
if n == 1 & ncurves > 1
    %  transpose Y if n == 1
    y = y';
    n = ncurves;
    ncurves = 1;
end

%  check basis

if ~isa_basis(basisobj)
    error('BASISOBJ is not a basis object.');
end

nbasis   = getnbasis(basisobj);
onebasis = ones(1,nbasis);
ncurves  = size(Dmat,2);

%  check LAMBDA

if lambda < 0
    warning ('Value of LAMBDA was negative, and 0 used instead.');
    lambda = 0;
end

%  transfer parameters from vector BVEC to cell objects
%  BWTCELL and AWTCELL

nderiv = size(bwtcell,length(size(bwtcell)));

%  calculate number of forcing functions

if isempty(awtcell) | isempty(ufdcell)
    nforce = 0;
    awtcell = {};
else
    nforce = size(ufdcell,2);
end

%  loop through derivatives to transfer coefficients from
%  BVEC to the appropriate derivative weight functions

m2 = 0;
for ideriv=1:nderiv
    fdPari = bwtcell{ideriv};
    fdobji = getfd(fdPari);
    if getestimate(fdPari)
        npsi_basis = getnbasis(getbasis(fdobji));
        m1 = m2 + 1;
        m2 = m2 + npsi_basis;
        fdobji = putcoef(fdobji, bvec(m1:m2));
    end
    bwtcell{ideriv} = putfd(fdPari, fdobji);
end

%  loop through forcing functions to transfer coefficients from
%  BVEC to the appropriate forcing function weight functions

for iforce=1:nforce
    afdPari = awtcell{iforce};
    afdobji = getfd(afdPari);
    if getestimate(afdPari)
        nbasisa = getnbasis(getbasis(afdobji));
        m1 = m2 + 1;
        m2 = m2 + nbasisa;
        afdobji = putcoef(afdobji, bvec(m1:m2));
    end
    awtcell{iforce} = putfd(afdPari, afdobji);
end
    
%  compute the penalty matrix and penalty vector

[penmat, penvec, DR, Ds] = ...
    eval_Rs(bwtcell, awtcell, ufdcell, basisobj, derivs);

%  check for ill conditioning in the linear equations

Bnorm   = sqrt(sum(sum(Bmat.^2)));
pennorm = sqrt(sum(sum(penmat.^2)));
condno  = pennorm/Bnorm;
if lambda*condno > 1e12
    lambda = 1e12/condno;
    warning(['LAMBDA reduced to ',num2str(lambda),...
            ' to prevent overflow']);
end

%  update the coefficient matrix with the penalty matrix

Mmat  = Bmat + lambda .* penmat;

%  update the right side vector with the penalty vector

if ~isempty(awtcell) & ~isempty(ufdcell)
    %  if the linear differential operator is nonhomogeneous
    %  use PENVEC to alter the right side of the equation.
    Dmat = Dmat + lambda.*(penvec*ones(1,ncurves));
end

%  compute inverse of Mmat

if is_diag(Mmat)
    Mmatinv = diag(1./diag(Mmat));
else
    Mmatinv = inv(Mmat);
end

%  compute degrees of freedom of smooth

df = sum(diag(Mmatinv*Bmat));

%  solve normal equations for each observation

coef = Mmatinv * Dmat;

%  compute error sum of squares

yhat = basismat * coef;
res  = y - yhat;
SSE  = sum(sum(res.^2));

%  compute gradient

if derivs
    PhiMmatinv = basismat * Mmatinv;
    m2 = 0;
    for j=1:nderiv
        bfdParj = bwtcell{j};
        if getestimate(bfdParj) 
            basisj  = getbasis(getfd(bfdParj));
            nbasisj = getnbasis(basisj);
            m1 = m2 + 1;
            m2 = m2 + nbasisj;
            for m = m1:m2
                DRm = squeeze(DR(:,:,m));
                if nforce > 0
                    DAm = lambda.*PhiMmatinv* ...
                        (DRm*coef - Ds(:,m)*ones(1,ncurves));
                else
                    DAm = lambda.*PhiMmatinv*DRm*coef;
                end
                DSSE(m) = 2.*trace(res'*DAm);
            end
        end
    end
    
    %  update penalized least squares for terms for 
    %  smoothing forcing function coefficients
    
    for k=1:nforce
        afdPark = awtcell{k};
        if getestimate(afdPark) 
            basisk  = getbasis(getfd(afdPark));
            nbasisk = getnbasis(basisk);
            m1 = m2 + 1;
            m2 = m2 + nbasisk;
            for m = m1:m2
                DAm = lambda.*PhiMmatinv*Ds(:,m);
                DSSE(m) = -2.*sum(res'*DAm);
            end
        end
    end
else
    DSSE = [];
end

%  compute  GCV index

if df < n
    gcv = (SSE/n)/((n - df)/n)^2;
else
    gcv = NaN;
end

fdobj = fd(coef, basisobj);

%  penalized least squares

PENSSE = SSE + lambda.*sum(diag(coef'*penmat*coef));

%  update least squares SSE for terms for 
%  smoothing derivative weight coefficients

m2 = 0;
for j=1:nderiv
    bfdParj = bwtcell{j};
    if getestimate(bfdParj) 
        basisj  = getbasis(getfd(bfdParj));
        nbasisj = getnbasis(basisj);
        m1 = m2 + 1;
        m2 = m2 + nbasisj;
        lambdaj = getlambda(bfdParj);
        if lambdaj > 0
            Lfdobjj = getLfd(bfdParj);
            penmatj = eval_penalty(basisj, Lfdobjj);
            bcoefj  = bvec(m1:m2);
            termj   = lambdaj.*bcoefj'*penmatj*bcoefj;
            SSE     = SSE + termj;
        end
    end
end

%  update penalized least squares for terms for 
%  smoothing forcing function coefficients

nforce = length(ufdcell);
for k=1:nforce
    afdPark = awtcell{k};
    if getestimate(afdPark) 
        basisk  = getbasis(getfd(afdPark));
        nbasisk = getnbasis(basisk);
        m1 = m2 + 1;
        m2 = m2 + nbasisk;
        lambdak = getlambda(afdPark);
        if lambdak > 0
            Lfdobjk = getLfd(afdPark);
            penmatk = eval_penalty(basisk, Lfdobjk);
            bcoefk  = bvec(m1:m2);
            SSE     = SSE + lambdak.*bcoefk'*penmatk*bcoefk;
        end
    end
end

%  update derivatives of least squares SSE for terms for 
%  smoothing derivative weight coefficients

if derivs
    
    m2 = 0;
    for j=1:nderiv
        bfdParj = bwtcell{j};
        if getestimate(bfdParj) 
            basisj  = getbasis(getfd(bfdParj));
            nbasisj = getnbasis(basisj);
            m1 = m2 + 1;
            m2 = m2 + nbasisj;
            lambdaj = getlambda(bfdParj);
            if lambdaj > 0
                Lfdobjj = getLfd(bfdParj);
                penmatj = eval_penalty(basisj, Lfdobjj);
                bcoefj  = bvec(m1:m2);
                termj   = lambdaj.*penmatj*bcoefj;
                DSSE(m1:m2) = DSSE(m1:m2) + 2*termj;
            end
        end
    end
    
    %  update penalized least squares for terms for 
    %  smoothing forcing function coefficients
    
    nforce = length(ufdcell);
    for k=1:nforce
        afdPark = awtcell{k};
        if getestimate(afdPark) 
            basisk  = getbasis(getfd(afdPark));
            nbasisk = getnbasis(basisk);
            m1 = m2 + 1;
            m2 = m2 + nbasisk;
            lambdak = getlambda(afdPark);
            if lambdak > 0
                Lfdobjk = getLfd(afdPark);
                penmatk = eval_penalty(basisk, Lfdobjk);
                bcoefk  = bvec(m1:m2);
                DSSE(m1:m2) = DSSE(m1:m2) + 2*lambdak.*penmatk*bcoefk;
            end
        end
    end
    
end

%  ----------------------------------------------------------------

function [penmat, penvec, DR, Ds] = ...
               eval_Rs(bwtcell, awtcell, ufdcell, basisobj, derivs)
% EVAL_RS computes the penalty matrix R(bvec) 
%  and penalty vector s(bvec) where BVEC is a vector of parameters
%  defining the weight coefficients in WFDCELL and AFDCELL
%  If DERVIVS is positive, derivatives of these vectors are also
%  computed

%  Last modified 26 May 2004

if nargin < 5
    derivs = 1;
end

nderiv  = length(bwtcell);
nbasis  = getnbasis(basisobj);
basisfd = fd(eye(nbasis),basisobj);
onebas  = ones(1,nbasis);

%  set up the homogeneous operator

Lfdobj = Lfd(nderiv, bwtcell);

%  compute the penalty matrix

%  retrieve quadrature points and weights
quadvals = getquadvals(basisobj);
tvalquad = quadvals(:,1);
quadwts  = quadvals(:,2);
%  compute L weighted basis values
Lbasismat = getvalues(basisobj,nderiv);
for j=1:nderiv
    Dbasismat  = getvalues(basisobj,j-1);
    bfdj       = getfd(bwtcell{j});
    bbasisj    = getbasis(bfdj);
    bcoefj     = getcoef(bfdj);
    bbasismatj = getvalues(bbasisj);
    bvecj      = bbasismatj*bcoefj;
    Lbasismat  = Lbasismat + (bvecj*onebas).*Dbasismat;
end
%  penalty matrix R is cross-product
penmat = Lbasismat'*Lbasismat;

% penmat2 = eval_penalty(basisobj, Lfdobj);

%  evaluate the penalty vector if required

forcewrd = ~isempty(awtcell) & ~isempty(ufdcell);
if forcewrd
    %  set up the part of the roughness penalty affected by the
    %  presence of forcing function(s)
    nforce = length(ufdcell);
    penvec = zeros(nbasis,1);
    for k=1:nforce
        afdk       = getfd(awtcell{k});
        abasisk    = getbasis(afdk);
        acoefk     = getcoef(afdk);
        abasismatk = getvalues(abasisk);
        aveck      = abasismatk*acoefk;
        ufdk       = ufdcell{k};
        uveck      = eval_fd(tvalquad, ufdk);
        auveck     = aveck.*uveck;
        penveck    = auveck'*Lbasismat;
%  The following code is for numerical integration using INPROD
%         acoefk  = getcoef(afdk);
%         ucoefk  = getcoef(ufdk);
%         abasisk = getbasis(afdk);
%         ubasisk = getbasis(ufdk);
%         aucellk{1} = fd(acoefk,abasisk);
%         aucellk{2} = fd(ucoefk,ubasisk);
%         penveck = inprod(aucellk, basisobj, 0, Lfdobj);
        penvec  = penvec - penveck'; 
        aveccell{k} = aveck;
        uveccell{k} = uveck;
        if derivs
            auveccell{k} = auveck;
%             ffdcell{k} = aucellk;
        end
    end
else
    nforce = 0;
    penvec = [];
    Ds     = [];
end

if derivs
    
    %  evaluate the derivatives

    m2 = 0;
    fdbasis = fd(eye(nbasis),basisobj);
    for j=1:nderiv
        Lfdobjj = int2Lfd(j-1);
        bfdParj = bwtcell{j};
        if getestimate(bfdParj) 
            basisj     = getbasis(getfd(bfdParj));
            bbasismatj = getvalues(basisj);
            nbasisj    = getnbasis(basisj);
            m1 = m2 + 1;
            m2 = m2 + nbasisj;
            for m=m1:m2
                bfdvecm      = bbasismatj(:,m-m1+1);
                Djbasismatm  = getvalues(basisobj, j-1);
                bDjbasismatm = Djbasismatm.*(bfdvecm*onebas);
                Dpenmatm     = bDjbasismatm'*Lbasismat;
%  The following code is for numerical integration using INPROD
%                 bDbasiscell{1} = fd(coefm, basisj);
%                 bDbasiscell{2} = basisfd;
%                 Dpenmatm = inprod(bDbasiscell, basisfd, Lfdobjj, Lfdobj);
                DR(1:nbasis,1:nbasis,m) = full(Dpenmatm + Dpenmatm');
                if forcewrd
                    %  set up the part of the roughness penalty a
                    %  ffected by the presence of forcing function(s)
                    Dpenvecm = zeros(nbasis,1);
                    for k=1:nforce
                        auveck    = auveccell{k};
                        Dpenvecmk = bDjbasismatm'*auveck;
%  The following code is for numerical integration using INPROD
%                         aucellk   = ffdcell{k};
%                         Dpenvecmk = inprod(bDbasiscell, aucellk, Lfdobjj, 0);
                        Dpenvecm = Dpenvecm - Dpenvecmk;
                    end
                    Ds(1:nbasis,m) = Dpenvecm;
                end
            end
        end
    end
    
    %  update penalized least squares for terms for 
    %  smoothing forcing function coefficients
    
    for k=1:nforce
        afdPark = awtcell{k};
        if getestimate(afdPark) 
            abasisk    = getbasis(getfd(afdPark));
            abasismatk = getvalues(abasisk);
            nbasisk    = getnbasis(abasisk);
            uveck      = uveccell{k};
            m1 = m2 + 1;
            m2 = m2 + nbasisk;
            for m=m1:m2
                Daveck     = abasismatk(:,m-m1+1);
                Dauveck    = Daveck.*uveck;
                Dpenvecm   = -Dauveck'*Lbasismat;
%  The following code is for numerical integration using INPROD
%                 acoefk = zeros(nbasisk,1);
%                 acoefk(m-m1+1) = 1;
%                 Daucellk{1} = fd(acoefk, abasisk);
%                 Daucellk{2} = fd(ucoefk, ubasisk);
%                 Dpenvecm = -inprod(Daucellk, basisfd, 0, Lfdobj);
                Ds(1:nbasis,m) = Dpenvecm';
            end
        end
    end   
else
    DR = [];
    Ds = [];
end
