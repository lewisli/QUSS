function [PENSSE, DSSE, coefcell, penmatcell, Dpenmatcell] = ...
                          profPDAm(bvec, fitcell, derivs)
                      
% profPDAm estimates a system of M homogeneous linear 
%  first order differential equations
%  directly from discrete noisy observations of a process.  
%  In this version forcing functions can be accommodated. 
%
%profPDAm works with the basis function expansions of the
%  estimates of the coefficient functions a_k(t) and b_j(t) 
%  in the possibly nonhomogeneous linear differential operator
%
%    L_ix_i(t) =  b_i1(t) x_1(t) + ... + b_iM x_m + D x_i(t) 
%                 - a_1(t)u_1(t) - ... - a_k(t)u_K(t)
%  that minimizes in a least squares sense the residual
%  functions f(t) = Lx(t).  
%
%  Arguments:
%  BVEC       ... A super-vector containing vectors of coefficients defining 
%                 the weight functions to be estimated.  
%                 Coefficients are stacked on top of each other as follows:
%                 coefficients for b_0 first, b_1 next, and so on, and
%                 then continues on to the forcing function coefficients
%                 to be estimated.  See BVEC2LFD for more details.
%  FITCELL    ... Cell array with one cell per variable.
%                 The contents of each cell are a struct containing:
%    Y        ... N by NCURVES matrix of curve values to be fit for this
%                 variable.  N and NCURVES can vary from one variable to
%                 another.
%    BASISOBJ ... A basis object for representing the curve(s)
%                   that are the solutions of the differential equations
%                   that fit the data
%    BASISMAT ... Matrix of basis function values at sampling points.
%    BMAT     ... Weighted cross-product matrix for basis matrix values.
%    DMAT     ... Weighted product of transposed basis matrix and Y.
%    LAMBDA   ... The smoothing parameter controlling the penalty on
%                   the roughness penalty.  In order to estimate the 
%                   DIFE, this should be large but not too large.  
%                   If it is too small, the roughness is ignored, and
%                   the differential equation is badly estimated.
%                   If it is too large, though, the data are ignored and
%                   the DIFE is also badly estimated.  
%    BWTCELL  ... Cell array for the weight functions for the 
%                   homogeneous part of the equation.
%    AWTCELL  ... Cell array for the weight functions for the 
%                   forcing functions
%    UFDCELL  ... Cell array containing functional data objects for
%                   the forcing functions
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

if nargin < 3
    derivs = 0;
end

%  check number of variables

nvar = length(fitcell);
nbasisvec = zeros(nvar,1);

%  -------------------------------------------------------------------
%                   check contents of each cell.
%  -------------------------------------------------------------------

N = 0;
for ivar = 1:nvar
    
    fitstruct = fitcell{ivar};

    %  check Y

    [n, ncurves] = size(fitstruct.y);
    if n == 1 & ncurves > 1
        %  transpose Y is n == 1
        y = y';
        fitstruct.y = y;
        fitcell{ivar} = fitstruct;
        temp    = n;
        n       = ncurves;
        ncurves = temp;
    end
    N = N + n;
        
    %  check basis
    
    basisobj  = fitstruct.basisobj;
    
    if ~isa_basis(fitstruct.basisobj)
        error(['BASISOBJ is not a basis object for FITCELL{', ...
                num2str(ivar), '}.']);
    end
    
    nbasis = getnbasis(basisobj);
    nbasisvec(ivar) = nbasis;
    
    %  check LAMBDA
    
    lambda = fitstruct.lambda;
    if lambda < 0
        warning ('Value of LAMBDA was negative, and 0 used instead.');
        fitstruct.lambda = 0;
        fitcell{ivar} = fitstruct;
    end
    
    %  check BWTCELL
    
    bwtcell = fitstruct.bwtcell;   
    if length(bwtcell) ~= nvar
        error(['BWTCELL is not of length NVAR for FITCELL{', ...
                num2str(ivar), '}.']);
    end
    
    %  check AWTCELL and UFDCELL
    
    awtcell = fitstruct.awtcell;   
    if isempty(awtcell) & isempty(fitstruct.ufdcell)
        nforce = 0;
    else
        nforce = length(awtcell);
        if length(fitstruct.ufdcell) ~= nforce
            error(['AWTCELL and UFDCELL are not of same length for FITCELL{', ...
                num2str(ivar), '}.']);
        end
    end
    
    %  check Bmat and Dmat
    
    Bmat    = fitstruct.Bmat;
    if size(Bmat,1) ~= nbasis | size(Bmat,2) ~= nbasis
        error(['BMAT has incorrect dimension(s) for FITCELL{', ...
                num2str(ivar), '}.']);
    end
    
    Dmat    = fitstruct.Dmat;
    if size(Dmat,1) ~= nbasis
        error(['DMAT has incorrect first dimension for FITCELL{', ...
                num2str(ivar), '}.']);
    end
    
end

ncoefs = sum(nbasisvec);
    
%  -------------------------------------------------------------------
%       transfer parameters from vector BVEC to cell objects
%                        BWTCELL and AWTCELL
%  -------------------------------------------------------------------
    
npar = zeros(nvar,2);  
m2 = 0;
for ivar = 1:nvar
    
    fitstruct = fitcell{ivar};
    bwtcell   = fitstruct.bwtcell;   
    awtcell   = fitstruct.awtcell;   
    if isempty(awtcell) | isempty(fitstruct.ufdcell)
        nforce = 0;
    else
        nforce = length(awtcell);
    end

    for jvar=1:nvar
        bfdParj = bwtcell{jvar};
        if getestimate(bfdParj)
            bfdobjj = getfd(bfdParj);
            bbasisj = getbasis(bfdobjj);
            bnbasis = getnbasis(bbasisj);
            m1 = m2 + 1;
            m2 = m2 + bnbasis;
            npar(ivar,1) = npar(ivar,1) + bnbasis;
            bfdobjj = putcoef(bfdobjj, bvec(m1:m2));
        end
        bwtcell{jvar} = putfd(bfdParj, bfdobjj);
    end
    fitstruct.bwtcell = bwtcell;
    
    %  loop through forcing functions to transfer coefficients from
    %  BVEC to the appropriate forcing function weight functions
    
    for k=1:nforce
        afdPark = awtcell{k};
        afdobjk = getfd(afdPark);
        if getestimate(afdPark)
            abasisk = getbasis(afdobjk);
            anbasis = getnbasis(abasisk);
            m1 = m2 + 1;
            m2 = m2 + anbasis;
            npar(ivar,2) = npar(ivar,2) + anbasis;
            afdobjk = putcoef(afdobjk, bvec(m1:m2));
        end
        awtcell{k} = putfd(afdPark, afdobjk);
    end
    fitstruct.awtcell = awtcell;
    fitcell{ivar} = fitstruct;
    
end

%  -------------------------------------------------------------------
%         Compute the penalty matrices and penalty vectors
%  -------------------------------------------------------------------

[penmatcell, Dpenmatcell] = eval_Rsm(nbasisvec, npar, fitcell, derivs);

%  -------------------------------------------------------------------
%          Define left and right sides of linear equation defining 
%                  the coefficient vector
%  -------------------------------------------------------------------

Cmat = zeros(ncoefs,ncoefs);
Dmat = zeros(ncoefs,1);
mi2 = 0;
for ivar=1:nvar
    mi1       = mi2 + 1;
    mi2       = mi2 + nbasisvec(ivar);
    indi      = mi1:mi2;
    
    %  extract data matrices
    
    fitstruct = fitcell{ivar};
    Bmati     = fitstruct.Bmat;
    Dmati     = fitstruct.Dmat;
    lambdai   = fitstruct.lambda;
    
    %  extract penalty matrices
    
    penstruct = penmatcell{ivar};
    Rmati     = penstruct.Rmat;
    Smati     = penstruct.Smat;
    Tmati     = penstruct.Tmat;
    Umati     = penstruct.Umat;
    Vmati     = penstruct.Vmat;
    
    if isempty(fitstruct.awtcell) | isempty(fitstruct.ufdcell)
        nforce = 0;
    else
        nforce = length(awtcell);
    end
    
    %  set up diagonal block with coefficient matrix BMAT 
    %  plus lambda*R
    
    Cmat(indi,indi) = Bmati + lambdai.*Rmati;
    
    %  update all blocks of Cmat with the Smat matices
    
    Cmat = Cmat + lambdai.*Smati;
    
    %  update off-diagonal blocks of coefficient matrix 
    %  for rows or columns in block IVAR
    
    Cmat(indi,:) = Cmat(indi,:) + lambdai.*Tmati;
    Cmat(:,indi) = Cmat(:,indi) + lambdai.*Tmati';
    
    %  set up right side vector DMAT
    
    Dmat(indi) = Dmat(indi) + Dmati;
    if nforce > 0
        Dmat(indi) = Dmat(indi) + lambdai.*Umati; 
        Dmat       = Dmat + lambdai.*Vmati;
    end
    
end

%  -------------------------------------------------------------------
%                    Solve the linear system
%  -------------------------------------------------------------------

%  compute inverse of Cmat

Cmatinv = inv(Cmat);

%  solve normal equations for each observation

coef = Cmatinv * Dmat;

%  -------------------------------------------------------------------
%                 Compute error sum of squares
%  -------------------------------------------------------------------

SSE = 0;
m2  = 0;
res = [];
for ivar=1:nvar
    m1 = m2 + 1;
    m2 = m2 + nbasisvec(ivar);
    indi = m1:m2;
    coefi = coef(indi);
    coefcell{ivar} = coefi;
    
    %  compute error sum of squares
    
    fitstruct = fitcell{ivar};
    basismati = fitstruct.basismat;
    yi        = fitstruct.y;
    yhati     = basismati * coefi;
    resi      = yi - yhati;
    res       = [res; resi];
    SSEi      = sum(sum(resi.^2));
    SSE       = SSE + SSEi;
end

%  -------------------------------------------------------------------
%       Update the error sum of squares 
%     for any roughness penalties on parameters
%  -------------------------------------------------------------------

PENSSE = SSE;

for ivar = 1:nvar
    
    fitstruct = fitcell{ivar};    
    bwtcell   = fitstruct.bwtcell;   
    awtcell   = fitstruct.awtcell;   

    if isempty(awtcell) | isempty(fitstruct.ufdcell)
        nforce = 0;
    else
        nforce = length(awtcell);
    end

    for jvar=1:nvar
        bfdParj = bwtcell{jvar};
        if getestimate(bfdParj)
            lambdaj = getlambda(bfdParj);
            if lambdaj > 0
                bfdobjj = getfd(bfdParj);
                bcoefj  = getcoef(bfdParj);
                bbasisj = getbasis(bfdobjj);
                Lfdobjj = getLfd(bfdParj);
                bRmatj  = eval_penalty(bbasisj, Lfdobjj);
                PENSSE  = PENSSE + lambdaj.*bcoefj'*bRmatj*bcoefj;
            end
        end
    end
    
    %  loop through forcing functions to transfer coefficients from
    %  BVEC to the appropriate forcing function weight functions
    
    for k=1:nforce
        afdPark = awtcell{k};
        if getestimate(afdPark)
            lambdak  = getlambda(afdPark);
            if lambdak > 0
                afdobjk = getfd(afdPark);
                acoefk  = getcoef(afdPark);
                abasisk = getbasis(afdobjk);
                aLfdk   = getLfd(afdPark);
                aRmatk  = eval_penalty(abasisk, aLfdk);
                PENSSE  = PENSSE + lambdak.*acoefk'*aRmatk*acoefk;
            end
        end
    end
    
end

%  --------------------------------------------------------------
%                       Compute gradient
%  --------------------------------------------------------------

if derivs
    
    %  set up supermatrix containing submatrices of basis function
    %  values on the diagonal blocks
    
    Phimat = zeros(N,ncoefs);
    mi2 = 0;
    Ni2 = 0;
    for ivar=1:nvar
        mi1 = mi2 + 1;
        mi2 = mi2 + nbasisvec(ivar);
        fitstruct = fitcell{ivar};
        basismat  = fitstruct.basismat;
        Ni        = length(fitstruct.y);
        Ni1 = Ni2 + 1;
        Ni2 = Ni2 + Ni;
        Phimat(Ni1:Ni2,mi1:mi2) = basismat;
    end
    
    %  multiply by the inverse of the coefficient supermatrix
    
    PhiCmatinv = Phimat*Cmatinv;
    
    %  the outer two loops go through parameters to be estimated
    
    %  loop through the equations
    
    mi2  = 0;
    mpar = 0;
    m2   = 0;
    for ivar=1:nvar
        mi1  = mi2 + 1;
        mi2  = mi2 + nbasisvec(ivar);
        indi = mi1:mi2;
        fitstruct  = fitcell{ivar};
        awtcell    = fitstruct.awtcell;
        bwtcell    = fitstruct.bwtcell;
        lambdai    = fitstruct.lambda;
        %  retrieve the derivatives of penalty matrices
        Dpenstruct = Dpenmatcell{ivar}; 
        DSmat      = Dpenstruct.DSmat;
        DTmat      = Dpenstruct.DTmat;
        DUmat      = Dpenstruct.DUmat;
        DaVmat     = Dpenstruct.DaVmat;
        DbVmat     = Dpenstruct.DbVmat;
        
        %  loop through variables J inside this equation I to
        %  compute derivatives with respect to weight functions
        %  variables
        
        mj2  = 0;
        nbpar = 0;
        for jvar = 1:nvar
            bfdParj = bwtcell{jvar};
            mj1  = mj2 + 1;
            mj2  = mj2 + nbasisvec(jvar);
            indj = mj1:mj2;
            if getestimate(bfdParj) 
                bfdobjj = getfd(bfdParj);
                bbasisj = getbasis(bfdobjj);
                bnbasis = getnbasis(bbasisj);
                %  
                m1 = m2 + 1;
                m2 = m2 + bnbasis;
                for m = m1:m2
                    nbpar = nbpar + 1;
                    %  set up derivative supermatrix for this
                    %  parameter
                    DCmat  = zeros(ncoefs,ncoefs);
                    %  update the supermatrix for DSmat
                    DCmat = DCmat + lambdai.*DSmat(:,:,nbpar);
                    %  update the supermatrix for DTmat
                    temp = lambdai.*DTmat(:,:,nbpar);
                    DCmat(indi,:) = DCmat(indi,:) + temp;
                    DCmat(:,indi) = DCmat(:,indi) + temp';
                    %  compute the derivative of the smoothing matrix
                    if nforce > 0
                        %  forcing functions present:
                        DDmatm = zeros(ncoefs,1);
                        DDmatm(indj) = DDmatm(indj) + ...
                                lambdai.*DbVmat(indj,nbpar);
                        DAmatm = PhiCmatinv* ...
                                 (-DCmat*coef + DDmatm*ones(1,ncurves));
                    else
                        %  no forcing functions
                        DAmatm = -PhiCmatinv*DCmat*coef;
                    end
                    mpar  = mpar  + 1;
                    DSSE(mpar) = -2.*trace(res'*DAmatm);
                end
            end
        end
        
        %  loop through forcing functions to compute derivatives
        %  with respect to forcing function coefficients
        
        napar = 0;
        for k=1:nforce
            afdPark = awtcell{k};
            if getestimate(afdPark) 
                abasisk  = getbasis(getfd(afdPark));
                nbasisk = getnbasis(abasisk);
                m1 = m2 + 1;
                m2 = m2 + nbasisk;
                for m = m1:m2
                    napar = napar + 1;
                    %  compute derivative or right side for
                    %  this parameter
                    DDmatm = zeros(ncoefs,1);
                    DDmatm(indi) = lambdai.*DUmat(:,napar);
                    DDmatm = DDmatm + lambdai.*DaVmat(:,napar);
                    DAmatm  = PhiCmatinv*DDmatm;
                    mpar = mpar  + 1;
                    DSSE(mpar) = -2.*sum(res'*DAmatm);
                end
            end
        end
    end
%  -------------------------------------------------------------------
%       Update the error sum of squares 
%     for any roughness penalties on parameters
%  -------------------------------------------------------------------

DPENSSE = DSSE;

m2 = 0;
for ivar = 1:nvar
    
    fitstruct = fitcell{ivar};    
    bwtcell   = fitstruct.bwtcell;   
    awtcell   = fitstruct.awtcell;   

    if isempty(awtcell) | isempty(fitstruct.ufdcell)
        nforce = 0;
    else
        nforce = length(awtcell);
    end

    for jvar=1:nvar
        bfdParj = bwtcell{jvar};
        if getestimate(bfdParj)
            lambdaj = getlambda(bfdParj);
            bfdobjj = getfd(bfdParj);
            bbasisj = getbasis(bfdobjj);
            bnbasis = getnbasis(bbasisj);
            m1      = m2 + 1;
            m2      = m2 + bnbasis;
            indm    = m1:m2;
            if lambdaj > 0
                bcoefj = getcoef(bfdParj);
                bLfdj  = getLfd(bfdParj);
                bRmatj = eval_penalty(bbasisj, bLfdj);
                DPENSSE(indm) = DPENSSE(indm) + ...
                                2*lambdaj.*bRmatj*bcoefj;
            end
        end
    end
    
    %  loop through forcing functions to transfer coefficients from
    %  BVEC to the appropriate forcing function weight functions
    
    for k=1:nforce
        afdPark = awtcell{k};
        if getestimate(afdPark)
            lambdak = getlambda(afdPark);
            afdobjk = getfd(afdPark);
            abasisk = getbasis(afdobjk);
            anbasis = getnbasis(abasisk);
            m1      = m2 + 1;
            m2      = m2 + anbasis;
            indm    = m1:m2;
            if lambdak > 0
                acoefk = getcoef(afdPark);
                aLfdk  = getLfd(afdPark);
                aRmatk = eval_penalty(abasisk, aLfdk);
                DPENSSE(indm) = DPENSSE(indm) + ...
                                2*lambdak.*aRmatk*acoefk;
            end
        end
    end
    
end

else
    DSSE = [];
end

%  ----------------------------------------------------------------

function [penmatcell, Dpenmatcell] = ...
              eval_Rsm(nbasisvec, npar, fitcell, derivs)
% EVAL_RS computes the matrices and vectors required to define the
%  penalty for an order one system of equations with forcing functions.
%  these are stored in cell arrays PENMATCELL and DPENMATCELL.
%  NCOEFS is the total number of coefficients defining the fits
%  for the M variables, that is NBASIS1 + ... + NBASISM.
%  If DERIVS is positive, derivatives of these vectors are also
%  computed

%  Last modified 15 June 2004

if nargin < 4
    derivs = 1;
end

nvar   = length(fitcell);
ncoefs = sum(nbasisvec);

%  loop through variables

mi2  = 0;
for ivar = 1:nvar
    fitstruct = fitcell{ivar};
    basisobj  = fitstruct.basisobj;
    nbasis    = getnbasis(basisobj);
    onebasi   = ones(1,nbasis);
    mi1       = mi2 + 1;
    mi2       = mi2 + nbasis;
    indi      = mi1:mi2;
    
    %  set up the weight and forcing function cells
    
    bwtcell = fitstruct.bwtcell;
    awtcell = fitstruct.awtcell;
    ufdcell = fitstruct.ufdcell;
    if isempty(awtcell) | isempty(ufdcell)
        nforce = 0;
    else
        nforce = length(awtcell);
    end
    
    %  set up arrays to be filled
    
    Rmat = zeros(nbasis,nbasis);
    Smat = zeros(ncoefs,ncoefs);
    Tmat = zeros(nbasis,ncoefs);
    if nforce > 0
        Umat = zeros(nbasis,1);
        Vmat = zeros(ncoefs,1);
        Wmat = 0;
    else
        Umat = [];
        Vmat = [];
        Wmat = [];
    end
    
    %  retrieve quadrature points and weights
    
    quadvals = getquadvals(basisobj);
    tvalquad = quadvals(:,1);
    
    %  retrieve weighted basis function values
    
    basismati  = getvalues(basisobj);
    Dbasismati = getvalues(basisobj,1);
    
    %  compute weight function values for this variable
    
    bfdi       = getfd(bwtcell{ivar});
    bbasisi    = getbasis(bfdi);
    bcoefi     = getcoef(bfdi);
    bbasismati = getvalues(bbasisi);
    bveci      = bbasismati*bcoefi;
    
    tempi      = basismati.*(bveci*onebasi);
    
    %  Compute penalty matrix Rmat, penalty vector Umat, and
    %  penalty scalar Wmat that depend only on variable IVAR
    
    Rmat   = Dbasismati'*Dbasismati;
    for k=1:nforce
        afdk       = getfd(awtcell{k});
        abasisk    = getbasis(afdk);
        acoefk     = getcoef(afdk);
        abasismatk = getvalues(abasisk);
        aveck      = abasismatk*acoefk;
        ufdk       = ufdcell{k};
        uveck      = eval_fd(tvalquad, ufdk);
        auveck     = aveck.*uveck.*sqrt(quadvals(:,2));
        Umat       = Umat + Dbasismati'*auveck;
        Wmat       = Wmat + sum(auveck.^2);
    end
    
    %  loop through other variables in this equation to compute
    %  penalty matrices Smat and Tmat
    
    mj2 = 0;
    for jvar=1:nvar
        fitstructj = fitcell{jvar};
        basisobjj  = fitstructj.basisobj;
        nbasisj    = getnbasis(basisobjj);
        onebasj    = ones(1,nbasisj);
        mj1        = mj2 + 1;
        mj2        = mj2 + nbasisj;
        indj       = mj1:mj2;
        basismatj  = getvalues(basisobjj);
        bfdj       = getfd(bwtcell{jvar});
        bbasisj    = getbasis(bfdj);
        bcoefj     = getcoef(bfdj);
        bbasismatj = getvalues(bbasisj);
        bvecj      = bbasismatj*bcoefj;
        tempj      = basismatj.*(bvecj*onebasj);
        %  loop through all over variables in this equation 
        %  to compute penalty matrix Smat entries for fixed J
        mk2 = 0;
        for kvar=1:jvar
            fitstructk = fitcell{kvar};
            basisobjk  = fitstructk.basisobj;
            nbasisk    = getnbasis(basisobjk);
            mk1        = mk2 + 1;
            mk2        = mk2 + nbasisk;
            indk       = mk1:mk2;
            basismatk  = getvalues(basisobjk);
            bfdk       = getfd(bwtcell{kvar});
            bbasisk    = getbasis(bfdk);
            bcoefk     = getcoef(bfdk);
            bbasismatk = getvalues(bbasisk);
            bveck      = bbasismatk*bcoefk;
            onebask    = ones(1,nbasisk);
            tempk      = basismatk.*(bveck*onebask);
            Smat(indk,indj) = tempk'*tempj;
            Smat(indj,indk) = Smat(indk,indj)';
        end
        %  Compute penalty matrix Tmat entries for this J
        Tmat(:,indj) = Dbasismati'*tempj;
        %  Loop through forcing functions to compute 
        %  penalty matrix Vmat for this J that sums over
        %  forcing functions
        for k=1:nforce
            afdk       = getfd(awtcell{k});
            abasisk    = getbasis(afdk);
            acoefk     = getcoef(afdk);
            abasismatk = getvalues(abasisk);
            aveck      = abasismatk*acoefk;
            ufdk       = ufdcell{k};
            uveck      = eval_fd(tvalquad, ufdk);
            auveck     = aveck.*uveck.*sqrt(quadvals(:,2));
            Vmat(indj) = Vmat(indj) + basismatj'*(auveck.*bvecj);
        end
    end
    penstruct.Rmat = Rmat;
    penstruct.Smat = Smat;
    penstruct.Tmat = Tmat;
    penstruct.Umat = Umat;
    penstruct.Vmat = Vmat;
    penstruct.Wmat = Wmat;
    
    penmatcell{ivar} = penstruct;
    
end

%  ------------------------------------------------------------
%                 evaluate the derivatives
%  ------------------------------------------------------------

if derivs
    
    %  loop through equations
    
    m2  = 0;
    mi2 = 0;
    for ivar=1:nvar
        fitstruct  = fitcell{ivar};
        penstruct  = penmatcell{ivar};
        basisobj   = fitstruct.basisobj;

        %  set up the weight and forcing function cells
    
        bwtcell    = fitstruct.bwtcell;
        awtcell    = fitstruct.awtcell;
        
        %  retrieve quadrature points and weights
    
        quadvals = getquadvals(basisobj);
        tvalquad = quadvals(:,1);
    
        %  retrieve weighted basis function values
    
        D0basismat = getvalues(basisobj);
        D1basismat = getvalues(basisobj, 1);
        nbasis     = nbasisvec(ivar);
        onebas     = ones(1,nbasis);

        %  set up arrays to be filled
    
        DSmat      = zeros(ncoefs,ncoefs,npar(ivar,1));
        DTmat      = zeros(nbasis,ncoefs,npar(ivar,1));
        DUmat      = zeros(nbasis,npar(ivar,2));
        DaVmat     = zeros(ncoefs,npar(ivar,2));
        DbVmat     = zeros(ncoefs,npar(ivar,1));
        nforce     = length(awtcell);
        
        mi1 = mi2 + 1;
        mi2 = mi2 + nbasis;
        
        %  loop through variables in this equation
        
        mj2   = 0;
        nbpar = 0;
        for jvar=1:nvar
            bfdParj = bwtcell{jvar};
            if getestimate(bfdParj) 
                fitstructj = fitcell{jvar};
                basisobjj  = fitstructj.basisobj;
                basismatj  = getvalues(basisobjj);
                onebasj    = ones(1,nbasisvec(jvar));
                mj1        = mj2 + 1;
                mj2        = mj2 + nbasisvec(jvar);
                indj       = mj1:mj2;
                bfdj       = getfd(bwtcell{jvar});
                bbasisj    = getbasis(bfdj);
                bbasismatj = getvalues(bbasisj);
                nbbasisj   = getnbasis(bbasisj);
                
                %  Compute DSmat by looping through all variables
                %    for fixed variable J
                
                m1  = m2 + 1;
                m2  = m2 + nbbasisj;
                mk2 = 0;
                for kvar=1:nvar
                    fitstructk = fitcell{kvar};
                    basisobjk  = fitstructk.basisobj;
                    basismatk  = getvalues(basisobjk);
                    onebask    = ones(1,nbasisvec(kvar));
                    mk1        = mk2 + 1;
                    mk2        = mk2 + nbasisvec(kvar);
                    indk       = mk1:mk2;
                    bfdk       = getfd(bwtcell{kvar});
                    bbasisk    = getbasis(bfdk);
                    bcoefk     = getcoef(bfdk);
                    bbasismatk = getvalues(bbasisk);
                    bveck      = bbasismatk*bcoefk;
                    tempk      = basismatk.*(bveck*onebask);
                    for m=m1:m2
                        indm         = nbpar+m-m1+1;
                        bfdvecm      = bbasismatj(:,m-m1+1);
                        bDjbasismatm = basismatj.*(bfdvecm*onebasj);
                        prodjkm      = bDjbasismatm'*tempk;
                        DSmat(indj,indk,indm) = DSmat(indj,indk,indm) + ...
                                                prodjkm;
                        DSmat(indk,indj,indm) = DSmat(indk,indj,indm) + ...
                                                prodjkm';
                    end
                end
                
                %  compute DTmat and DbVmat for variable J
                
                for m=m1:m2
                    bfdvecm      = bbasismatj(:,m-m1+1);
                    bDjbasismatm = basismatj.*(bfdvecm*onebasj);
                    prod         = D1basismat'*bDjbasismatm;
                    DTmat(:,indj,nbpar+m-m1+1) = prod;
                    %  set up the part of the roughness penalty 
                    %  affected by the presence of forcing function(s)
                    Dpenvecm = zeros(nbasisvec(jvar),1);
                    if nforce > 0
                        for k=1:nforce
                            afdk       = getfd(awtcell{k});
                            abasisk    = getbasis(afdk);
                            acoefk     = getcoef(afdk);
                            abasismatk = getvalues(abasisk);
                            aveck      = abasismatk*acoefk;
                            ufdk       = ufdcell{k};
                            uveck      = eval_fd(tvalquad, ufdk);
                            auveck     = aveck.*uveck.*sqrt(quadvals(:,2));
                            Dpenvecmk  = bDjbasismatm'*auveck;
                            Dpenvecm   = Dpenvecm + Dpenvecmk;
                        end
                        DbVmat(indj,nbpar+m-m1+1) = Dpenvecm;
                    end
                end
                nbpar = nbpar + nbbasisj;
            end
        end
        
        %  update penalized least squares for terms for 
        %  smoothing forcing function coefficients
        
        napar = 0;
        for k=1:nforce
            afdPark = awtcell{k};
            if getestimate(afdPark) 
                abasisk    = getbasis(getfd(afdPark));
                abasismatk = getvalues(abasisk);
                anbasis    = getnbasis(abasisk);
                ufdk       = ufdcell{k};
                uveck      = eval_fd(tvalquad, ufdk);
                m1 = m2 + 1;
                m2 = m2 + anbasis;
                for m=m1:m2
                    napar    = napar + 1;
                    Daveck   = abasismatk(:,m-m1+1);
                    Dauveck  = Daveck.*uveck.*sqrt(quadvals(:,2));
                    Dpenvecm = D1basismat'*Dauveck;
                    DUmat(:,napar) = Dpenvecm;
                    %  compute DaVmat by looping through variables
                    %  for fixed forcing function K
                    mk2 = 0;
                    for kvar=1:nvar
                        fitstructk = fitcell{kvar};
                        basisobjk  = fitstructk.basisobj;
                        basismatk  = getvalues(basisobjk);
                        onebask    = ones(1,nbasisvec(kvar));
                        mk1        = mk2 + 1;
                        mk2        = mk2 + nbasisvec(kvar);
                        indk       = mk1:mk2;
                        bfdk       = getfd(bwtcell{kvar});
                        bbasisk    = getbasis(bfdk);
                        bcoefk     = getcoef(bfdk);
                        bbasismatk = getvalues(bbasisk);
                        bveck      = bbasismatk*bcoefk;
                        tempk      = basismatk.*(bveck*onebask);
                        tempk      = tempk'*Dauveck;
                        DaVmat(indk,napar) = tempk;
                    end
                end
            end
        end  
        Dpenstruct.DSmat  = DSmat;
        Dpenstruct.DTmat  = DTmat;
        Dpenstruct.DaVmat = DaVmat;
        Dpenstruct.DbVmat = DbVmat;
        Dpenstruct.DUmat  = DUmat;
        Dpenmatcell{ivar} = Dpenstruct;
        
    end
    
else
    Dpenmatcell = {};
end


