function [PENSSE, DSSE, PEN, coefcell] = ...
                      profPDA_MIMO(bvec, fitcell, gradwrd)
                      
% profPDAm estimates a system of M homogeneous linear 
%  differential equations directly from discrete noisy 
%  observations of a process.  
%  In this version forcing functions can be accommodated. 
%
%profPDAm works with the basis function expansions of the
%  estimates of the coefficient functions a_k(t) and b_j(t) 
%  in the possibly nonhomogeneous linear differential operator
%
%   L_i x_i(t) =  
%     b_{i01}(t) x_1(t) + ... + b_{i,m-1,1}(t) D^{m-1}x_1(t) +
%                          .
%                          .
%                          .
%     b_{iM1}(t) x_M(t) + ... + b_{i,m-1,M}(t) D^{m-1}x_M(t) + 
%     D^m x_i(t) - a_{1i}(t)u_1(t) - ... - a_{k_i,i}(t)u_K(t) 
%
%  that minimizes the least squares fitting criterion
%
%   \sum\i \{ \sigma_i^{-2} \| [y_i - x_i(t)] \|^2 + 
%             \lambda_i     \| L_i x_i    \|^2 \}
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
%    DDORDER  ... Order of the differential equation.
%    BWTCELL  ... Cell array of dimensions M and m for the weight functions  
%                   for the homogeneous part of the equation.  
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

%  Last modified 4 December 2004

if nargin < 3,  gradwrd = 0;  end

%  check number of variables

nvar = length(fitcell);

nbasisvec = zeros(nvar,1);

%  -------------------------------------------------------------------
%       check contents of each cell for cell array FITCELL.
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
        
    %  check BASISOBJ
    
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
    
    %  check DORDER
    
    Dorderi = fitstruct.Dorder;
    if ivar == 1
        Dorder = Dorderi;
        if Dorder < 1
            error('Order of differential equation less than one.');
        end
    else
        if Dorder ~= Dorderi
            error('DORDER is not the same for all variables.');
        end
    end
        
    %  check BWTCELL
    
    bwtcell = fitstruct.bwtcell;  
    
    if size(bwtcell,1) ~= nvar
        error(['First dimension of BWTCELL is not of length NVAR for FITCELL{', ...
                num2str(ivar), '}.']);
    end
    if size(bwtcell,2) ~= Dorder
        error('Second dimension of BWTCELL is not equal to DORDER.');
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
%       transfer parameters from vector BVEC to cell arrays
%                        BWTCELL and AWTCELL
%  -------------------------------------------------------------------

[fitcell, npar] = bvec2fitcell(bvec, fitcell);

%  -------------------------------------------------------------------
%         Compute the penalty matrices and penalty vectors
%  -------------------------------------------------------------------

[penmatcell, Dpenmatcell] = eval_Rsm(npar, fitcell, gradwrd);

%  -------------------------------------------------------------------
%          Define left and right sides of linear equation defining 
%                  the coefficient vector
%  -------------------------------------------------------------------

%  set up matrices for the linear equation

Cmat = zeros(ncoefs,ncoefs);
Pmat = zeros(ncoefs,ncoefs);
Dmat = zeros(ncoefs,1);

%  fill the matrices

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
    
    if isempty(fitstruct.awtcell) | ...
       isempty(fitstruct.ufdcell)
        nforce = 0;
    else
        nforce = length(awtcell);
    end
    
    %  set up diagonal block with coefficient matrix BMAT 
    %  plus lambda*R
    
    Cmat(indi,indi) = Bmati;
    Pmat(indi,indi) = Pmat(indi,indi) + lambdai.*Rmati;
    
    %  update all blocks of Cmat with the Smat matices
    
    Pmat = Pmat + lambdai.*Smati;
    
    %  update off-diagonal blocks of coefficient matrix 
    %  for rows or columns in block IVAR
    
    lamTmati = lambdai.*Tmati;
    Pmat(indi,:) = Pmat(indi,:) + lamTmati;
    Pmat(:,indi) = Pmat(:,indi) + lamTmati';
    
    %  set up right side vector DMAT
    
    Dmat(indi) = Dmat(indi) + Dmati;
    if nforce > 0
        Dmat(indi) = Dmat(indi) + lambdai.*Umati; 
        Dmat       = Dmat       + lambdai.*Vmati;
    end
    
end

%  -------------------------------------------------------------------
%                    Solve the linear system
%  -------------------------------------------------------------------

%  compute inverse of Cmat

Cmatinv = inv(Cmat + Pmat);

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

%  compute penalty term

PEN = coef' * Pmat * coef;

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
        for jderiv=1:Dorder
            bfdParj = bwtcell{jvar,jderiv};
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

if gradwrd
    
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
            mj1  = mj2 + 1;
            mj2  = mj2 + nbasisvec(jvar);
            indj = mj1:mj2;
            for jderiv=1:Dorder
                bfdParj = bwtcell{jvar,jderiv};
                if getestimate(bfdParj) 
                    bfdj    = getfd(bfdParj);
                    bbasisj = getbasis(bfdj);
                    bnbasis = getnbasis(bbasisj);
                    m1 = m2 + 1;
                    m2 = m2 + bnbasis;
                    for m = m1:m2
                        nbpar = nbpar + 1;
                        %  set up derivative supermatrix for this  parameter
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
                    DDmatm(indi) =          lambdai.*DUmat(:,napar);
                    DDmatm       = DDmatm + lambdai.*DaVmat(:,napar);
                    DAmatm       = PhiCmatinv*DDmatm;
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
        for jderiv=1:Dorder
            bfdParj = bwtcell{jvar};
            if getestimate(bfdParj)
                lambdaj = getlambda(bfdParj);
                bfdj = getfd(bfdParj);
                bbasisj = getbasis(bfdj);
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

%  -------------------------------------------------------------------

function [penmatcell, Dpenmatcell] = ...
              eval_Rsm(npar, fitcell, gradwrd)
% EVAL_RS computes the matrices and vectors required to define the
%  penalty for an order one system of equations with forcing functions.
%  these are stored in cell arrays PENMATCELL and DPENMATCELL.
%  NCOEFS is the total number of coefficients defining the fits
%  for the M variables, that is NBASIS1 + ... + NBASISM.
%  If DERIVS is positive, derivatives of these vectors are also
%  computed

%  Last modified 13 August 2004

if nargin < 4
    gradwrd = 1;
end

nvar   = length(fitcell);
nbasisvec = zeros(nvar,1);
for ivar=1:nvar
    fitstruct = fitcell{ivar};
    basisobji = fitstruct.basisobj;
    nbasis    = getnbasis(basisobji);
    nbasisvec(ivar) = nbasis;
end
ncoefs = sum(nbasisvec);

%  loop through variables

mi2  = 0;
for ivar = 1:nvar
    fitstruct = fitcell{ivar};
    basisobji = fitstruct.basisobj;
    nbasis    = getnbasis(basisobji);
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
    
    quadvals = getquadvals(basisobji);
    tvalquad = quadvals(:,1);
    
    %  retrieve weighted basis function values
    
    basismati = getvalues(basisobji);
    if ivar == 1
        Dorder = fitstruct.Dorder;
    end

    %  Compute penalty matrix Rmat, penalty vector Umat, and
    %  penalty scalar Wmat that depend only on variable IVAR
    %  compute L_phi  values for this variable
    
    Dmbasismati = getvalues(basisobji, Dorder);
    Rmat = Dmbasismati'*Dmbasismati;
    
    for k=1:nforce
        auveck = eval_au(quadvals, awtcell{k}, ufdcell{k});
        Umat   = Umat + Dmbasismati'*auveck;
        Wmat   = Wmat + sum(auveck.^2);
    end
    
    %  loop through other variables in this equation to compute
    %  penalty matrices Smat and Tmat
    
    mj2 = 0;
    for jvar=1:nvar
        fitstructj = fitcell{jvar};
        basisobjj  = fitstructj.basisobj;
        nbasisj    = getnbasis(basisobjj);
        mj1  = mj2 + 1;
        mj2  = mj2 + nbasisj;
        indj = mj1:mj2;
        
        %  compute L_phi  values for this variable
    
        tempj = eval_Lphi(jvar, bwtcell, basisobjj);
       
        %  Compute penalty matrix Tmat entries for this J
        
        Tmat(:,indj) = Dmbasismati'*tempj;
        
        %  loop through all over variables in this equation 
        %  to compute penalty matrix Smat entries for fixed J

        mk2 = 0;
        for kvar=1:jvar
            fitstructk = fitcell{kvar};
            basisobjk  = fitstructk.basisobj;
            nbasisk    = getnbasis(basisobjk);
            mk1  = mk2 + 1;
            mk2  = mk2 + nbasisk;
            indk = mk1:mk2;
            
            %  compute L_phi  values for this variable
            
            tempk = eval_Lphi(kvar, bwtcell, basisobjk);
            
            Smatkj = tempk'*tempj;
            Smat(indk,indj) = Smatkj;
            if kvar ~= jvar
                Smat(indj,indk) = Smatkj';
            end
        end
        
        %  Loop through forcing functions to compute 
        %  penalty matrix Vmat for this J that sums over
        %  forcing functions
        
        for k=1:nforce
            auveck = eval_au(quadvals, awtcell{k}, ufdcell{k});
            tempj = zeros(nbasisj,1);
            for jderiv = 1:Dorder
                bfdj       = getfd(bwtcell{jvar,jderiv});
                bbasisj    = getbasis(bfdj);
                bcoefj     = getcoef(bfdj);
                bbasismatj = getvalues(bbasisj);
                bvecj      = bbasismatj*bcoefj;
                Dbasismatj = getvalues(basisobjj,jderiv-1);
                tempj      = tempj + Dbasismatj'*(auveck.*bvecj);
            end
            Vmat(indj) = Vmat(indj) + tempj;
        end
    end
    
    %  place these matrices in struct PENSTRUCT
    
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

if gradwrd
    
    %  loop through equations
    
    m2  = 0;
    mi2 = 0;
    for ivar=1:nvar
        fitstruct = fitcell{ivar};
        penstruct = penmatcell{ivar};
        basisobj  = fitstruct.basisobj;
        nbasis = nbasisvec(ivar);
        onebas = ones(1,nbasis);
        Dmbasismat = getvalues(basisobj, Dorder);
        
        %  set up the weight and forcing function cells
        
        bwtcell = fitstruct.bwtcell;
        awtcell = fitstruct.awtcell;
        
        %  retrieve quadrature points and weights
        
        quadvals = getquadvals(basisobj);
        tvalquad = quadvals(:,1);
        
        %  set up arrays of derivative values
        
        nbcoefs = sum(npar(ivar,1:Dorder));
        nacoefs = npar(ivar,Dorder+1);
        DSmat   = zeros(ncoefs,ncoefs,nbcoefs);
        DTmat   = zeros(nbasis,ncoefs,nbcoefs);
        DUmat   = zeros(nbasis,       nacoefs);
        DaVmat  = zeros(ncoefs,       nacoefs);
        DbVmat  = zeros(ncoefs,       nbcoefs);
        nforce  = length(awtcell);
        
        mi1 = mi2 + 1;
        mi2 = mi2 + nbasis;
        
        %  loop through derivatives to be calculated within this equation
        
        mj2   = 0;
        nbpar = 0;
        for jvar=1:nvar
            fitstructj = fitcell{jvar};
            basisobjj  = fitstructj.basisobj;
            onebasj    = ones(1,nbasisvec(jvar));
            mj1        = mj2 + 1;
            mj2        = mj2 + nbasisvec(jvar);
            indj       = mj1:mj2;
            for jderiv=1:Dorder
                bfdParj = bwtcell{jvar,jderiv};
                if getestimate(bfdParj) 
                    bfdj       = getfd(bfdParj);
                    bbasisj    = getbasis(bfdj);
                    bbasismatj = getvalues(bbasisj);
                    nbbasisj   = getnbasis(bbasisj);
                    Dbasismatj = getvalues(basisobjj,jderiv-1);
                    
                    %  Compute DSmat by looping through all variables and
                    %    derivatives for fixed variable J
                    
                    m1  = m2 + 1;
                    m2  = m2 + nbbasisj;
                    mk2 = 0;
                    for kvar=1:nvar
                        fitstructk = fitcell{kvar};
                        basisobjk  = fitstructk.basisobj;
                        basismatk  = getvalues(basisobjk);
                        onebask    = ones(1,nbasisvec(kvar));
                        mk1   = mk2 + 1;
                        mk2   = mk2 + nbasisvec(kvar);
                        indk  = mk1:mk2;
                        tempk = eval_Lphi(kvar, bwtcell, basisobjk);
                        for m=m1:m2
                            indm = nbpar+m-m1+1;
                            bfdvecm      = bbasismatj(:,m-m1+1);
                            bDjbasismatm = Dbasismatj.*(bfdvecm*onebasj);
                            prodjkm      = bDjbasismatm'*tempk;
                            DSmat(indj,indk,indm) = ...
                                DSmat(indj,indk,indm) + prodjkm;
                            DSmat(indk,indj,indm) = ...
                                DSmat(indk,indj,indm) + prodjkm';
                        end
                    end
                    
                    %  compute DTmat and DbVmat for variable J
                    
                    for m=m1:m2
                        indm = nbpar+m-m1+1;
                        bfdvecm      = bbasismatj(:,m-m1+1);
                        bDjbasismatm = Dbasismatj.*(bfdvecm*onebasj);
                        prodjm       = Dmbasismat'*bDjbasismatm;
                        DTmat(:,indj,indm) = prodjm;
                        %  set up the part of the roughness penalty 
                        %  affected by the presence of forcing function(s)
                        Dpenvecm = zeros(nbasisvec(jvar),1);
                        if nforce > 0
                            for k=1:nforce
                                auveck    = eval_au(quadvals, awtcell{k}, ufdcell{k});
                                Dpenvecmk = bDjbasismatm'*auveck;
                                Dpenvecm  = Dpenvecm + Dpenvecmk;
                            end
                            DbVmat(indj,indm) = Dpenvecm;
                        end
                    end
                    nbpar = nbpar + nbbasisj;
                end
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
                    Dpenvecm = Dmbasismat'*Dauveck;
                    DUmat(:,napar) = Dpenvecm;
                    %  compute DaVmat by looping through variables
                    %  for fixed forcing function K
                    mk2 = 0;
                    for kvar=1:nvar
                        fitstructk = fitcell{kvar};
                        basisobjk  = fitstructk.basisobj;
                        mk1   = mk2 + 1;
                        mk2   = mk2 + nbasisvec(kvar);
                        indk  = mk1:mk2;
                        tempk = eval_Lphi(kvar, bwtcell, basisobjk);
                        tempk = tempk'*Dauveck;
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

%  ----------------------------------------------------------------

function temp = eval_Lphi(jvar, bwtcell, basisobj)
%  Evaluates the linear combination of basis function derivatives 
%  for variable JVAR in differential operator L
%  defined by weight cell array BWTCELL.
%  The evaluation is at the quadrature points.

%  Last modified 10 August 2004

Dorder    = size(bwtcell, 2);
nbasis    = getnbasis(basisobj);
onebas    = ones(1,nbasis);

%  the function term in the operator

bfd       = getfd(bwtcell{jvar,1});
bbasis    = getbasis(bfd);
bcoef     = getcoef(bfd);
bbasismat = getvalues(bbasis);
bmat      = (bbasismat*bcoef)*onebas;
basismat  = getvalues(basisobj);
temp      = basismat.*bmat;

%  the terms for the derivatives

for jderiv = 2:Dorder
    Dbfd       = getfd(bwtcell{jvar, jderiv});
    Dbbasis    = getbasis(Dbfd);
    Dbcoef     = getcoef(Dbfd);
    Dbbasismat = getvalues(Dbbasis);
    Dbmat      = (Dbbasismat*Dbcoef)*onebas;
    Dbasismat  = getvalues(basisobj,jderiv-1);
    temp       = temp + Dbasismat.*Dbmat;
end

%  --------------------------------------------------------------

function auvec = eval_au(quadvals, awtcell, ufd)
%  evaluate weight function times forcing function
afd       = getfd(awtcell);
abasis    = getbasis(afd);
acoef     = getcoef(afd);
abasismat = getvalues(abasis);
avec      = abasismat*acoef;
uvec      = eval_fd(quadvals(:,1), ufd);
auvec     = avec.*uvec.*sqrt(quadvals(:,2));

