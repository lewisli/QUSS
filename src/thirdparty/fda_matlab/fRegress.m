function [betaestcell, yhatfdobj, betastderrcell, bvar, c2bmap] = ...
    fRegress(yfdPar, xfdcell, betacell, y2cmap, sigmae)
%  FREGRESS  Fits a functional linear model using multiple 
%  functional independent variables with the dependency being 
%  pointwise or concurrent.  
%  The case of a scalar independent variable is included by treating
%  it as a functional independent variable with a constant basis
%  and a unit coefficient.
%
%  Arguments:
%  YFDPAR   ... an object for the dependent variable, 
%               which may be: 
%                   a functional data object, 
%                   a functional parameter (fdPar) object, or 
%                   a vector
%  XFDCELL  ... a cell object of length p with each cell 
%               containing an object for an independent variable.
%               the object may be: 
%                   a functional data object or
%                   a vector
%  BETACELL ... a cell object of length p with each cell
%               containing a functional parameter object for 
%               the corresponding regression function.
%  Y2CMAP   ... the matrix mapping from the vector of observed values
%               to the coefficients for the dependent variable.  
%               This is output by function SMOOTH_BASIS.  If this is 
%               supplied, confidence limits are computed, otherwise not.
%  SIGMAE   ... Estimate of the covariances among the residuals.  This
%               can only be estimated after a preliminary analysis.
%               If this is not provided, it is automatically estimated.
%  
%  Returns:  
%  BETAESTCELL    ... a cell object, each cell containing a 
%                     functional parameter object
%                     for a regression function.
%  YHATFDOBJ      ... a functional data object for the predicted values
%                     of the response variable.
%  BETASTDERRCELL ... a cell object, each cell containing a fdPar object
%                     for the standard error of a regression function.
%  BVARIANCE      ... the symmetric matrix of sampling variances and
%                     covariances for the matrix of regression coefficients
%                     for the regression functions.  These are stored
%                     column-wise in defining BVARIANCE.  This is only 
%                     computed if the argument Y2CMAP is supplied.
%  C2BMAP         ... the matrix mapping from response variable coefficients 
%                     to coefficients for regression coefficients

%  Last modified 15 February 2005

if nargin < 2
    error('Less than two arguments supplied.');
end

%  Get sample size and check YFDPAR.  

if strcmp(class(yfdPar), 'fdpar') | ...
   strcmp(class(yfdPar), 'fd')
    %  ----------------------------------------------------------------
    %                   YFDPAR is functional
    %  ----------------------------------------------------------------
    
    if strcmp(class(yfdPar), 'fd')
        yfdPar = fdPar(yfdPar);
    end
    yfdobj  = getfd(yfdPar);
    ylambda = getlambda(yfdPar);
    yLfdobj = getLfd(yfdPar);
    ycoef   = getcoef(yfdobj);
    if length(size(ycoef)) > 2
        error('YFDOBJ from YFDPAR is not univariate.');
    end
    N         = size(ycoef,2);
    ybasisobj = getbasis(yfdobj);
    rangeval  = getbasisrange(ybasisobj);
    ynbasis   = getnbasis(ybasisobj);
    nfine     = max(501,10*ynbasis+1);
    tfine     = linspace(rangeval(1), rangeval(2), nfine)';
    ywtvec    = ones(nfine,1);
    deltat    = tfine(2) - tfine(1);
    ymat      = eval_fd(tfine, yfdobj);
    
    %  get number of independent variables 
    
    p = length(xfdcell);
    
    %  check each cell.  If the object is a vector of length N,
    %  it is converted to a functional data object with a 
    %  constant basis
    
    onebasis = create_constant_basis(rangeval);
    
    for j=1:p
        xfdj = xfdcell{j};
        if isa_fd(xfdj)
            xcoef = getcoef(xfdj);
            if length(size(xcoef)) > 2
                error(['Covariate ',num2str(j),' is not univariate.']);
            end
            rangevalx  = getbasisrange(getbasis(xfdj));
            if any(rangevalx ~= rangeval)
                error(['Range for covariate ',num2str(j), ...
                        ' does not match that of YFDOBJ.']);
            end
        elseif strcmp(class(xfdj), 'double')
            xfdcell{j} = fd(xfdj(:)', onebasis);
        else
            error(['Covariate ',num2str(j), ...
                    ' is neither a functional nor a multivariate object.']);
        end
        %   set up default betacell argument if required
        if nargin < 3
            betacell{j} = fdPar(onebasis);
        end
        %  check size of coefficient array
        coefj = getcoef(xfdcell{j});
        Nj = size(coefj, 2);
        if Nj ~= N
            error('Incorrect number of replications in XFDCELL');
        end
    end
    
    if length(betacell) ~= p
        error(['Number of regression coefficients does not match', ...
                ' number of independent variables.']);
    end
    
    %  set up a matrix of values of covariates over a fine mesh
    
    xmat = zeros(nfine, N, p);
    for j=1:p
        xmatj          = eval_fd(tfine, xfdcell{j});
        xmat(:,:,j)    = xmatj;
        betabasisj     = getbasis(getfd(betacell{j}));
        betamatcell{j} = eval_basis(tfine, betabasisj);
    end
    
    %  -----------------------------------------------------------
    %          set up the linear equations for the solution
    %  -----------------------------------------------------------
    
    %  compute the total number of coefficients to be estimated
    
    ncoef = 0;
    for j=1:p
        betafdParj = betacell{j};
        betafdj    = getfd(betafdParj);
        ncoefj     = size(getcoef(betafdj),1);
        ncoef      = ncoef + ncoefj;
    end
    
    Cmat = zeros(ncoef,ncoef);
    Dmat = zeros(ncoef,1);
    
    %  loop through rows of CMAT
    
    mj2 = 0;
    for j=1:p
        betafdParj = betacell{j};
        betafdj    = getfd(betafdParj);
        ncoefj     = size(getcoef(betafdj),1);
        betabasisj = getbasis(betafdj);
        xfdj       = xfdcell{j};
        if ~isa_fd(xfdj)
            if strcmp(class(xfdj), 'double')
                xfdj = fd(xfdj(:)', onebasis);
            else
                error(['Independent variable ', num2str(j), ...
                        ' is neither functional nor scalar.']);
            end
        end
        mj1 = mj2 + 1;
        mj2 = mj2 + ncoefj;
        indexj = mj1:mj2;
        %  compute right side of equation DMAT
        xywtvec  = sum(squeeze(xmat(:,:,j)).*ymat,2);
        betamatj = betamatcell{j};
        Dmatj    = deltat.*sum(betamatj'.*(xywtvec*ones(1,ncoefj))',2);
        Dmat(indexj) = Dmatj;
        %  loop through columns of CMAT
        mk2 = 0;
        for k=1:j
            betafdPark = betacell{k};
            betafdk    = getfd(betafdPark);
            ncoefk     = size(getcoef(betafdk),1);
            betabasisk = getbasis(betafdk);
            xfdk       = xfdcell{k};
            if ~isa_fd(xfdk)
                if strcmp(class(xfdk), 'double')
                    xfdk = fd(xfdk(:)', onebasis);
                else
                    error(['Independent variable ', num2str(k), ...
                            ' is neither functional nor scalar.']);
                end
            end
            mk1 = mk2 + 1;
            mk2 = mk2 + ncoefk;
            indexk = mk1:mk2;
            %  set up two weight functions
            xxwtvec  = sum(xmat(:,:,j).*xmat(:,:,k),2);
            betamatk = betamatcell{k};
            Cmatjk   = deltat.*((betamatj.*(xxwtvec*ones(1,ncoefj)))' ...
                              *betamatk);
            Cmat(indexj,indexk) = Cmatjk;
            Cmat(indexk,indexj) = Cmatjk';
        end
        %  attach penalty term to diagonal block
        lambda = getlambda(betafdParj);
        if lambda > 0
            Lfdj  = getLfd(betafdParj);
            Rmatj = eval_penalty(betafdParj);
            Cmat(indexj,indexj) = Cmat(indexj,indexj) + lambda.*Rmatj;
        end
    end
    
    %  check Cmat for singularity
    
    eigval = sort(eig(Cmat));
    if (eigval(1) < 0)
        disp('Smallest 10 eigenvalues:')
        for ieig=1:10
            fprintf('  %g \n',eigval(ieig));
        end
        disp('Largest  10 eigenvalues:')
        for ieig=ncoef-9:ncoef
            fprintf('  %g \n',eigval(ieig));
        end
        error('Negative eigenvalue of coefficient matrix.');
    end
    if (eigval(1) == 0)
        error('Zero eigenvalue of coefficient matrix.');
    end
    logcondition = log10(eigval(ncoef)) - log10(eigval(1));
    if logcondition > 12
        warning('Near singularity in coefficient matrix.');
        warning(['Eigenvalues range from ',num2str(eigval(1)), ...
                 ' to ',num2str(eigval(ncoef)),'.']);
    end
    
%     fprintf(['Log10 condition number = ',...
%               num2str(logcondition),'\n']);
          
    %  solve for coefficients defining BETA
    
    Cmatinv  = inv(Cmat);
    betacoef = Cmatinv*Dmat;
    
    %  set up fdPar object for BETAFDPAR
    
    betaestcell = betacell;
    mj2 = 0;
    for j=1:p
        betafdParj     = betacell{j};
        betafdj        = getfd(betafdParj);
        ncoefj = size(getcoef(betafdj),1);
        mj1    = mj2 + 1;
        mj2    = mj2 + ncoefj;
        indexj = mj1:mj2;
        betaestfdj     = putcoef(betafdj, betacoef(indexj));
        betaestfdPar   = putfd(betafdParj, betaestfdj);
        betaestcell{j} = betaestfdPar;
    end
    
    %  set up fd object for predicted values
    
    yhatmat = zeros(nfine,N);
    for j=1:p
        xmat    = eval_fd(tfine, xfdcell{j});
        betafdj = getfd(betaestcell{j});
        betavec = eval_fd(tfine, betafdj);
        yhatmat = yhatmat + xmat.*(betavec*ones(1,N));
    end
    yhatfdobj = data2fd(yhatmat, tfine, ybasisobj);
    
    %  In order save time, the following computation only happens 
    %    if at least argument Y2CMAP is supplied.
    
    if nargin >= 4
        
        %  compute linear mapping c2bmap takinging coefficients for 
        %  response into coefficients for regression functions
        
        betabasisobj = getbasis(getfd(betacell{1}));
        Jmatphipsi = inprod_basis(ybasisobj,betabasisobj);
        
        ybasismat    = eval_basis(tfine, ybasisobj);
        basisprodmat = zeros(ncoef,ynbasis*N);
        mj2 = 0;
        for j=1:p
            betafdParj = betacell{j};
            betafdj    = getfd(betafdParj);
            ncoefj     = size(getcoef(betafdj),1);
            betabasisj = getbasis(betafdj);
            bnbasis    = getnbasis(betabasisj);
            bbasismat  = eval_basis(tfine, betabasisj);
            xfdj       = xfdcell{j};
            xcoefmat   = getcoef(xfdj);
            xbasisj    = getbasis(xfdj);
            xbasismat  = eval_basis(tfine, xbasisj);
            xnbasis    = getnbasis(xbasisj);
            mj1 = mj2 + 1;
            mj2 = mj2 + ncoefj;
            indexj = mj1:mj2;
            %  inner products of beta basis and response basis
            %    weighted by covariate basis functions
            temparray = zeros(bnbasis,ynbasis,xnbasis);
            for k=1:ynbasis
                wtmatk = ybasismat(:,k)*ones(1,bnbasis);
                temparray(:,k,:) = deltat.*(bbasismat.*wtmatk)'*xbasismat;
            end
            mk2 = 0;
            for k=1:ynbasis
                mk1 = mk2 + 1;
                mk2 = mk2 + N;
                basisprodmat(indexj,mk1:mk2) = ...
                    squeeze(temparray(:,k,:))*xcoefmat;
            end
        end    
        
        %  estimate SIGMAE if required
        
        if nargin == 4
            tfine    = linspace(rangeval(1), rangeval(2), size(y2cmap,2));
            ymat     = eval_fd(tfine, yfdobj);
            yhatmat  = eval_fd(tfine, yhatfdobj);
            residmat = ymat - yhatmat;
            sigmae   = cov(residmat');
        end
        
        %  check dimensions of Y2CMAP
        
        y2cdim = size(y2cmap);
        if y2cdim(1) ~= ynbasis | y2cdim(2) ~= size(sigmae,1)
            error('Dimensions of Y2CMAP not correct.');
        end
        
        c2bmap    = Cmatinv*basisprodmat;
        VarCoef   = y2cmap*sigmae*y2cmap';
        CVariance = kron(VarCoef,eye(N));
        bvar = c2bmap*CVariance*c2bmap';
        mj2 = 0;
        for j=1:p
            betafdParj = betacell{j};
            betafdj    = getfd(betafdParj);
            ncoefj     = size(getcoef(betafdj),1);
            mj1 = mj2 + 1;
            mj2 = mj2 + ncoefj;
            indexj = mj1:mj2;
            betabasisj = getbasis(betafdj);
            bnbasis    = getnbasis(betabasisj);
            bbasismat  = eval_basis(tfine, betabasisj);
            bvarj      = bvar(indexj,indexj);
            bstderrj   = sqrt(diag(bbasismat*bvarj*bbasismat'));
            bstderrfdj = data2fd(bstderrj, tfine, betabasisj);
            betastderrcell{j} = bstderrfdj;
        end
    else
        c2bmap         = [];
        bvar           = [];
        betastderrcell = {};
    end
    
elseif strcmp(class(yfdPar),'double')
    
    %  ----------------------------------------------------------------
    %                   YFDPAR is scalar or multivariate
    %  ----------------------------------------------------------------
    
    ymat = yfdPar;
    N = size(ymat,1);
    
    %  get number of independent variables 
    
    p = length(xfdcell);
    
    if length(betacell) ~= p
        error(['Number of regression coefficients does not match', ...
                ' number of independent variables.']);
    end
    
    %  check each cell.  If the object is a functional data object,
    %  it is converted to a multivariate object
    
    Zmat  = [];
    Rmat  = [];
    pjsum = 0;
    pjvec = zeros(p,1);
    for j=1:p
        xfdj = xfdcell{j};
        if strcmp(class(xfdj), 'fd')
            if ~(strcmp(class(betacell{j}),'fdPar') | ...
                 strcmp(class(betacell{j}),'fdpar'))
                error(['BETACELL{',num2str(j),'} does not contain' ...
                       ' an fdPar object when XFDCELL{',num2str(j), ...
                       '} is functional.']);
            end
            xcoef  = getcoef(xfdj);
            Nj = size(xcoef,2);
            if Nj ~= N
                error(['Coefficient matrix ',num2str(j), ...
                       ' has the wrong number of columns.']);
            end
            xbasis     = getbasis(xfdj);
            betafdParj = betacell{j};
            bbasis     = getbasis(getfd(betafdParj));
            bnbasis    = getnbasis(bbasis);
            pjvec(j)   = bnbasis;
            Jpsithetaj = inprod_basis(xbasis,bbasis);
            Zmat       = [Zmat,xcoef'*Jpsithetaj];
            lambdaj    = getlambda(betafdParj);
            Lfdj       = getLfd(betafdParj);
            if getestimate(betafdParj) & lambdaj > 0
                Rmatj = lambdaj.*eval_penalty(betafdParj);
            else
                Rmatj = zeros(pjvec(j));
            end
            Rmat  = [ [Rmat,            zeros(pjsum,bnbasis)]; ...
                      [zeros(bnbasis,pjsum), Rmatj] ];
            pjsum = pjsum + bnbasis;
        elseif strcmp(class(xfdj), 'double')
            Zmatj    = xfdj;
            [Nj,pj]  = size(Zmatj);
            pjvec(j) = pj;
            if Nj ~= N
                error(['Covariate matrix ',num2str(j), ...
                       ' has the wrong number of rows.']);
            end
            Zmat  = [Zmat,Zmatj];
            Rmatj = zeros(pj);
            Rmat  = [ [Rmat,  zeros(pjsum,pj)]; ...
                      [zeros(pj,pjsum), Rmatj] ];
            pjsum = pjsum + pj;
        else
            error(['Covariate ',num2str(j), ...
                   ' is neither a functional nor a multivariate object.']);
        end
    end
    
    %  -----------------------------------------------------------
    %          set up the linear equations for the solution
    %  -----------------------------------------------------------
    
    %  solve for coefficients defining BETA
    
    Cmat     = Zmat'*Zmat + Rmat;
    Cmatinv  = inv(Cmat);
    Dmat     = Zmat'*ymat;
    
    %  compute and print degrees of freedom measure
    
    df = trace(Zmat*Cmatinv*Zmat');
%     fprintf('  Degrees of freedom = %6.1f\n', df);
    
    betacoef = Cmatinv*Dmat;
    
    %  set up fdPar object for BETAFDPAR
    
    betaestcell = betacell;
    onebasis    = create_constant_basis([0,1]);
    mj2 = 0;
    for j=1:p
        mj1 = mj2 + 1;
        mj2 = mj2 + pjvec(j);
        indexj = mj1:mj2;
        betacoefj  = betacoef(indexj);
        betafdParj = betacell{j};
        if strcmp(class(xfdj), 'fd')        
            betafdj        = getfd(betafdParj);
            betaestfdj     = putcoef(betafdj, betacoefj);
            betaestfdParj  = putfd(betafdParj, betaestfdj);
            betaestcell{j} = betaestfdParj;
        else
            betaestfdj     = fd(betacoefj',onebasis);
            betaestfdParj  = putfd(betafdParj, betaestfdj);
            betaestcell{j} = betaestfdParj;
        end
    end
    
    %  set up fd object for predicted values
    
    yhatmat = zeros(N,1);
    for j=1:p
        xfdj = xfdcell{j};
        if strcmp(class(xfdj), 'fd')    
            xbasis  = getbasis(xfdj);
            xnbasis = getnbasis(xbasis);
            xrng    = getbasisrange(xbasis);
            nfine   = max(501,10*xnbasis+1);
            tfine   = linspace(xrng(1), xrng(2), nfine)';
            deltat  = tfine(2)-tfine(1);
            xmat    = eval_fd(tfine, xfdj);
            betafdj = getfd(betaestcell{j});
            betamat = eval_fd(tfine, betafdj);
            yhatmat = yhatmat + deltat.*(xmat'*betamat    - ...
                      0.5.*(xmat(1,    :)'*betamat(1)     + ...
                            xmat(nfine,:)'*betamat(nfine)));
        else
            betaj   = getcoef(getfd(betaestcell{j}));
            yhatmat = yhatmat + xfdj*betaj';
        end
    end
    yhatfdobj = yhatmat;
    
    %  In order save time, the following computation only happens 
    %    if the complete set of arguments is supplied.
    
    if nargin >= 4
        
        %  compute sigmae if not input
        
        if nargin == 4
            resmat = ymat - yhatmat;
            sigmae = cov(yhatmat);
        end
        
        %  compute linear mapping c2bmap takinging coefficients for 
        %  response into coefficients for regression functions
        
        c2bmap = Cmatinv*Zmat';
        y2bmap = c2bmap;
        bvar   = y2bmap*sigmae*y2bmap';
        mj2 = 0;
        for j=1:p
            mj1 = mj2 + 1;
            mj2 = mj2 + pjvec(j);
            indexj = mj1:mj2;
            bvarj      = bvar(indexj,indexj);
            if strcmp(class(xfdj), 'fd')        
                betafdParj = betacell{j};
                betafdj    = getfd(betafdParj);
                ncoefj     = size(getcoef(betafdj),1);
                betabasisj = getbasis(betafdj);
                bnbasis    = getnbasis(betabasisj);
                betarng    = getbasisrange(betabasisj);
                nfine      = max(501,10*bnbasis+1);
                tfine      = linspace(betarng(1), betarng(2), nfine)';
                bbasismat  = eval_basis(tfine, betabasisj);
                bstderrj   = sqrt(diag(bbasismat*bvarj*bbasismat'));
                bstderrfdj = data2fd(bstderrj, tfine, betabasisj);
                betastderrcell{j} = bstderrfdj;
            else
                bstderrj   = sqrt(diag(bvarj));
                bstderrfdj = fd(bstderrj', onebasis);
                betastderrcell{j} = bstderrfdj;
            end
        end
    else
        c2bmap         = [];
        bvar           = [];
        betastderrcell = {};
    end
    
else
    %  YFDOBJ is neither functional nor multivariate
    error('YFDOBJ is neither functional nor multivariate.');
end
