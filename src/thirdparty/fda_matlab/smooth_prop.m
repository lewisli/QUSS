function [fdobj, df, gcv, coef, SSE, penmat, y2cMap] = ...
    smooth_prop(argvals, freq, N, basisobj, wtvec, Lfdobj, lambda, fdnames)
%SMOOTH_PROP  Smooths binomial observations over ARGVALS using penalized 
%  basis expansion of the log-odds of expected proportions.
%  Like SMOOTH_MONOTONE, this function only smooths one relationship
%     at a time. 
%  Arguments for this function:
%
%  FREQ     ... an array containing frequencies.
%  N        ... an array of sample sizes for the binomial observations.
%  ARGVALS  ... A set of argument values, set by default to equally spaced on
%               the unit interval (0,1).
%  BASISOBJ ... A basis.fd object created by function create_basis.fd.
%  WTVEC    ... A vector of N weights, set to one by default, that can
%               be used to differentially weight observations.
%  LFDOBJ   ... The order of derivative or a linear differential
%               operator to be penalized.
%  LAMBDA   ... The smoothing parameter determining the weight to be
%               placed on the roughness penalty.  
%               This is 0 by default.
%  FDNAMES  ... A cell of length 3 with names for
%               1. argument domain, such as 'Time'
%               2. replications or cases
%               3. the function.
%  Returns:
%    FDOBJ  ...  an object of class fd containing coefficients
%    DF     ...  a degrees of freedom measure
%    GCV    ...  a measure of lack of fit discounted for df.
%    COEF   ...  the coefficient matrix for the basis function
%                  expansion of the smoothing function
%    SSE    ...  the error sums of squares
%    PENMAT ...  the penalty matrix
%    y2cMat ...  the matrix mapping the data to the coefficients

%  last modified 25 November 2003

if nargin < 4
    error('There is not at least four arguments.');
end

n = length(argvals);

%  set default argument values

if nargin < 8
    fdnames{1} = 'time';
    fdnames{2} = 'reps';
    fdnames{3} = 'values';
end

if nargin < 7, lambda = 0;           end
if nargin < 6, Lfdobj = int2Lfd(2);  end;
if nargin < 5, wtvec  = ones(n,1);   end

%  check LFD

Lfdobj = int2Lfd(Lfdobj);
nderiv = getnderiv(Lfdobj);

%  check BASIS

if ~isa_basis(basisobj)
    error('BASIS is not a basis object.');
end

nbasis   = getnbasis(basisobj);
onebasis = ones(1,nbasis);

%  check WTVEC

sizew = size(wtvec);
if (length(sizew) > 1 & sizew(1) > 1 & sizew(2) > 1) | ...
      length(sizew) > 2
    error ('WTVEC must be a vector.');
end
if length(sizew) == 2 & sizew(1) == 1
    wtvec = wtvec';
end
if length(wtvec) ~= n
    error('WTVEC of wrong length');
end
if min(wtvec) <= 0
    error('All values of WTVEC must be positive.');
end

%  check LAMBDA

if lambda < 0
    warning ('Value of LAMBDA was negative, and 0 used instead.');
    lambda = 0;
end

%  check N 

sizeN = size(N);
if sizeN(1) ~= n
    error(['Number of arguments differs from ', ...
           'number of sample sizes.']);
end
Nfloor = floor(N);
if any(Nfloor ~= N)
    error('Samples sizes are not all integers.');
end
if any(Nfloor <= 0)
    error('Sample sizes are not all positive.');
end

%  check FREQ

sizefreq = size(freq);
if sizefreq(1) ~= n
    error(['Number of arguments differs from ', ...
           'number of frequencies.']);
end
freqfloor = floor(freq);
if any(freqfloor ~= freq)
    error('Frequencies are not all integers.');
end
if any(freq > N)
    error('Frequencies are not less than or equal to N.');
end
if any(freq < 0)
    error('Negative frequencies are found.');
end

basismat  = getbasismatrix(argvals, basisobj);

if n >= nbasis | lambda > 0
    
    %  The following code is for the coefficients completely determined
    
    basisw = basismat .* (wtvec * ones(1,nbasis));
    Bmat   = basisw' * basismat;
    Bmat0  = Bmat;
    
    Dmat = basisw' * y;
    
    if lambda > 0
        %  smoothing required, set up coefficient matrix for normal equations
        afdcell = getafdcell(Lfdobj);  % multiplier(s) of forcing function(s)
        ufdcell = getufdcell(Lfdobj);  % forcing function(s)
        if ~isempty(afdcell) & ~isempty(ufdcell)
            %  the linear differential operator is not homogeneous.  
            %  first set up the homogeneous counterpart
            wtfdcell = getwfdcell(Lfdobj);
            Lfdhom = Lfd(nderiv, wtfdcell);
            %  evaluate the penalty matrix for the homogeneous operator
            penmat = eval_penalty(basisobj, Lfdhom);
            %  set up the part of the roughness penalty affected by the
            %  presence of forcing function(s)
            penvec = zeros(nbasis,1);
            nforce = length(ufdcell);
            for k=1:nforce
                afdk = getfd(afdcell{k});
                ufdk = ufdcell{k};
                ffdk = times(afdk,ufdk);
                fvec = inprod(basisobj, ffdk, Lfdhom, int2Lfd(0));
                penvec = penvec - mean(fvec,2);
            end
        else
            %  here the linear differential operator is homogeneous
            %  only the penalty matrix is needed.
            penmat = eval_penalty(basisobj, Lfdobj);
            penvec = zeros(nbasis,1);
        end
        Bnorm   = sqrt(sum(sum(Bmat.^2)));
        pennorm = sqrt(sum(sum(penmat.^2)));
        condno  = pennorm/Bnorm;
        if lambda*condno > 1e12
            lambda = 1e12/condno;
            warning(['lambda reduced to ',num2str(lambda), ...
                    ' to prevent overflow']);
        end
        Bmat = Bmat + lambda .* penmat;
        if ~all(penvec == 0)
            %  if the linear differential operator is nonhomogeneous
            %  use PENVEC to alter the right side of the equation.
            Dmat = Dmat + lambda.*(penvec*ones(1,ncurves));
        end
        
    else
        
        penmat = zeros(nbasis);
        
    end
    
    %  compute inverse of Bmat
    
    if is_diag(Bmat)
        Bmatinv = diag(1./diag(Bmat));
    else
        Bmatinv = inv(Bmat);
    end
    
    %  compute map from y to c
    
    y2cMap = Bmatinv * basisw';
    
    %  compute degrees of freedom of smooth
    
    df = sum(diag(Bmatinv * Bmat0));
    
    %  solve normal equations for each observation
    
    if ndim < 3
        coef = Bmatinv * Dmat;
    else
        for ivar = 1:nvar
            coef(:,:,ivar) = Bmatinv * Dmat(:,ivar);
        end
    end
    
else
    
    %  The following code is for the underdetermined coefficients:
    %     the number of basis functions exceeds the number of argument values.
    %  No smoothing is used.  
    
    [Qmat,Rmat] = qr(basismat');
    Q1mat  = Qmat(:,1:n);
    Q2mat  = Qmat(:,((n+1):nbasis));
    Hmat   = eval_penalty(basisobj);
    Q2tHmat   = Q2mat' * Hmat;
    Q2tHQ2mat = Q2tHmat * Q2mat;
    Q2tHQ1mat = Q2tHmat * Q1mat;
    if ndim < 3
        z1mat = symsolve(Rmat, y);
        z2mat = symsolve(Q2tHQ2mat, Q2tHQ1matz1mat);
        coef = Q1mat * z1mat + Q2mat * z2mat;
    else
        for ivar = 1:nvar
            z1mat = symsolve(Rmat, y(:,:,ivar));
            z2mat = symsolve(Q2tHQ2mat, Q2tHQ1mat*z1mat);
            coef(:,:,ivar) = Q1mat * z1mat + Q2mat * z2mat;
        end
    end
    y2cMap = eye(n);
    df = n;
end

%  compute error sum of squares

if ndim < 3
    yhat = basismat * coef;
    SSE = sum(sum((y - yhat).^2));
else
    SSE = 0;
    for ivar = 1:nvar
        yhat = basismat * coef(:,:,ivar);
        SSE = SSE + sum(sum((y(:,:,ivar) - yhat).^2));
    end
end

%  compute  GCV index

if df < n
    gcv = (SSE/n)/(nvar*(n - df)/n)^2;
else
    gcv = NaN;
end

fdobj = fd(coef, basisobj, fdnames);


