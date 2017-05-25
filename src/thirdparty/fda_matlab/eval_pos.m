function evalarray = eval_pos(evalarg, fdobj, Lfdobj)
%  Evaluates a derivative of a single positive functional data 
%  observation. 
%  A positive functional data object h  is = the form
%           h(x) = (exp fdobj)(x)
%  Note that the linear differential operator object LFDOBJ
%  MUST be an integer in the range 0 to 1.
%  Note that the first two arguments may be interchanged.
%
%  Arguments:
%  EVALARG ... A vector of values at which all functions are to 
%              evaluated.
%  FDOBJ   ... Functional data object.  It must define a single
%              functional data observation.
%  LFDOBJ  ... A linear differential operator object
%              applied to the functions that are evaluated.
%
%  Returns:  An array of function values corresponding to the evaluation
%              arguments in EVALARG

%  Last modified 30 January 2003

if nargin < 2
    error('Number of arguments is less than 2.');
end

%  check LFDOBJ and convert an integer to Lfd if needed.

if nargin < 3 
    %  set default LFDOBJ to 0
    Lfdobj = int2Lfd(0); 
else
    %  check LFDOBJ
    if isnumeric(Lfdobj)
        % if integer, check and convert to Lfd
        nderiv = Lfdobj;
        if nderiv ~= round(Lfdobj)
            error('LFDOBJ numeric but not an integer.');
        end
        if nderiv < 0
            error('LFDOBJ an integer but negative.');
        end
        Lfdobj = int2Lfd(nderiv);
    else
        if ~isa_Lfd(Lfdobj)
            error (['Argument LFDOBJ is neither a functional data object', ...
                    ' nor an integer.']);
        end
        nderiv = getnderiv(Lfdobj);
        if ~isinteger(Lfdobj)
            error('LFD is not a D^m operator.');
        end        
    end
end

%  Exchange the first two arguments if the first is an FD object
%    and the second numeric

if isnumeric(fdobj) & isa_fd(evalarg)
    temp    = fdobj;
    fdobj   = evalarg;
    evalarg = temp;
end

%  Check the arguments

if ~(isnumeric(evalarg))
    error('Argument EVALARG is not numeric.');
end

%  transpose EVALARG if necessary to make it a column vector

evaldim = size(evalarg);
if evaldim(1) == 1 & evaldim(2) > 1  
    evalarg = evalarg';  
    evaldim = size(evalarg);
end

%  check EVALARG

sizeevalarg = size(evalarg);
if sizeevalarg(1) > 1 & sizeevalarg(2) > 1
    error('Argument EVALARG is not a vector.');
end
evalarg = evalarg(:);

%  check FDOBJ

if ~isa_fd(fdobj)
    error('Argument FD is not a functional data object.');
end

%  check LFDOBJ

if ~(isa_Lfd(Lfdobj))
    error('LFDOBJ is not linear differential operator object.');
end

%  Extract information about the basis

basisfd  = getbasis(fdobj);
nbasis   = getnbasis(basisfd);
rangeval = getbasisrange(basisfd);
onerow   = ones(1,nbasis);

%  determine the highest order of derivative NDERIV required

nderiv = getnderiv(Lfdobj);

if nderiv > 1
    error('EVAL_POS not programmed for NDERIV > 1.');
end

%  Set up coefficient array for FD

coef  = getcoef(fdobj);

%  Case where EVALARG is a vector of values to be used for all curves

evalarg(evalarg < rangeval(1)-1e-10) = NaN;
evalarg(evalarg > rangeval(2)+1e-10) = NaN;
basismat = getbasismatrix(evalarg, basisfd);
fdvec    = exp(basismat*coef);

%  If a differential operator has been defined in LFDOBJ, compute
%  the weighted combination of derivatives

if nderiv > 0
    basismat = eval_basis(evalarg, basisfd, Lfdobj);
    afd  = getafd(Lfdobj);
    ufd  = getufd(Lfdobj);
    avec = eval_fd(evalarg, afd);
    uvec = eval_fd(evalarg, ufd);
    fvec = avec.*uvec;
    evalarray = fdvec.*(basismat*coef + fvec);
else
    evalarray = fdvec;
end

