%  test of smoothing by estimating a linear differential operator

addpath('c:\Matlab\smooth_pda')

%  ------------------------------------------------------------
%    A forced second order constant coefficient DIFE
%        (r^2 + w^2)x  
%    Annihilating operator is L = (w^2 + r^2)x + 2rDx + D^2x
%    Range is [0,4 pi], forcing function is unit step function
%    with step at 2 pi.  
%  ------------------------------------------------------------

%  set order of the differential operator

DIFEorder = 2;

%  generate data

N =  1;  %  number of curves

n = 101;  %  number of sampling values per curve

Tlim  = 4*pi;
range = [0,Tlim];
tval  = linspace(range(1),range(2),n)';  %  sampling points
delta = Tlim/(n-1);

%  set spring constant 

w = 2;

%  set decay constant

r = 0.2;

%  coefficients for sine and cosine

c(1) = 1;
c(2) = 1;

%  coefficient for forcing function

alpha = 4;

%  step time

steptime = 2*pi;

%  compute errorless curves Y0 and their NORDER derivatives

[y0, xhom, uvec, fvec] = harmfn(tval, c, w, r, steptime, alpha);

%  plot errorless functions and their derivatives and 
%  the result of applying the linear differential operator to them

subplot(1,1,1)
plot(tval, y0,'b-', tval, xhom, 'b--', ...
     [0,        steptime], [0,0], 'b:', ...
     [steptime, steptime], [0,1], 'b:', ...
     [steptime, 4*pi    ], [1,1], 'b:')
xlabel('\fontsize{19} t')
ylabel('\fontsize{19} x(t)')
axis([0,4*pi,-1,1.5])

%  add some error

% sigma = 0.0;
sigma = 0.002;
y = y0*ones(1,N) + randn(n,N)*sigma;

%  plot some of the curves and data

plot(tval, y, '.', tval, y0, '-', ...
     tval, uvec, 'r--')
axis([0,4*pi,-1.5,2])
legend('\fontsize{16} data', 'x(t)', 'u(t)')

%  basis for the analysis of the second order forced problem
%  put three coincident knots at 2*pi to allow for derivative
%  discontinuity

norder =  4;
breaks = [linspace(0, 2*pi, 25), 2*pi, linspace(2*pi, 4*pi, 25)];
nbasis = length(breaks) + norder - 2;
basisobj = create_bspline_basis(range, nbasis, norder);

%  generate quadrature values in BASISOBJ
%  set up Simpson's rule over [0,2*pi] and [2*pi,4*pi]

quadpts = [linspace(0,2*pi,481)'; linspace(2*pi,4*pi,481)'];
quadwts = ones(481,1);
quadwts(2:2:480) = 4;
quadwts(3:2:479) = 2;
quadwts = (2*pi/480).*quadwts/3;
quadwts = [quadwts; quadwts];
quadvals = [quadpts, quadwts];

basisobj = putquadvals(basisobj, quadvals);

%  compute basis fn & deriv. values and load into BASISOBJ

for ivalue=1:(DIFEorder+1)
    basismat       = eval_basis(quadpts, basisobj, ivalue-1);
    values{ivalue} = basismat.*(sqrt(quadwts)*ones(1,nbasis));
end

basisobj = putvalues(basisobj, values);

%  check penalty matrix for accuracy

penmat1 = values{3}'*values{3};
penmat2 = eval_penalty(basisobj, 2);

max(max(abs(full(penmat1-penmat2))))/max(max(abs(full(penmat2))))

%  constant basis for second order forced operator

%  constant basis

nbasisL = 1;
basisL  = create_constant_basis(range);  

%  load quadrature values into BASISL

basisL = putquadvals(basisL,quadvals);

%  compute basis fn & deriv. values and load into BASISOBJ

basisLmat  = eval_basis(quadvals(:,1), basisL);
valuesL{1} = basisLmat;

basisL = putvalues(basisL, valuesL);

%  set up BWTCELL to specify weight fn b(t)
%  first weight function

clear bwtcell

wfd0       = fd(0, basisL);
bwtcell{1} = fdPar(wfd0);
bwtcell{2} = fdPar(wfd0);

%  set up AWTCELL to specify weight fn. a(t)

nbasisA = 1;
basisA  = create_constant_basis(range);

%  load quadrature values into BASISA

basisA = putquadvals(basisA,quadvals);

%  compute basis fn & deriv. values and load into BASISOBJ

basisAmat  = eval_basis(quadvals(:,1), basisA);
basisAmat  = basisAmat.*sqrt(quadwts);
valuesA{1} = basisAmat;

basisA = putvalues(basisA, valuesA);

afd0    = fd(-alpha, basisA);
awtcell{1} = fdPar(afd0);

%  set up UFDCELL to specify a single forcing function.  

nbasisU = 2;
norderU = 1;
basisU  = create_bspline_basis(range,nbasisU,norderU);
ufd     = data2fd(uvec,tval,basisU);
ufdcell{1} = ufd;     

%  -----------------------------------------------------------
%                  do PDA using PROFPDA_SISO
%  -----------------------------------------------------------

%  set up the analysis using profPDA_SISO
%  these set up a weighted cross-product matrix Bmat
%  corresponding to the basis functions
%  and a weighted cross product Dmat of the basis functions
%  with the observed function.  Defining these outside 
%  of function profPDA saves considerable computation since
%  profPDA is designed to be called many times by an 
%  optimization function

basismat = getbasismatrix(tval, basisobj);
Bmat     = basismat' * basismat;

%  set up options for FMINUNC

options = optimset('LargeScale', 'off', 'Display', ...
                   'iter', 'MaxIter', 20, 'GradObj', 'on', ...
                   'TolFun', 1e-5, 'TolCon', 1e-5, ...
                   'TolX',   1e-5, 'TolPCG', 1e-5);
               
%  define initial values

bvec0 = zeros(3,1);
bvec0(1) = w^2 + r^2;
bvec0(2) = 2.0*r;
bvec0(3) = -alpha;

%  set up range of log10 lambda's

loglam = (-2:2)';
nlam   = length(loglam);
Nsam   = 1;
df     = zeros(nlam,Nsam);
gcv    = zeros(nlam,Nsam);
SSE    = zeros(nlam,Nsam);
PENSSE = zeros(nlam,Nsam);
bmat   = zeros(nlam,length(bvec0),Nsam);
% hessdg = zeros(nlam,length(bvec0),N);

% for i=1:10
    Dmat     = basismat' * y;
    for ilam=1:nlam
        lambda = 10^loglam(ilam);
        [bvec, fval, exitflag, output, grad] = ...
            fminunc(@profPDA, bvec0, options, ...
                    y, basisobj, basismat, Bmat, Dmat, ...
                    bwtcell, awtcell, ufdcell, lambda, 1);
        bmat(ilam,:,1) = bvec';
        [SSE(ilam,1), DSSE, PENSSE(ilam,1), fdobj, ...
                df(ilam,1), gcv(ilam,1)] = ...
            profPDA(bvec, y, basisobj, basismat, Bmat, Dmat, ...
                    bwtcell, awtcell, ufdcell, lambda);
%         hessdg(ilam,:) = diag(inv(hessian))';
        [loglam(ilam), df(ilam,1), gcv(ilam,1)]
    end
    disp(['Curve ',num2str(1)])
    [loglam, df(:,1), gcv(:,1)]
%     pause
% end

%  ----------------------------------------------------------
%                 set up the analysis using PDACELL
%  ----------------------------------------------------------

%  use profPDA

load secondorder

bvec0 = [(r^2+w^2)*ones(nbasisL,1); ...
         2*r*ones(nbasisL,1);       ...
         -alpha];

bvec0 = [zeros(nbasisL,1); zeros(nbasisL,1); 0];

%  set up the analysis using profPDA

basismat = getbasismatrix(tval, basisobj);
basisw   = basismat .* (wtvec * ones(1,nbasis));
Bmat     = basisw' * basismat;
Dmat     = basisw' * y;

%  A trial run with profPDA

lambda = 2;

[SSE, DSSE, PENSSE, fdobj, df, gcv] = ...
    profPDA(bvec0, y, basisobj, basismat, Bmat, Dmat, ...
            bwtcell, awtcell, ufdcell, lambda, 1);

[SSE, PENSSE, DSSE] 
[df, gcv]

% bvec =
%      grad    ~grad
%      
%     4.1316  4.1308
%     0.4312  0.4330
%    -2.0071 -2.0070

%  A run with a grid of parameter values

lambda = 2;

bpar0 =  3.5: 1.0:  4.5;
bpar1 =  0.3: 0.1:  0.5;
avec  = -3.0: 1.0: -1.0;
SSEstore = zeros(length(bpar0),length(bpar1),length(avec));

i=0;
for bval0=bpar0
    i = i + 1;
    j = 0;
    for bval1=bpar1
        j = j + 1;
        k = 0;
        for aval=avec;
            k = k + 1;
            bvec = [bval0, bval1, aval];
            [SSE, DSSE, PENSSE, fdobj, df, gcv] = ...
                profPDA(bvec, y, basisobj, basismat, Bmat, Dmat, ...
                bwtcell, awtcell, ufdcell, lambda);
            full([bval0, bval1, aval, SSE, PENSSE, df, gcv])
            SSEstore(i,j,k) = SSE;
        end
    end
end

SSEstore
 
%  set options for fminunc

options = optimset('LargeScale', 'off', 'Display', ...
                   'iter', 'MaxIter', 5, 'GradObj', 'on');

%  set up range of log10 lambda's

loglam = 2;
loglam = (2:0.5:4)';
nlam   = length(loglam);

df     = zeros(nlam,1);
gcv    = zeros(nlam,1);
SSE    = zeros(nlam,1);
PENSSE = zeros(nlam,1);
bmat   = zeros(nlam,length(bvec0));
hessdg = zeros(nlam,length(bvec0));

for ilam=1:nlam
    lambda = 10^loglam(ilam);
    [bvec, fval, exitflag, output, grad] = ...
        fminunc(@profPDA, bvec0, options, ...
        y, basisobj, basismat, Bmat, Dmat, ...
        bwtcell, awtcell, ufdcell, lambda, 1);
    bmat(ilam,:) = bvec';
    [SSE(ilam), DSSE, PENSSE(ilam), fdobj, ...
            df(ilam), gcv(ilam)] = ...
        profPDA(bvec, y, basisobj, basismat, Bmat, Dmat, ...
        bwtcell, awtcell, ufdcell, lambda);
    hessdg(ilam,:) = diag(inv(hessian))';
    [loglam(ilam), df(ilam), gcv(ilam)]
end

%  display RMS and penalized RMS

sqrt(SSE./((N*(n-df))))

%  display DF and GCV values

[df, gcv]

%  display coefficient vectors

[bmat]

bindex = 1:nlam;

bmat(bindex,:,:)

% output some results

for ilam=bindex    
    bvec = bmat(ilam,:)';
    
    %  set up linear differential operator
    
    [Lfdobj, bfdcell, afdcell] = ...
        bvec2Lfd(bvec, bwtcell, awtcell, ufdcell);
    
    lambda   = 10^loglam(ilam);
    fdParobj = fdPar(basisobj, Lfdobj, lambda);

     %  plot the operator
    
    plot(Lfdobj)
    title(['Lambda = ',num2str(lambda)]);
    pause;
    
    %  smooth the data with the operator
    
    [fdobj, dfi, gcvi, coef, SSEi, penmati] = ...
                                smooth_basis(tval, y, fdParobj);
    
    %  plot the data and fit
    
    subplot(1,1,1)
    yhat = eval_fd(tval, fdobj);
    yres = y - yhat;
    plot(tval, y, '.', tval, yhat, 'b-')
    xlabel('\fontsize{16} t')
    ylabel('\fontsize{16} x(t)')
    title(['Lambda = ',num2str(lambda)]);
    legend('\fontsize{16} Data', 'Estimate')
    pause;
    
    plot(tval, eval_fd(tval, fdobj, 1), 'b-')
    xlabel('\fontsize{16} t')
    ylabel('\fontsize{16} Dx(t)')
    title(['Lambda = ',num2str(lambda)]);
    pause;
end

%  ------------------------------------------------------------
%    A monotone function  
%  ------------------------------------------------------------

%  set order of the differential operator

DIFEorder = 2;

%  generate data

N =  1;  %  number of curves

n = 101;  %  number of sampling values per curve

range = [0,1];
tval  = linspace(range(1),range(2),n)';  %  sampling points
delta = 1/(n-1);
wtvec = ones(n,1);  %  vector of weights

%  set up function W(t)

Wval = 2.*sin(4*pi*tval);

%  linear coefficients 

c(1) = 0
c(2) = 1.5;

%  compute errorless curves Y0 and their NORDER derivatives

y0 = c(1) + c(2).*delta.*cumtrapz(exp(Wval));

%  plot errorless functions and their derivatives and 
%  the result of applying the linear differential operator to them

subplot(1,1,1)
plot(tval, y0,'-')
ylabel('x')

%  add some error

sigma = 0.2;
y = y0 + randn(n,N)*sigma;

%  plot some of the curves and data

subplot(1,1,1)
for i=1:N
    plot(tval, y(:,i), '.', tval, y0(:,i), '-')
    %axis([0,2,-.5,1.5])
    title(['record ',num2str(i)])
    pause;
end

%  basis for the analysis of the problem

nbasis = 24;
basisobj  = create_bspline_basis(range,nbasis,5);

%  basis for operator

%  B-spline basis

nbasisL = 13;
norderL = 4
basisL  = create_bspline_basis(range, nbasisL, norderL);  
wfd0    = fd(zeros(nbasisL,1), basisL);

%  set up the analysis using PDACELL

%  set up BWTCELL to specify weight fn b(t)
%  first weight function

fd       = wfd0;
estimate = 0;
lambda   = 0;
Lfdobj   = int2Lfd(1);
bwtcell{1} = fdPar(wfd0, Lfdobj, lambda, estimate);

%  second weight function

fd       = wfd0;
estimate = 1;
lambda   = 0;
Lfdobj   = int2Lfd(1);
bwtcell{2} = fdPar(wfd0, Lfdobj, lambda, estimate);

%  set up AWTCELL and UFDCELL

awtcell = {};
ufdcell = {};

%  set up the analysis using profPDA

basismat = getbasismatrix(tval, basisobj);
basisw   = basismat .* (wtvec * ones(1,nbasis));
Bmat     = basisw' * basismat;
Dmat     = basisw' * y;

%  set up initial values of coefficients

wfd0  = data2fd(c(2).*cos(2*pi*tval), tval, basisL);
bvec0 = -getcoef(wfd0);

options = optimset('LargeScale', 'off', 'Display', ...
                   'iter', 'MaxIter', 10, 'GradObj', 'on');

%  set up range of log10 lambda's

loglam = 2;
loglam = (0:1:3)';
nlam   = length(loglam);

df     = zeros(nlam,1);
gcv    = zeros(nlam,1);
SSE    = zeros(nlam,1);
PENSSE = zeros(nlam,1);
bmat   = zeros(nlam,length(bvec0));
% hessdg = zeros(nlam,length(bvec0));

for ilam=1:nlam
    lambda = 10^loglam(ilam);
    [bvec, fval, exitflag, output, grad] = ...
        fminunc(@profPDA, bvec0, options, ...
        y, basisobj, basismat, Bmat, Dmat, ...
        bwtcell, awtcell, ufdcell, lambda, 1);
    bmat(ilam,:) = bvec';
    [SSE(ilam), DSSE, PENSSE(ilam), fdobj, ...
            df(ilam), gcv(ilam)] = ...
        profPDA(bvec, y, basisobj, basismat, Bmat, Dmat, ...
        bwtcell, awtcell, ufdcell, lambda);
    hessdg(ilam,:) = diag(inv(hessian))';
    [loglam(ilam), df(ilam), gcv(ilam)]
end

%  test of smoothing by estimating a linear differential operator

%  ------------------------------------------------------------
%    A unforced second order constant coefficient DIFE
%        (r^2 + w^2)x  
%    Annihilating operator is L = (w^2 + r^2)x + 2rDx + D^2x
%    Range is [0,8 pi], forcing function is unit step function
%    with step at 2 pi.  
%  ------------------------------------------------------------

%  set order of the differential operator

DIFEorder = 2;

%  generate data

N =  1;  %  number of curves

n = 201;  %  number of sampling values per curve

Tlim  = 8*pi;
range = [0,Tlim];
tval  = linspace(range(1),range(2),n)';  %  sampling points
delta = Tlim/(n-1);

%  set decay constant

r = 0.05;

%  set spring constant 

w = sqrt(1 - r^2);

%  coefficients for sine and cosine

c(1) = 1;
c(2) = 0;

%  compute errorless curves Y0 and their NORDER derivatives

[y0, xhom, uvec, fvec] = harmfn(tval, c, w, r);

%  plot errorless functions and their derivatives and 
%  the result of applying the linear differential operator to them

subplot(1,1,1)
plot(tval, y0,'b-', tval, xhom, 'b--', ...
     [0,        steptime], [0,0], 'b:', ...
     [steptime, steptime], [0,1], 'b:', ...
     [steptime, 4*pi    ], [1,1], 'b:')
xlabel('\fontsize{19} t')
ylabel('\fontsize{19} x(t)')
axis([0,Tlim,-1,1.5])

%  add some error

% sigma = 0.0;
sigma = 0.1;
y = y0*ones(1,N) + randn(n,N)*sigma;

%  plot some of the curves and data

plot(tval, y, '.', tval, y0, '-', ...
     tval, uvec, 'r--')
axis([0,Tlim,-1.5,2])
legend('\fontsize{16} data', 'x(t)', 'u(t)')

%  basis for the analysis of the second order forced problem
%  put three coincident knots at 2*pi to allow for derivative
%  discontinuity

norder = 6;
nbasis = norder + 14;
basisobj = create_bspline_basis(range, nbasis, norder);

%  generate quadrature values in BASISOBJ
%  set up Simpson's rule over [0,2*pi] and [2*pi,Tlim]

nquad = 2001;
quadpts = [linspace(0,T,nquad)'];
quadwts = ones(nquad,1);
quadwts(2:2:nquad-1) = 4;
quadwts(3:2:nquad-2) = 2;
quadwts = (T/(nquad-1)).*quadwts/3;
quadvals = [quadpts, quadwts];

basisobj = putquadvals(basisobj, quadvals);

%  compute basis fn & deriv. values and load into BASISOBJ

for ivalue=1:(DIFEorder+1)
    basismat       = eval_basis(quadpts, basisobj, ivalue-1);
    values{ivalue} = basismat.*(sqrt(quadwts)*ones(1,nbasis));
end

basisobj = putvalues(basisobj, values);

%  constant basis for second order forced operator

%  constant basis

nbasisL = 1;
basisL  = create_constant_basis(range);  

%  load quadrature values into BASISL

basisL = putquadvals(basisL,quadvals);

%  compute basis fn & deriv. values and load into BASISOBJ

basisLmat  = eval_basis(quadvals(:,1), basisL);
valuesL{1} = basisLmat;

basisL = putvalues(basisL, valuesL);

%  set up BWTCELL to specify weight fn b(t)
%  first weight function

clear bwtcell

wfd0       = fd(0, basisL);
bwtcell{1} = fdPar(wfd0);
bwtcell{2} = fdPar(wfd0);

%  set up AWTCELL

awtcell = {};

%  set up UFDCELL 

ufdcell = {};     

%  set up cross product matrix

basismat = getbasismatrix(tval, basisobj);
Bmat     = basismat' * basismat;
Dmat     = basismat' * y;

%  set up fitcell

fitstruct.y = y;
fitstruct.basisobj = basisobj;
fitstruct.basismat = basismat;
fitstruct.Bmat = Bmat;
fitstruct.Dmat = Dmat;
fitstruct.bwtcell = bwtcell;
fitstruct.awtcell = {};
fitstruct.ufdcell = {};
fitstruct.Dorder  = 2;
fitstruct.lambda  = 1e1;

fitcell{1} = fitstruct;

%  define initial values

bvec0 = zeros(2,1);
bvec0(1) = w^2 + r^2;
bvec0(2) = 2.0*r;

%  compute penalty matrix for these values

gradwrd = 0;
fitstruct = fitcell{1};
bwtcell = fitstruct.bwtcell;
penmat_SISO = full(eval_Rs(bwtcell, awtcell, ufdcell, basisobj, gradwrd));


%  compute solution for these values

gradwrd = 0;
[SSE, DSSE, PENSSE, fdobj, df, gcv] = ...
    profPDA_SISO(bvec0, fitcell, gradwrd);

%  set up options for FMINUNC

options = optimset('LargeScale', 'off', 'Display', ...
                   'iter', 'MaxIter', 20, 'GradObj', 'on', ...
                   'TolFun', 1e-5, 'TolCon', 1e-5, ...
                   'TolX',   1e-5, 'TolPCG', 1e-5);
               
%  define initial values

bvec0 = zeros(2,1);
bvec0(1) = w^2 + r^2;
bvec0(2) = 2.0*r;

%  set up range of log10 lambda's

loglam = 1;
nlam   = length(loglam);
Nsam   = 1;
df     = zeros(nlam,Nsam);
gcv    = zeros(nlam,Nsam);
SSE    = zeros(nlam,Nsam);
PENSSE = zeros(nlam,Nsam);
bmat   = zeros(nlam,length(bvec0),Nsam);

for ilam=1:nlam
    lambda = 10^loglam(ilam);
    fitstruct.lambda = lambda;
    fitcell{1} = fitstruct;
    [bvec, fval, exitflag, output, grad] = ...
        fminunc(@profPDA_SISO, bvec0, options, ...
        fitcell, 1);
    bmat(ilam,:,1) = bvec';
    [SSE(ilam,1), DSSE, PENSSE(ilam,1), fdobj, ...
            df(ilam,1), gcv(ilam,1)] = ...
        profPDA(bvec, y, basisobj, basismat, Bmat, Dmat, ...
        bwtcell, awtcell, ufdcell, lambda);
    [loglam(ilam), df(ilam,1), gcv(ilam,1)]
end
disp(['Curve ',num2str(1)])
[loglam, df(:,1), gcv(:,1)]


