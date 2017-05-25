%    test profPDAM function

%  -----------------------------------------------------------
%  set up a two-variable exponential problem
%  -----------------------------------------------------------

t = linspace(0,4,101)';

%  exponential functions, no forcing, two initial values

sigma = 0.2;
y10 =     exp(-2.*t);
y1  = y10 + sigma.*randn(101,1);
y20 =     exp(-t);
y2  = y20 + sigma.*randn(101,1);

%  exponential functions, unit forcing, two gain values

sigma = 0.2;
y10 = 1.*(1-exp(-2.*t));
y1  = y10 + sigma.*randn(101,1);
y20 = 1.*(1-exp(-t));
y2  = y20 + sigma.*randn(101,1);

subplot(2,1,1)
plot(t, y1, '.', t, y10, 'b-')
subplot(2,1,2)
plot(t, y2, '.', t, y20, 'b-')

%  set up basis for fitting first set of data

nbasis1   = 5;
norder1   = 4;
basisobj1 = create_bspline_basis([0,4],nbasis1,norder1);
% nbasis1   = 3;
% basisobj1 = create_monomial_basis([0,4],nbasis1);
basismat1 = eval_basis(t, basisobj1);

Bmat1     = basismat1'*basismat1;
Dmat1     = basismat1'*y1;

quadpts = [linspace(0,4,201)'];
quadwts = ones(201,1);
quadwts(2:2:200) = 4;
quadwts(3:2:199) = 2;
quadwts = (4/200).*quadwts/3;
quadvals = [quadpts, quadwts];

basisobj1 = putquadvals(basisobj1, quadvals);

for ivalue=1:2
    basisvalues1    = eval_basis(quadpts, basisobj1, ivalue-1);
    values1{ivalue} = basisvalues1.*(sqrt(quadwts)*ones(1,nbasis1));
end

basisobj1 = putvalues(basisobj1, values1);

%  set up basis for fitting second set of data

nbasis2   = 5;
norder2   = 4;
basisobj2 = create_bspline_basis([0,4],nbasis2,norder2);
% nbasis2   = 3;
% basisobj2 = create_monomial_basis([0,4],nbasis2);
basismat2 = eval_basis(t, basisobj2);

Bmat2     = basismat2'*basismat2;
Dmat2     = basismat2'*y2;

basisobj2 = putquadvals(basisobj2, quadvals);

for ivalue=1:2
    basisvalues2    = eval_basis(quadpts, basisobj2, ivalue-1);
    values2{ivalue} = basisvalues2.*(sqrt(quadwts)*ones(1,nbasis2));
end

basisobj2 = putvalues(basisobj2, values2);

%  set up bases for weight functions for output 1

%  variable 1:

% bbasis11 = create_monomial_basis([0,4],2);
bbasis11 = create_constant_basis([0,4]);

bbasis11 = putquadvals(bbasis11, quadvals);

%  compute basis fn & deriv. values and load into BASISOBJ
%  note:  do not weight basis function values

clear bvalues11
for ivalue=1:1
    bbasisvalues11    = eval_basis(quadpts, bbasis11);
    bvalues11{ivalue} = bbasisvalues11;
end

bbasis11 = putvalues(bbasis11, bvalues11);

%  variable 2:

bbasis12 = create_constant_basis([0,4]);

bbasis12 = putquadvals(bbasis12, quadvals);

%  compute basis fn & deriv. values and load into BASISOBJ
%  note:  do not weight basis function values

clear bvalues12
for ivalue=1:1
    bbasisvalues12    = eval_basis(quadpts, bbasis12);
    bvalues12{ivalue} = bbasisvalues12;
end

bbasis12 = putvalues(bbasis12, bvalues12);

% bwtcell1{1} = fdPar(fd([1;0], bbasis11));
bwtcell1{1} = fdPar(fd(1, bbasis12));
bwtcell1{2} = fdPar(fd(0, bbasis12));

%  set up bases for weight functions for output 2

%  variable 1:

bbasis21 = create_constant_basis([0,4]);

bbasis21 = putquadvals(bbasis21, quadvals);

%  compute basis fn & deriv. values and load into BASISOBJ
%  note:  do not weight basis function values

clear bvalues21
for ivalue=1:1
    bbasisvalues21    = eval_basis(quadpts, bbasis21);
    bvalues21{ivalue} = bbasisvalues21;
end

bbasis21 = putvalues(bbasis21, bvalues21);

%  variable 2:

% bbasis22 = create_monomial_basis([0,4],2);
bbasis22 = create_constant_basis([0,4]);

bbasis22 = putquadvals(bbasis22, quadvals);

%  compute basis fn & deriv. values and load into BASISOBJ
%  note:  do not weight basis function values

clear bvalues22
for ivalue=1:1
    bbasisvalues22    = eval_basis(quadpts, bbasis22);
    bvalues22{ivalue} = bbasisvalues22;
end

bbasis22 = putvalues(bbasis22, bvalues22);

bwtcell2{1} = fdPar(fd(1, bbasis21));
% bwtcell2{2} = fdPar(fd([1;0], bbasis22));
bwtcell2{2} = fdPar(fd(2, bbasis22));

%  set up basis for weight functions for input

%  Variable 1:

abasis1 = create_constant_basis([0,4]);

abasis1 = putquadvals(abasis1, quadvals);

%  compute basis fn & deriv. values and load into BASISOBJ
%  note:  these are not weighted

clear avalues1
abasisvalues1    = eval_basis(quadpts, abasis1);
avalues1{ivalue} = abasisvalues1;

abasis1 = putvalues(abasis1, avalues1);

awtcell1{1} = fdPar(fd(1, abasis1));

%  Variable 2:

% abasis2 = create_monomial_basis([0,4],2);
abasis2 = create_constant_basis([0,4]);

abasis2 = putquadvals(abasis2, quadvals);

%  compute basis fn & deriv. values and load into BASISOBJ
%  note:  these are not weighted

clear avalues2
abasisvalues2    = eval_basis(quadpts, abasis2);
avalues2{ivalue} = abasisvalues2;

abasis2 = putvalues(abasis2, avalues2);

% awtcell2{1} = fdPar(fd([2;0], abasis2));
awtcell2{1} = fdPar(fd(2, abasis2));

%  set up forcing function

ufdcell{1} = fd(1,abasis1);

clear fitcell

%  first fitcell

fitstruct1.y = y1;
fitstruct1.basisobj = basisobj1;
fitstruct1.basismat = basismat1;
fitstruct1.Bmat = Bmat1;
fitstruct1.Dmat = Dmat1;
fitstruct1.lambda = 1e0;
fitstruct1.bwtcell = bwtcell1;
% fitstruct1.awtcell = {};
% fitstruct1.ufdcell = {};
fitstruct1.awtcell = awtcell1;
fitstruct1.ufdcell = ufdcell;
fitcell{1} = fitstruct1;

%  second fitcell

fitstruct2.y = y2;
fitstruct2.basisobj = basisobj2;
fitstruct2.basismat = basismat2;
fitstruct2.Bmat = Bmat2;
fitstruct2.Dmat = Dmat2;
fitstruct2.lambda = 1e0;
fitstruct2.bwtcell = bwtcell2;
% fitstruct2.awtcell = {};
% fitstruct2.ufdcell = {};
fitstruct2.awtcell = awtcell2;
fitstruct2.ufdcell = ufdcell;
fitcell{2} = fitstruct2;

derivs = 1;

bvec = [2, 0, 0, 1]';  %  initial values for unforced example

bvec = [ [1, 0, 0],  1, ...
         [0, 2, 0], [2, 0] ]';  %  initial values for forced example

bvec = [ [2, 0],  2, ...
         [0, 1],  1 ]';         %  initial values for forced example

[PENSSE, DPENSSE, coefcell] = ...
                 profPDAm(bvec, fitcell, derivs);

coefcell{1}

coefcell{2}

fdobj1 = fd(coefcell{1}, basisobj1);
fdobj2 = fd(coefcell{2}, basisobj2);

subplot(2,1,1)
plotfit_fd(y1, t, fdobj1)
subplot(2,1,2)
plotfit_fd(y2, t, fdobj2)

PENSSE

DPENSSE

%  check out derivatives

PENSSE0  = PENSSE;
bvec0 = bvec;

for i=1:length(bvec);
    bvec = bvec0;
    bvec(i) = bvec(i) + 0.00001;
    PENSSE = profPDAm(bvec, fitcell, derivs);
    [(PENSSE - PENSSE0)/0.00001, DPENSSE(i)]
end

%  set up options for FMINUNC

options = optimset('LargeScale', 'off', 'Display', ...
                   'iter', 'MaxIter', 20, 'GradObj', 'off', ...
                   'TolFun', 1e-5, 'TolCon', 1e-5, ...
                   'TolX',   1e-5, 'TolPCG', 1e-5);

% optimize the fit

bvec0 = [2, 0, 0, 1]';  %  initial values for unforced example
bvec0 = [ [2, 0.1],  2, [0.1, 1],  1 ]'; %  initial values for forced example


[bvec, fval, exitflag, output, grad] = ...
            fminunc(@profPDAm, bvec0, options, ...
                    fitcell, 0);

bvec

%     1.7992
%     0.2370
%     0.1322
%     0.9671

[PENSSE, DPENSSE, coefcell, penmatcell, Dpenmatcell] = ...
                 profPDAm(bvec, fitcell, derivs);

fdobj1 = fd(coefcell{1}, basisobj1);
fdobj2 = fd(coefcell{2}, basisobj2);

subplot(2,1,1)
plotfit_fd(y1, t, fdobj1)
subplot(2,1,2)
plotfit_fd(y2, t, fdobj2)

PENSSE

DPENSSE

for ivar=1:2
    disp(ivar)
    Dpenstruct = Dpenmatcell{ivar}; 
    disp(Dpenstruct.DSmat);
    disp(Dpenstruct.DTmat);
    disp(Dpenstruct.DUmat);
    disp(Dpenstruct.DaVmat);
    disp(Dpenstruct.DbVmat);
end

