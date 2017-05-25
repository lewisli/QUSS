%  test fRegressFM and fRegressFFtt using weather data

%  Last modified 30 November 2003

addpath ('c:\matlab\fdaM')
addpath ('c:\matlab\fdaM\examples\weather')

%  -----------------------------------------------------------------------
%                     Daily Weather Data
%  -----------------------------------------------------------------------

%    load dailydata dailydata.mat available

load dailydata

%  ------------------------  input the data  -----------------------

fid = fopen('dailtemp.dat','rt');
tempav = fscanf(fid,'%f');
tempav = reshape(tempav, [365,35]);

fid = fopen('dailprec.dat','rt');
precav = fscanf(fid,'%f');
precav = reshape(precav, [365,35]);

daytime  = (1:365)-0.5;
weektime = linspace(0,365,53)';   %  day values roughly in weeks

place = [ ...
'arvida  '; 'bagottvi'; 'calgary '; 'charlott'; 'churchil'; 'dawson  '; ...
'edmonton'; 'frederic'; 'halifax '; 'inuvik  '; 'iqaluit '; 'kamloops'; ...
'london  '; 'montreal'; 'ottawa  '; 'princeal'; 'princege'; 'princeru'; ...
'quebec  '; 'regina  '; 'resolute'; 'scheffer'; 'sherbroo'; 'stjohns '; ...
'sydney  '; 'thepas  '; 'thunderb'; 'toronto '; 'uraniumc'; 'vancouvr'; ...
'victoria'; 'whitehor'; 'winnipeg'; 'yarmouth'; 'yellowkn'];

%  indices for weather stations in each of four climate zones

atlindex = [1,2,4,8,9,13,14,15,19,22,23,24,25,28,34];pacindex = [12,17,18,30,31];conindex = [3,5,6,7,16,20,26,27,29,32,33,35];artindex = [10,11,21];
load dailydata

%  set up the harmonic acceleration operator

Lbasis  = create_constant_basis([0,365]);  %  create a constant basis
Lcoef   = [0,(2*pi/365)^2,0];    %  set up three coefficients
wfd     = fd(Lcoef,Lbasis);      % define an FD object for weight functions
wfdcell = fd2cell(wfd);          % convert the FD object to a cell object
harmaccelLfd = Lfd(3, wfdcell);  %  define the operator object

%  -----------------------------------------------------------------------
%       predict log precipitation from temperature and Dtemperature
%  -----------------------------------------------------------------------

%  set up precipitation data for logging ... fill 0's with 0.05's

temp = precav;
for j=1:35
    index = find(temp(:,j)==0);
    temp(index,j) = 0.05;
end

logprec = log10(temp);

%  smooth log10 precipitation

nbasis   = 365;
daybasis = create_fourier_basis([0,365], nbasis);

%  Set up the weight vector ... this will be revised below

Ywtvec = ones(365,1);  

lambda = 5e6;  %  minimizes GCV

[logprecfd, df, gcv, coef, SSE, penmat, y2cMap] = ...
    smooth_basis(daytime, logprec, daybasis, Ywtvec, harmaccelLfd, lambda);

[df, gcv]
%   9.4612    1.4211

plotfit_fd(logprec, daytime, logprecfd)

%  set up logprecfd with a smaller basis

logprecbasis = create_fourier_basis([0,365], 11);
logprecfd = data2fd(logprec, daytime, logprecbasis);

basismat = eval_basis(daytime, logprecbasis);
y2cMap = inv(basismat'*basismat)*basismat';

%  smooth temperature

% lambda = 1e-2;  %  minimizes GCV

lambda = 1e2;
[tempfd, df, gcv, coef, SSE, penmat, y2cMap] = ...
       smooth_basis(daytime, tempav, daybasis, Ywtvec, harmaccelLfd, lambda);

[df, gcv]
%  254.8718    7.9204  %  minimizes GCV
%   52.7789   13.7022  %  for lambda = 1e2

plotfit_fd(tempav, daytime, tempfd)

%  set up temperature with a smaller basis

tempbasis = create_fourier_basis([0,365], 11);
tempfd    = data2fd(tempav, daytime, tempbasis);

Dtempfd = deriv(tempfd, 1);

onefd = fd(ones(1,35),create_constant_basis([0,365]));

%  set up covariate cells

p = 3;
clear xfdcell
xfdcell{1} = onefd;
xfdcell{2} = tempfd;    
xfdcell{3} = Dtempfd;

for j=1:p
    plot(xfdcell{j})
    pause
end

load dailydata

%    set up cells for regression functions  

nbasis    = 11;
betabasis = create_fourier_basis([0,365], nbasis);
betafd0   = fd(zeros(nbasis,1), betabasis);
estimate  = 1;
lambda    = 1e2;
betafdPar = fdPar(betafd0, estimate, lambda, harmaccelLfd);

clear betacell
betacell{1} = betafdPar;
betacell{2} = betafdPar;
betacell{3} = betafdPar;

%  fit the linear model

betaestcell = fRegress(logprecfd, xfdcell, betacell);

plotbeta(betaestcell)

beta1fd = getfd(betaestcell{1});
beta2fd = getfd(betaestcell{2});
beta3fd = getfd(betaestcell{3});

plot(beta1fd)
title('Constant')

plot(beta2fd)
title('Temperature')

plot(beta3fd)
title('Temperature Derivative')

for j=1:p
    betafdj = getfd(betaestcell{j});
    plot(betafdj)
    pause
end

beta1mat = eval_fd(daytime, getfd(betaestcell{1}))*ones(1,35);
beta2mat = eval_fd(daytime, getfd(betaestcell{2}))*ones(1,35);
beta3mat = eval_fd(daytime, getfd(betaestcell{3}))*ones(1,35);

x1mat = eval_fd(daytime, xfdcell{1});
x2mat = eval_fd(daytime, xfdcell{2});
x3mat = eval_fd(daytime, xfdcell{3});

logprechat = x1mat.*beta1mat + x2mat.*beta2mat + x3mat.*beta3mat;
% logprechat = x1mat.*beta1mat + x2mat.*beta2mat;
% logprechat = x1mat.*beta1mat;

%  compute error mean squared errors for each station

% logprechat = eval_fd(daytime, precfdhat);

MSE = mean((logprec - logprechat).^2)';
sum(MSE)

MSE123 = MSE;  % full model       sum = 3.6464
MSE12  = MSE;  % no Dtemp         sum = 3.8288
MSE1   = MSE;  % no temp or Dtemp sum = 5.6884

% weeks = linspace(0,365,53)';   %  day values roughly in weeks
logprecmat = eval_fd(daytime, logprecfd);
logprecmn  = mean(logprecmat,2);
for i=1:35
    plot(daytime, logprecmat(:,i), '-',  ...
         daytime, logprechat(:,i), '--', ...
         daytime, logprecmn,       ':');
    axis([0,365,-1.1,1.1])
    title([place(i,:),' RMSE = ',num2str(sqrt(MSE(i)))])
    pause
end

%  plot Vancouver's fit

i = 30;
plot(daytime, logprecmat(:,i), 'k-',  ...
     daytime, logprec(:,i),    'k.',  ...
     daytime, logprechat(:,i), 'k--');
xlabel('\fontsize{16} Day')
ylabel('\fontsize{16} Log_{10} Precipitation (mm)')
axis([0,365,-1.4,1.0])

print -dpsc2 'c:/MyFiles/fdabook1/revision/figs.dir/VanLogPrec.ps'

%  compute residual covariance matrix

rmat    = logprecmat - logprechat;
SigmaE  = cov(rmat');

contour(SigmaE)
colorbar

plot(daytime,sqrt(diag(SigmaE)))

% re-analyze data with full set of arguments

[betaestcell, betastderrcell, yhatfdobj, BVariance, c2bMap] = ...
    fRegress(logprecfd, xfdcell, betacell, y2cMap, SigmaE);

%  set up mapping from B-matrix to beta functions

weeks     = linspace(0,365,53)';   %  day values roughly in weeks
plotbeta(betaestcell, betastderrcell, weeks)

%  ---------------------------------------------------------------
%    Test fRegress by predicting temperature from climate
%  ---------------------------------------------------------------

%  set up temperature with a smaller basis

smallnbasis    = 65;
smalltempbasis = create_fourier_basis([0,365], smallnbasis);
yfdobj         = data2fd(tempav, daytime, smalltempbasis);

tempbasismat = eval_basis(daytime, smalltempbasis);
y2cMap       = inv(tempbasismat'*tempbasismat)*tempbasismat';

%  Set up a design matrix having a column for the grand mean, and
%    a column for each climate zone effect. Add a dummy contraint
%    observation
zmat = zeros(35,5);zmat(:       ,1) = 1;zmat(atlindex,2) = 1;zmat(pacindex,3) = 1;zmat(conindex,4) = 1;zmat(artindex,5) = 1;
%  attach a row of 0, 1, 1, 1, 1
z36    = ones(1,5);z36(1) = 0;zmat   = [zmat; z36];
%  revise YFDOBJ by adding a zero function
coef   = getcoef(yfdobj);  
coef36 = [coef,zeros(smallnbasis,1)];  yfdobj = putcoef(yfdobj, coef36);  
p = 5;
clear xfdcell
for j=1:p
    xfdcell{j} = zmat(:,j);
end

%  set up the basis for the regression functions

nbetabasis = 13;
betabasis  = create_fourier_basis([0,365], nbetabasis);

%  set up the functional parameter object for the regression fns.

betafd    = fd(zeros(nbetabasis,p), betabasis);
estimate  = 1;
lambda    = 0;
betafdPar = fdPar(betafd, estimate, lambda, harmaccelLfd);
clear betacell
for j=1:p
    betacell{j} = betafdPar;
end

yfdPar = fdPar(yfdobj, 1, 0, harmaccelLfd);

%  First call:

[betaestcell, betastderrcell, yhatfdobj] = ...
    fRegress(yfdPar, xfdcell, betacell);

plotbeta(betaestcell)

plot(yhatfdobj)

%  compute residual matrix and get covariance of residuals

yhatmat = eval_fd(daytime, yhatfdobj);
ymat    = eval_fd(daytime, yfdobj);
rmat    = ymat(:,1:35) - yhatmat(:,1:35);
sigmae  = cov(rmat');

contour(sigmae)
colorbar

%  explore variation in residuals

lambda = 1e5;
rfdobj = smooth_basis(daytime, rmat, smalltempbasis, ones(365,1), ...
                      harmaccelLfd, lambda);

rvarobj = var(rfdobj);

rvarmat = eval_bifd(rvarobj, weektime, weektime);

surf(rvarmat)

contour(rvarmat)
colorbar

nharm   = 3;
rpcastr = pca(rfdobj, nharm, lambda, harmaccelLfd);
rpcastr = varmx_pca(rpcastr);

%  plot harmonics

subplot(1,1,1)
plot_pca(rpcastr)

%  three factors account for 98% of the variance

%  plot log eigenvalues ... this favors four, but three
%  gives a sigmae closer to the sigmae from the raw residuals

rharmeigval = rpcastr.eigvals;
x = ones(17,2);
x(:,2) = reshape((4:20),[17,1]);
y = log10(rharmeigval(4:20));
c = x\y;
subplot(1,1,1)
plot(1:20,log10(rharmeigval(1:20)),'-o', ...
     1:20, c(1)+ c(2).*(1:20), ':')
xlabel('Eigenvalue Number')
ylabel('Log10 Eigenvalue')

%  construct sigmae from first four harmonics

harmmat = eval_fd(daytime, rpcastr.harmfd);
sigmae = harmmat*diag(rharmeigval(1:nharm))*harmmat';

%  do the same using weektime for plotting purposes

harmmat = eval_fd(weektime, rpcastr.harmfd);
sigmae = harmmat*diag(rharmeigval(1:nharm))*harmmat';
surf(sigmae)
plot(weektime, sqrt(diag(sigmae)), '-', ...
     weektime, sqrt(diag(rvarmat)), '--')
 
%  Second call: to get confidence intervals

[betaestcell, betastderrcell, yhatfdobj, ncoef, bvar, c2bmap] = ...
    fRegress(yfdPar, xfdcell, betacell, y2cMap, sigmae);

plotbeta(betaestcell, betastderrcell, weektime)


%  add rfdobj to the set of predictors

xfdcell{6} = rfdobj;
betacell{6} = betafdPar;

p = 6;

%  -----------------------------------------------------------------------
% predict log precipitation from climate zone, temperature and Dtemperature
%  -----------------------------------------------------------------------

prectmp = precav;
for j=1:35
    index = find(prectmp(:,j)==0);
    prectmp(index,j) = 0.05;
end

logprec = log10(prectmp);

%  smooth log10 precipitation

nbasis   = 365;
daybasis = create_fourier_basis([0,365], nbasis);

%  Set up the weight vector ... this will be revised below

Ywtvec = ones(365,1);  

lambda = 5e6;  %  minimizes GCV

[logprecfd, df, gcv, coef, SSE, penmat, y2cMap] = ...
    smooth_basis(daytime, logprec, daybasis, Ywtvec, harmaccelLfd, lambda);

[df, gcv]
%   9.4612    1.4211

plotfit_fd(logprec, daytime, logprecfd)

%  set up logprecfd with a smaller basis

ynbasis = 11;
logprecbasis = create_fourier_basis([0,365], ynbasis);
logprecfd = data2fd(logprec, daytime, logprecbasis);

basismat = eval_basis(daytime, logprecbasis);
y2cMap = inv(basismat'*basismat)*basismat';

yfdobj = logprecfd;

%  Set up a design matrix having a column for the grand mean, and
%    a column for each climate zone effect. Add a dummy contraint
%    observation

zmat = zeros(35,5);
zmat(:       ,1) = 1;
zmat(atlindex,2) = 1;
zmat(pacindex,3) = 1;
zmat(conindex,4) = 1;
zmat(artindex,5) = 1;

%  attach a row of 0, 1, 1, 1, 1

z36    = ones(1,5);
z36(1) = 0;
zmat   = [zmat; z36];

%  revise YFDOBJ by adding a zero function

coef   = getcoef(yfdobj);  
coef36 = [coef,zeros(ynbasis,1)];  
yfdobj = putcoef(yfdobj, coef36);  

coef   = getcoef(rfdobj);  
coef36 = [coef,zeros(size(coef,1),1)];  
rfdobj = putcoef(rfdobj, coef36);  

p = 5
clear xfdcell
for j=1:p
    xfdcell{j} = zmat(:,j);
end

%  set up the basis for the regression functions

nbetabasis = 13;
betabasis  = create_fourier_basis([0,365], nbetabasis);

%  set up the functional parameter object for the regression fns.

betafd    = fd(zeros(nbetabasis,p), betabasis);
estimate  = 1;
lambda    = 0;
betafdPar = fdPar(betafd, estimate, lambda, harmaccelLfd);
clear betacell
for j=1:p
    betacell{j} = betafdPar;
end

yfdPar = fdPar(yfdobj, 1, 0, harmaccelLfd);

%  First call:

[betaestcell, betastderrcell, yhatfdobj] = ...
    fRegress(yfdPar, xfdcell, betacell);

plotbeta(betaestcell)

%  Second call:

[betaestcell, betastderrcell, yhatfdobj, ncoef, bvar, c2bmap] = ...
    fRegress(yfdPar, xfdcell, betacell, y2cMap, sigmae);

plotbeta(betaestcell, betastderrcell, weektime)

plot(yhatfdobj)

%  compute residual matrix and get covariance of residuals

yhatmat = eval_fd(daytime, yhatfdobj);
ymat    = eval_fd(daytime, yfdobj);
rmat    = ymat(:,1:35) - yhatmat(:,1:35);
sigmae  = cov(rmat');

contour(sigmae)
colorbar

%  add rfdobj to the set of predictors

xfdcell{6}  = rfdobj;
betacell{6} = betafdPar;

p = 6;

%  add Drfdobj to the set of predictors

Drfdobj = deriv(rfdobj, 1);

xfdcell{7}  = Drfdobj;
betacell{7} = betafdPar;

p = 7;

%  -----------------------------------------------------------------------
%          model log10 total precipitation using temperature  
%  -----------------------------------------------------------------------

load dailydata

%  set up the harmonic acceleration operator

Lbasis  = create_constant_basis([0,365]);  %  create a constant basis
Lcoef   = [0,(2*pi/365)^2,0];    %  set up three coefficients
wfd     = fd(Lcoef,Lbasis);      % define an FD object for weight functions
wfdcell = fd2cell(wfd);          % convert the FD object to a cell object
harmaccelLfd = Lfd(3, wfdcell);  %  define the operator object

annualprec = log10(sum(precav))';

%  set up temperature with a smaller basis

tempbasis = create_fourier_basis([0,365], 65);
tempfd    = data2fd(tempav, daytime, tempbasis);

%  set up covariate cells

p = 2;
clear xfdcell
xfdcell{1} = ones(35,1);
xfdcell{2} = tempfd;    

%    set up cells for regression functions  

onefd = fd(ones(1,35),create_constant_basis([0,365]));

nbetabasis = 11;
betabasis  = create_fourier_basis([0,365], nbetabasis);
betafd0    = fd(zeros(nbetabasis,1), betabasis);
estimate   = 1;
lambda     = 3e3;  %  minimizer of CV
Lfdobj     = int2Lfd(2);
betafdPar  = fdPar(betafd0, estimate, lambda, Lfdobj);

clear betacell
betacell{1} = fdpar(onefd, estimate, 0, 0);
betacell{2} = betafdPar;

%  fit the linear model

[betaestcell, betastderrcell, yhatfd] = ...
    fRegress(annualprec, xfdcell, betacell);

plotbeta(betaestcell)

getcoef(getfd(betaestcell{1}))

annualprechat = getcoef(yhatfd)';

plot(annualprechat, annualprec, 'o', ...
     annualprechat, annualprechat, '.')

sigmae = cov(annualprec - annualprechat);
sqrt(sigmae)

[betaestcell, betastderrcell, yhatfd, bvar, c2bmap] = ...
    fRegress(annualprec, xfdcell, betacell, ones(35,1), 1, sigmae);

plotbeta(betaestcell, betastderrcell)

annualprechat = getcoef(yhatfd)';

plot(annualprechat, annualprec, 'o', ...
     annualprechat, annualprechat, '.')
for i=1:35
text(annualprechat(i)+.01, annualprec(i), num2str(i))
end

%  Use cross-validation to select smoothing parameter

tfine = linspace(0,365,501);
delta = tfine(2)-tfine(1);
p = 2;
clear xfdcelli
xfdcelli{1} = ones(34,1);

lambda    = 3e3;
betafdPar = fdPar(betafd0, estimate, lambda, Lfdobj);
betacell{2} = betafdPar;

hwait = waitbar(0,'Please wait...');

SSE = 0;
for i=1:35
    waitbar(i/35,hwait);
    indexi = find((1:35) ~= i);
    yveci = annualprec(indexi);
    xfdcelli{2} = tempfd(indexi);
    betaestcelli = fRegress(yveci, xfdcelli, betacell);
    beta1i   = getcoef(getfd(betaestcelli{1}));
    beta2fdi =         getfd(betaestcelli{2});
    zmati    = eval_fd(tfine,tempfd(i));
    beta2i   = eval_fd(tfine,beta2fdi);
    yhati    = beta1i + delta.*sum(zmati.*beta2i);
    SSE = SSE + (annualprec(i)-yhati)^2;
end
SSE

%  -------------------------------------------------------------
%                   Analyses of lipemg data 
%    In this version we lag EMG by 50 milliseconds, sacrificing
%     the first 50 milliseconds of the lip acceleration records
%  -------------------------------------------------------------

addpath ('c:\matlab\lipemg')

load lipemgdata

%  sample size

N = 32;

%  upper time limit

tn = 0.69;

%  -------------------------------------------------------------
%  set up FD objects for lip position and acceleration and emg curves
%  -------------------------------------------------------------

emgmat    = load ('EMG.dat');
lipaccmat = load ('Acc.dat');
lipposmat = load ('Pos.dat');

timevec = linspace(0,tn,501)';   %  sampling points for each curve

save lipemgdata

%  ----------  set up the b-spline basis object  ------------
%       use order 6 splines so we can look at acceleration

nbasis = 100;
norder =   6;
basis  = create_bspline_basis([0,tn], nbasis, norder);

%  -----  create the fd objects  ---------

emgfd = data2fd(emgmat, timevec, basis);
emgfd_fdnames{1} = 'Time';
emgfd_fdnames{2} = 'Replications';
emgfd_fdnames{3} = 'Millivolts';
emgfd = putnames(emgfd, emgfd_fdnames);

lipaccfd = data2fd(lipaccmat, timevec, basis);
lipaccfd_fdnames{1} = 'Time';
lipaccfd_fdnames{2} = 'Replications';
lipaccfd_fdnames{3} = 'Meters/sec/sec';
lipaccfd = putnames(lipaccfd, lipaccfd_fdnames);

lipposfd = data2fd(lipposmat, timevec, basis);
lipposfd_fdnames{1} = 'Time';
lipposfd_fdnames{2} = 'Replications';
lipposfd_fdnames{3} = 'Meters';
lipposfd = putnames(lipposfd, lipposfd_fdnames);

subplot(2,1,1)
plot(lipaccfd)
subplot(2,1,2)
plot(emgfd)

%  -------------------------------------------------------------
%                      Set up the lagged records
%  -------------------------------------------------------------

%  fine set of time values

tfine = (0:0.001:tn)';
nfine = length(tfine);

%  unlagged function values

emgmat    = eval_fd(tfine, emgfd);
lipaccmat = eval_fd(tfine, lipaccfd);
lipposmat = eval_fd(tfine, lipposfd);

%  lagged function values

nlag = 50;

emgmatlag    =    emgmat(1:(nfine-nlag),:);
lipaccmatlag = lipaccmat((nlag+1):nfine,:);
lipposmatlag = lipposmat((nlag+1):nfine,:);

%  lagged time values

nfinelag  = nfine - nlag;
tfinelag  = tfine(1:nfinelag);

%  set up spline basis for representing data

tnlag  = tn - tfine(nlag+1);
nbasis = 93;
norder =  6;
basis  = create_bspline_basis([0,tnlag], nbasis, norder);

%  create functional data objects

lipacclagfd = data2fd(lipaccmatlag, tfinelag, basis);
emglagfd    = data2fd(emgmatlag,    tfinelag, basis);

basismat = eval_basis(tfinelag, basis);
y2cmap   = inv(basismat'*basismat)*basismat';

%  set up the covariates

onebasis = create_constant_basis([0,tnlag]);
onefd    = fd(ones(1,N), onebasis);

xfdcell{1} = onefd;
xfdcell{2} = emglagfd;

p = 2;

%  set up the functional parameter object for the regression fns.

nbetabasis = 21;
norder     =  6;
betabasis  = create_bspline_basis([0,tnlag], nbetabasis, norder);

betafd    = fd(zeros(nbetabasis,1), betabasis);
estimate  = 1;
lambda    = 0;
Lfd       = int2Lfd(2);
betafdPar = fdPar(betafd, estimate, lambda, Lfd);

clear betacell
for j=1:p
    betacell{j} = betafdPar;
end

%  First call:

[betaestcell, betastderrcell, yhatfdobj] = ...
    fRegress(lipacclagfd, xfdcell, betacell);

%  Second call:

[betaestcell, betastderrcell, yhatfdobj, bvar, c2bmap] = ...
    fRegress(lipaccfd, xfdcell, betacell, y2cmap, sigmae);

figure(2)
plotbeta(betaestcell)

plotbeta(betaestcell, betastderrcell)

yhatmat = eval_fd(tfinelag, yhatfdobj);
ymat    = eval_fd(tfinelag, yfdobj);
rmat    = ymat - yhatmat;
sigmae  = cov(rmat');

contour(sigmae)
colorbar

StdErrE = sqrt(diag(sigmae));
CorrE = diag(1./StdErrE)*sigmae*diag(1./StdErrE);

contour(CorrE)
colorbar

MSE = mean(rmat.^2)';
sqrt(mean(MSE))          %   0.7364

%  ------------------------------------------------------------
%        long term changes in trend and seasonality
%               in the nondurable goods index
%  ------------------------------------------------------------

%  attach FDA functions

addpath 'c:/Matlab/fdaM'

%  ------------------------------------------------------------------------
%                    Set up the data for analysis
%  ------------------------------------------------------------------------

load goodsdata  %  use this command if goodsdata.mat file available

load goodslmfun2fun

%  This is input of the data 

fid = fopen('nondurprod.dat','rt');
temp = fscanf(fid,'%f');
temp = reshape(temp, [18, 81]);
tempmat = temp(2:13,:);
tempmat(12,81) = 0;
nondurables = reshape(tempmat, [12*81, 1]);
nondurables = nondurables(1:971);
ndur = 971;

%  for completeness, make dec 99 equal to dec 98, jan 00 equal to jan 99

nondurables = [nondurables; nondurables(961,:)];
nondurables = [nondurables; nondurables(962,:)];
ndur = 973;

%  set up time values

durtime = (0:(ndur-1))'./12 + 1919;

%  compute log nondurables

lognondur = log10(nondurables);

save goodsdata

%  select data for years >= startyear;

startyear = 1952;
index     = find(durtime >= startyear);
durtime   = durtime(index) - startyear;
lognondur = lognondur(index);
n         = length(durtime);

%  plot index

plot(durtime+startyear, lognondur, '-')

%  smooth log10 goods index

norder     = 4;
nbasis     = n + norder - 2;
goodsrange = [durtime(1),durtime(n)];
goodsbasis = create_bspline_basis(goodsrange, nbasis);

%  Set up the weight vector ... this will be revised below

Ywtvec = ones(n,1);  
Lfdobj = int2Lfd(2);

% try smoothing with a range of lambda's

loglam = -7:0.25:-2;
gcvsave = zeros(length(loglam),1);
for i=1:length(loglam)
    lambdai = 10^loglam(i);
    [lognondurfd, df, gcv] = ...
             smooth_basis(durtime, lognondur, goodsbasis, Ywtvec, ...
                  Lfdobj, lambdai);
    disp([loglam(i), df, log10(gcv)])
    gcvsave(i) = gcv;
end

%  plot GCV function

plot(loglam, gcvsave)

%  Do the final smooth

% lambda = 1e-4;  %  minimizes GCV, df = 338

lambda = 5e-6;  %  gets the mid-year shoulder much better

[lognondurfd, df, gcv, coef, SSE, penmat, y2cMap] = ...
             smooth_basis(durtime, lognondur, goodsbasis, Ywtvec, ...
                  Lfdobj, lambda);

%  plot fit to data in blocks

nplots = 8;
nyears = floor(n/nplots/12);

goodsplotfit(lognondurfd, lognondur, [], ...
                      nplots, nyears, durtime, startyear)
                  
%  get values of smooth

lognondursmth = eval_fd(durtime, lognondurfd);

%  plot fits for selected years

years = [1965, 1966, 1995, 1996];
lognsave = zeros(101,length(years));
Xmat = ones(101,2);  Xmat(:,2) = (-50:50)';
m = 0;
for year = years
    m = m + 1;
    index = find(durtime+startyear >= year & ...
                 durtime+startyear <= year+1);
    durtemp = linspace(year-startyear,year-startyear+1,101)';
    lognsave(:,m) = eval_fd(durtemp, lognondurfd);
    lineartrend = Xmat*(Xmat\lognsave(:,m));
    lognsave(:,m) = lognsave(:,m) - lineartrend;
end

plot((0:100)/100, lognsave(:,1), 'k--', ...
     (0:100)/100, lognsave(:,2), 'k--', ...
     (0:100)/100, lognsave(:,3), 'k-', ...
     (0:100)/100, lognsave(:,4), 'k-')
xlabel('\fontsize{16} Time in years')
ylabel('\fontsize{16} Periodic trend')
% legend('\fontsize{12} 1965', '1966', '1995', '1996')

%  this is a string for labelling axes by month
\fontsize{16} J    F   M    A    M    J    J    A    S   O   N    D

%  set up fdPar object for lognondurfd

lognondurfdPar = fdPar(lognondurfd, 1, 1e-6, Lfdobj);

%  set up covariate cells

p           = 11;    
width       = 1;
ltnbasis = (2000 - startyear)/width + 3;
% ltnbasis    = 5;
ltlambda    = 1e-2;
ltpernbasis = 7;
ltperlambda = 0;

[betaestcell, betastderrcell, yhatfd, ncoef] = ...
    goodsfitfn(p, width, ltnbasis, ltlambda,  ...
               ltpernbasis, ltperlambda, ...
               startyear, lognondurfdPar);

ncoef
res = lognondursmth - eval_fd(durtime, yhatfd);
sqrdres = res.^2;
sqrt(mean(sqrdres)*(n/(n-ncoef)))

%  plot regression functions

plotbeta(betaestcell)

%  plot fit

goodsplotfit(yhatfd, lognondursmth, betaestcell, ...
             nplots, nyears, durtime, startyear)

%  plot periodic trend

yhatvec = eval_fd(durtime, yhatfd);
beta1vec = eval_fd(durtime, getfd(betaestcell{1}));

plot(durtime+startyear, yhatvec - beta1vec, '-')

years = [1965, 1966, 1995, 1996];
lognsave = zeros(101,length(years));
yhatsave = lognsave;
Xmat = ones(101,2);  Xmat(:,2) = (-50:50)';
m = 0;
for year = years
    m = m + 1;
    index = find(durtime+startyear >= year & ...
                 durtime+startyear <= year+1);
    durtemp = linspace(year-startyear,year-startyear+1,101)';
    lognsave(:,m) = eval_fd(durtemp, lognondurfd);
    yhatsave(:,m) = eval_fd(durtemp, yhatfd);
    lineartrend = Xmat*(Xmat\lognsave(:,m));
    lognsave(:,m) = lognsave(:,m) - lineartrend;
    yhatsave(:,m) = yhatsave(:,m) - lineartrend;
end

plot((0:100)/100, lognsave(:,2), 'b--', ...
     (0:100)/100, yhatsave(:,2), 'b-.', ...
     (0:100)/100, lognsave(:,4), 'm-', ...
     (0:100)/100, yhatsave(:,4), 'm.')
xlabel('\fontsize{16} Time in years')
ylabel('\fontsize{16} Periodic trend')

lambda = 1e0;  %  minimizes GCV, df = 338
[sqrdresfd, df, gcv, coef, SSE, penmat, y2cMap] = ...
             smooth_basis(durtime, sqrdres, goodsbasis, Ywtvec, ...
                  Lfdobj, lambda);

plotfit_fd(sqrdres, durtime, sqrdresfd)

sigmae = diag(eval_fd(durtime, sqrdresfd));

%  not worth doing this ... confidence limits are so tight
%    that they can't be seen in the plots

[betaestcell, betastderrcell, yhatfd, bvar, c2bmap] = ...
    fRegress(lognondurfdPar, xfdcell, betacell, y2cMap, sigmae);

plotbeta(betaestcell, betastderrcell)


