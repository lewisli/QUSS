Windows:

addpath ('c:\matlab\fdaM')
addpath ('c:\matlab\fdaM\examples\weather')

%  Last modified 6 September 2004

%  -----------------------------------------------------------------------
%                     Daily Weather Data
%  -----------------------------------------------------------------------

%  ------------------------  input the data  -----------------------

fid    = fopen('dailtemp.dat','rt');
tempav = fscanf(fid,'%f');
tempav = reshape(tempav, [365,35]);

fid    = fopen('dailprec.dat','rt');
precav = fscanf(fid,'%f');
precav = reshape(precav, [365,35]);

%  set up the times of observation at noon

daytime = (1:365)-0.5;

%  day values roughly in weeks

weeks = linspace(0,365,53)';   

%  define 8-character names for stations

place = [ ...
'arvida  '; 'bagottvi'; 'calgary '; 'charlott'; 'churchil'; 'dawson  '; ...
'edmonton'; 'frederic'; 'halifax '; 'inuvik  '; 'iqaluit '; 'kamloops'; ...
'london  '; 'montreal'; 'ottawa  '; 'princeal'; 'princege'; 'princeru'; ...
'quebec  '; 'regina  '; 'resolute'; 'scheffer'; 'sherbroo'; 'stjohns '; ...
'sydney  '; 'thepas  '; 'thunderb'; 'toronto '; 'uraniumc'; 'vancouvr'; ...
'victoria'; 'whitehor'; 'winnipeg'; 'yarmouth'; 'yellowkn'];

%  define 1-character names for months

monthletter = ['J'; 'F'; 'M'; 'A'; 'M'; 'J'; 'J'; 'A'; 'S'; 'O'; 'N'; 'D'];

%  -------------  set up fourier basis  ---------------------------
%  Here it was decided that 65 basis functions captured enough of
%  the detail in the temperature data: about one basis function
%  per week.  However, see below for smoothing with a saturated
%  basis (365 basis functions) where smoothing is defined by the
%  GCV criterion.

nbasis     = 65;
daybasis65 = create_fourier_basis([0,365], nbasis);

%  ---------  create fd objects for temp. and prec. ---------------

daytempfd = data2fd(tempav, daytime, daybasis65);
daytempfd_fdnames{1} = 'Day';
daytempfd_fdnames{2} = 'Station';
daytempfd_fdnames{3} = 'Deg C';
daytempfd = putnames(daytempfd, daytempfd_fdnames);

dayprecfd2 = data2fd(precav, daytime, daybasis65);
dayprecfd_fdnames{1} = 'Day';
dayprecfd_fdnames{2} = 'Station';
dayprecfd_fdnames{3} = 'mm';
dayprecfd = putnames(dayprecfd, dayprecfd_fdnames);

%  Plot temperature curves and values

plotfit_fd(tempav, daytime, daytempfd, place)

%  Plot residuals for three best fits and three worst fits

casenames = place;
varnames  = 'Temperature';
rng       = [0,365];
index     = [1,2,3,33,34,35];
residual  = 1;
sortwrd   = 1;

plotfit_fd(tempav, daytime, daytempfd, casenames, varnames, ...
           residual, sortwrd, rng, index)

%  Plot precipitation curves and values

plotfit_fd(precav, daytime, dayprecfd, place)

%  Assessment: the functions are definitely too rough with
%  this many basis functions, and especially for precip. which
%  has a much higher noise level.

%  these smoothing parameter values probably undersmooth the data,
%  but we can impose further smoothness on the results of analyses

%  set up the harmonic acceleration operator

Lbasis  = create_constant_basis([0,365]);  %  create a constant basis
Lcoef   = [0,(2*pi/365)^2,0];    %  set up three coefficients
wfd     = fd(Lcoef,Lbasis);      % define an FD object for weight functions
wfdcell = fd2cell(wfd);          % convert the FD object to a cell object
harmaccelLfd = Lfd(3, wfdcell);  %  define the operator object

%  set up the functional parameter objects to define smoothing

templambda = 1e5;
preclambda = 1e5;

tempfdPar  = fdPar(daybasis65, harmaccelLfd, templambda);
precfdPar  = fdPar(daybasis65, harmaccelLfd, preclambda);

daytempfd = smooth_fd(daytempfd, tempfdPar);
dayprecfd = smooth_fd(dayprecfd2, precfdPar);

%  plot each pair of functions along with raw data

tempmat = eval_fd(daytempfd, daytime);
precmat = eval_fd(dayprecfd, daytime);

index = 1:35
for i = index
  subplot(2,1,1)
  plot(daytime,tempav(:,i),'bo',daytime,tempmat(:,i),'r-','LineWidth',2)
  axis([0 365 -5 20])
  xlabel('Day','fontsize',20)
  ylabel('Temperature (deg. C)','fontsize',14)
  title(place(i,:),'fontsize',20)
  subplot(2,1,2)
  plot(daytime,precav(:,i),'bo',daytime,precmat(:,i),'r-','LineWidth',2)
  axis([0 365 0 10])
  xlabel('Day','fontsize',20)
  ylabel('Precipitation (mm)','fontsize',14)
  pause
end

%  plot all the functions

subplot(1,1,1)
plot(daytempfd);
  axis([0 365 -35 25])
xlabel('\fontsize{12} Day')
title('\fontsize{16} Mean Temperature')

plot(dayprecfd);
  axis([0 365 0 13])
xlabel('\fontsize{12} Day')
title('\fontsize{16} Mean Precipitation')

%  -------------------------------------------------------------
%                 Choose level of smoothing using 
%          the generalized cross-validation criterion
%              with smoothing function smooth_basis.
%  -------------------------------------------------------------

wtvec = ones(365,1);

% set up a saturated basis capable of interpolating the data

nbasis      = 365;  
daybasis365 = create_fourier_basis([0,365], nbasis);

%  --------------------  smooth temperature  ------------------

%  set up range of smoothing parameters in log_10 units

loglam = (-5:1)';
nlam   = length(loglam);

dfsave  = zeros(nlam,1);
gcvsave = zeros(nlam,1);

%  loop through smoothing parameters

for ilam=1:length(loglam)
    lambda = 10^loglam(ilam)
    fdParobj = fdPar(daybasis365, harmaccelLfd, lambda);
    [fdobj, df, gcv] = smooth_basis(daytime, tempav, fdParobj);
    dfsave(ilam)  = df;
    gcvsave(ilam) = gcv;
end

%  display and plot degrees of freedom and GCV criterion

[loglam, dfsave, gcvsave]

subplot(2,1,1)
plot(loglam, gcvsave, 'o-')
ylabel('\fontsize{16} GCV Criterion')
title('\fontsize{16} Temperature Smoothing')
subplot(2,1,2)
plot(loglam, dfsave, 'o-')
xlabel('\fontsize{16} log_{10} \lambda')
ylabel('\fontsize{16} Degrees of Freedom')

%  Do final smooth with minimum GCV value

lambda   = 0.01;  %  minimum GCV estimate, corresponding to 255 df
fdParobj = fdPar(daybasis365, harmaccelLfd, lambda);

[daytempfd, df, gcv, coef, SSE] = ...
             smooth_basis(daytime, tempav, fdParobj);

%  estimate standard error of fit

stderr = sqrt(SSE/(35*(365-df)));  %  0.26 deg C

%  plot data and fit

subplot(1,1,1)
plotfit_fd(tempav, daytime, daytempfd, place)

%  --------------------  smooth precipitation  ------------------

%  set up range of smoothing parameters in log_10 units

loglam = (4:9)';
nlam = length(loglam);

dfsave  = zeros(nlam,1);
gcvsave = zeros(nlam,1);

%  loop through smoothing parameters
for ilam=1:nlam
    lambda = 10^loglam(ilam)
    fdParobj = fdPar(daybasis365, harmaccelLfd, lambda);
    [fdobj, df, gcv] = smooth_basis(daytime, precav, fdParobj);
    dfsave(ilam)  = df;
    gcvsave(ilam) = gcv;
end

%  display and plot degrees of freedom and GCV criterion

[loglam, dfsave, gcvsave]

subplot(2,1,1)
plot(loglam, gcvsave, 'o-')
ylabel('\fontsize{16} GCV Criterion')
title('\fontsize{16} Precipitation Smoothing')
subplot(2,1,2)
plot(loglam, dfsave, 'o-')
xlabel('\fontsize{16} log_{10} \lambda')
ylabel('\fontsize{16} Degrees of Freedom')

%  Do final smooth with minimum GCV value

lambda = 1e7;  %  minimum GCV estimate, corresponding to 8.5 df
fdParobj = fdPar(daybasis365, harmaccelLfd, lambda);

[dayprecfd, df, gcv, coef, SSE] = ...
            smooth_basis(daytime, precav, fdParobj);

%  estimate standard error of fit

stderr = sqrt(SSE/(35*(365-df)));  %  0.94 mm

%  plot data and fit

plotfit_fd(precav, daytime, dayprecfd, place)

%  Assessment: the temperature curves are still pretty rough,
%  although the data themselves show that there are very
%  high frequency effects in the mean temperature, especially
%  early in the year. 
%  The precip. curves may be oversmoothed for some weather
%  stations. 

%  smooth precipitation in Prince Rupert

PRprecfd = smooth_basis(daytime, precav(:,18), fdParobj);

PRprecvec = eval_fd(daytime, PRprecfd);

subplot(1,1,1)

plot(daytime, precav(:,18), '.', ...
     daytime, PRprecvec,    '-')
ylabel('\fontsize{19} Precipitation (mm)')
axis([0,365,0,18])

%  -----------------------------------------------------------------
%                     PCA of temperature  
%  -----------------------------------------------------------------

nharm  = 4;
Lfd    = harmaccelLfd;
lambda = 1e4;

daytemppcastr = pca(daytempfd, nharm, lambda, Lfd);
daytemppcastr = varmx_pca(daytemppcastr);

%  plot harmonics

subplot(1,1,1)
plot_pca(daytemppcastr)

%  plot log eigenvalues

daytempharmeigval = daytemppcastr.eigvals;
x = ones(16,2);
x(:,2) = reshape((5:20),[16,1]);
y = log10(daytempharmeigval(5:20));
c = x\y;
subplot(1,1,1)
plot(1:20,log10(daytempharmeigval(1:20)),'-o', ...
     1:20, c(1)+ c(2).*(1:20), ':')
xlabel('Eigenvalue Number')
ylabel('Log10 Eigenvalue')

%  plot factor scores

harmscr = daytemppcastr.harmscr;

plot(harmscr(:,1), harmscr(:,2), 'o')
xlabel('Harmonic I')
ylabel('Harmonic II')
text(harmscr(:,1), harmscr(:,2), place)

%  ------------------------------------------------------------------
%               Functional linear models 
%  ------------------------------------------------------------------

%  ---------------------------------------------------------------
%             Predicting temperature from climate zone 
%  ---------------------------------------------------------------

%  set up a smaller basis using only 65 Fourier basis functions
%  to save some computation time

smallnbasis = 65;
smallbasis  = create_fourier_basis([0,365], smallnbasis);
tempfd      = data2fd(tempav, daytime, smallbasis);

smallbasismat = eval_basis(daytime, smallbasis);
y2cMap = inv(smallbasismat'*smallbasismat)*smallbasismat';

%  names for climate zones

zonenames = [ ...
'Canada  '; 'Atlantic'; 'Pacific '; 'Contintl'; 'Arctic  '];

%  indices for weather stations in each of four climate zones

atlindex = [1,2,4,8,9,13,14,15,19,22,23,24,25,28,34];
pacindex = [12,17,18,30,31];
conindex = [3,5,6,7,16,20,26,27,29,32,33,35];
artindex = [10,11,21];

%  Set up a design matrix having a column for the grand mean, and
%    a column for each climate zone effect. Add a dummy contraint
%    observation

zmat = zeros(35,5);
zmat(:       ,1) = 1;
zmat(atlindex,2) = 1;
zmat(pacindex,3) = 1;
zmat(conindex,4) = 1;
zmat(artindex,5) = 1;

%  attach a row of 0, 1, 1, 1, 1 to force zone
%  effects to sum to zero, and define first regression
%  function as grand mean for all stations

z36    = ones(1,5);
z36(1) = 0;
zmat   = [zmat; z36];

%  revise YFDOBJ by adding a zero function

coef   = getcoef(tempfd);  
coef36 = [coef,zeros(smallnbasis,1)];  
tempfd = putcoef(tempfd, coef36);  

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
betafdPar = fdPar(betafd, harmaccelLfd, lambda, estimate);
clear betacell
for j=1:p
    betacell{j} = betafdPar;
end

%  compute regression coefficient functions and 
%  predicted functions

[betaestcell, yhatfdobj] = fRegress(tempfd, xfdcell, betacell);

%  plot regression functions

for j=1:p
    plot(getfd(betaestcell{j}))
    title(['\fontsize{16} ',zonenames(j,:)])
    pause
end

%  plot predicted functions

plot(yhatfdobj)

%  compute residual matrix and get covariance of residuals

yhatmat  = eval_fd(daytime, yhatfdobj);
ymat     = eval_fd(daytime, tempfd);
temprmat = ymat(:,1:35) - yhatmat(:,1:35);
SigmaE   = cov(temprmat');

%  plot covariance surface for errors

contour(SigmaE)
colorbar

%  Repeat regression, this time outputting results for
%  confidence intervals

[betaestcell, yhatfdobj, betastderrcell, bvar, c2bmap] = ...
    fRegress(tempfd, xfdcell, betacell, y2cMap, SigmaE);

%  plot regression functions with confidence limits

for j=1:p
    plotbeta(betaestcell{j}, betastderrcell{j}, weeks)
    title(['\fontsize{16} ',zonenames(j,:)])
    pause
end

%  -----------------------------------------------------------------------
%         predict log precipitation from climate zone and temperature
%  -----------------------------------------------------------------------

%  set up functional data object for log precipitation

logprecmat = log10(eval_fd(daytime,dayprecfd));

lnprecfd = data2fd(logprecmat, daytime, smallbasis);

lnprecfd_fdnames{1} = 'Months';
lnprecfd_fdnames{2} = 'Station';
lnprecfd_fdnames{3} = 'log_{10} mm';
lnprecfd = putnames(lnprecfd, lnprecfd_fdnames);

%  plot log precipitation functions

plot(lnprecfd);
title('Log Precipitation Functions')

%  revise LOGPREDFD by adding a zero function

coef   = getcoef(lnprecfd);  
nbasis = getnbasis(smallbasis);
coef36 = [coef,zeros(nbasis,1)];  
lnprecfd = putcoef(lnprecfd, coef36);  

p = 5;
clear xfdcell
for j=1:p
    xfdcell{j} = zmat(:,j);
end

%  set up a FD object for temperature residuals

lambda     = 1e5;
fdParobj   = fdPar(smallbasis, harmaccelLfd, lambda);
temprfdobj = smooth_basis(daytime, temprmat, fdParobj);

%  plot temperature residuals

plot(temprfdobj)

%  extend temperature residual functions to include
%  zero function

coef   = getcoef(temprfdobj); 
nbasis = size(coef,1);
coef36 = [coef,zeros(nbasis,1)];  
temprfdobj = putcoef(temprfdobj, coef36);  

%  add TEMPRFDOBJ to the set of predictors

xfdcell{6}  = temprfdobj;
betacell{6} = betafdPar;
p = 6;

%  set up the basis for the regression functions

nbetabasis = 13;
betabasis  = create_fourier_basis([0,365], nbetabasis);

%  set up the functional parameter object for the regression fns.

betafd    = fd(zeros(nbetabasis,p), betabasis);
estimate  = 1;
lambda    = 0;
betafdPar = fdPar(betafd, harmaccelLfd, lambda, estimate);
clear betacell
for j=1:p
    betacell{j} = betafdPar;
end

%  compute regression coefficient functions and 
%  predicted functions

[betaestcell, yhatfdobj] = fRegress(lnprecfd, xfdcell, betacell);

%  plot regression functions

prednames = [zonenames; 'tempres '];
for j=1:p
    plot(getfd(betaestcell{j}))
    title(['\fontsize{16} ',prednames(j,:)])
    pause
end

%  plot predicted functions

plot(yhatfdobj)

%  compute residual matrix and get covariance of residuals

yhatmat    = eval_fd(daytime, yhatfdobj);
ymat       = eval_fd(daytime, lnprecfd);
lnprecrmat = ymat(:,1:35) - yhatmat(:,1:35);
SigmaE  = cov(lnprecrmat');

contour(SigmaE)
colorbar

%  repeat regression analysis to get confidence intervals

[betaestcell, yhatfdobj, betastderrcell, bvar, c2bmap] = ...
    fRegress(lnprecfd, xfdcell, betacell, y2cMap, SigmaE);

%  plot regression functions

for j=1:p
    plotbeta(betaestcell{j}, betastderrcell{j}, weeks)
    title(['\fontsize{16} ',prednames(j,:)])
    pause
end

%  ---------------------------------------------------------------
%      log annual precipitation predicted by temperature profile
%  ---------------------------------------------------------------

%  set up log10 total precipitation 

annualprec = log10(sum(precav))';

%  set up a smaller basis using only 65 Fourier basis functions
%  to save some computation time

smallnbasis = 65;
smallbasis  = create_fourier_basis([0,365], smallnbasis);
tempfd      = data2fd(tempav, daytime, smallbasis);

smallbasismat = eval_basis(daytime, smallbasis);
y2cMap = inv(smallbasismat'*smallbasismat)*smallbasismat';

%  set up the covariates, the first the constant, and the second
%  temperature

p = 2;
constantfd = fd(ones(1,35), create_constant_basis([0,365]));
clear xfdcell
xfdcell{1} = constantfd;
xfdcell{2} = tempfd;

%  set up the functional parameter object for the regression fns.
%  the smoothing parameter for the temperature function
%  is obviously too small here, and will be revised below by
%  using cross-validation.

clear betacell
%  set up the first regression function as a constant
betabasis1 = create_constant_basis([0,365]);
betafd1    = fd(0, betabasis1);
betafdPar1 = fdPar(betafd1);
betacell{1} = betafdPar1;
%  set up the second with same basis as for temperature
%  35 basis functions would permit a perfect fit to the data
nbetabasis  = 35;
betabasis2  = create_fourier_basis([0,365], nbetabasis);
betafd2     = fd(zeros(nbetabasis,1), betabasis2);
lambda      = 1;
betafdPar2  = fdPar(betafd2, harmaccelLfd, lambda);
betacell{2} = betafdPar2;

%  carry out the regression analysis

[betaestcell, annualprechat] = ...
                   fRegress(annualprec, xfdcell, betacell);

%  constant term

getcoef(getfd(betaestcell{1}))

%  plot the coefficient function for temperature
               
plot(getfd(betaestcell{2}))
title('\fontsize{16} Regression coefficient for temperature')

%  plot the fit

plot(annualprechat, annualprec,    'o', ...
     annualprechat, annualprechat, '--')

%  compute cross-validated SSE's for a range of smoothing parameters

loglam = (5:0.5:15)';
nlam   = length(loglam);
SSE_CV = zeros(nlam,1);
for ilam = 1:nlam;
    lambda       = 10^loglam(ilam);
    betacelli    = betacell;
    betacelli{2} = putlambda(betacell{2}, lambda);
    SSE_CV(ilam) = fRegress_CV(annualprec, xfdcell, betacelli);
    fprintf('%3.f %6.2f %10.4f\n', ilam, loglam(ilam), SSE_CV(ilam));
end

plot(loglam, SSE_CV, 'bo-')
xlabel('\fontsize{19} log_{10} smoothing parameter \lambda')
ylabel('\fontsize{19} Cross-validation score')

%  analysis with minimum CV smoothing

lambda      = 10^12.5;
betafdPar2  = fdPar(betafd2, harmaccelLfd, lambda);
betacell{2} = betafdPar2;

%  carry out the regression analysis

[betaestcell, annualprechat] = ...
                   fRegress(annualprec, xfdcell, betacell);

%  constant term

getcoef(getfd(betaestcell{1}))

%  plot the coefficient function for temperature
               
plot(getfd(betaestcell{2}))
title('\fontsize{16} Regression coefficient for temperature')

%  plot the fit

plot(annualprechat, annualprec, 'o', ...
     annualprechat, annualprechat, '--')
xlabel('\fontsize{16} Predicted log precipitation')
ylabel('\fontsize{16} Observed log precipitation')

%  compute squared multiple correlation

covmat = cov([annualprec, annualprechat]);
Rsqrd = covmat(1,2)^2/(covmat(1,1)*covmat(2,2))
%   0.7540

%  compute sigmae

resid = annualprec - annualprechat;
sigmae = mean(resid.^2);

%  recompute the analysis to get confidence limits

[betaestcell, annualprechat, betastderrcell] = ...
                   fRegress(annualprec, xfdcell, betacell, [], sigmae);

%  constant  coefficient standard error:

getcoef(betastderrcell{1})

%  plot the temperature coefficient function
               
plotbeta(betaestcell{2}, betastderrcell{2})
title('\fontsize{16} Regression coefficient for temperature')

%  ---------------------------------------------------------------
%         predict log precipitation from temperature
%  ---------------------------------------------------------------

%  change 0's to 0.05 mm in precipitation data

prectmp = precav;
for j=1:35
    index = find(prectmp(:,j)==0);
    prectmp(index,j) = 0.05;
end

%  work with log base 10 precipitation

logprec = log10(prectmp);

%  set up functional data object for log precipitation

logprecmat = log10(eval_fd(daytime,dayprecfd));

lnprecfd = data2fd(logprecmat, daytime, smallbasis);
lnprecfd_fdnames{1} = 'Months';
lnprecfd_fdnames{2} = 'Station';
lnprecfd_fdnames{3} = 'log_{10} mm';
lnprecfd = putnames(lnprecfd, lnprecfd_fdnames);

%  plot precipitation functions

plot(lnprecfd);
title('Log Precipitation Functions')

%  set up smoothing levels for s (xLfd) and for t (yLfd)

xLfd = harmaccelLfd;
yLfd = harmaccelLfd;
xlambda = 1e9;
ylambda = 1e7;

%  compute the linear model

wtvec = ones(35,1);
linmodstr = linmod(daytempfd, lnprecfd, wtvec, ...
                   xLfd, yLfd, xlambda, ylambda);

afd = linmodstr.alpha;   %  The intercept function
bfd = linmodstr.reg;     %  The bivariate regression function

%  plot the intercept function

plot(afd);

%  plot the regression function as a surface

bfdmat = eval_bifd(bfd, weeks, weeks);

subplot(1,1,1)
surf(weeks, weeks, bfdmat)
xlabel('\fontsize{12} Day(t)')
ylabel('\fontsize{12} Day(s)')

%  Get fitted functions

lnprechatfd = linmodstr.yhat;

% Compute mean function as a benchmark for comparison

lnprecmeanfd = mean(lnprecfd);
lnprechat0 = eval_fd(weeks, lnprecmeanfd);

%  Plot actual observed, fitted, and mean log precipitation for
%      each weather station, 

for i=1:35
    lnpreci    = eval_fd(lnprecfd(i),    weeks);
    lnprechati = eval_fd(lnprechatfd(i), weeks);
    SSE = sum((lnpreci-lnprechati).^2);
    SSY = sum((lnpreci-lnprechat0).^2);
    RSQ = (SSY-SSE)/SSY;
    plot(weeks, lnpreci, 'o', weeks, lnprechati, '-')
    xlabel('\fontsize{12} Day')
    ylabel('\fontsize{12} Log Precipitation')
    title(['\fontsize{16}', place(i,:),'  R^2 = ',num2str(RSQ)])
    pause
end

%  -------------------------------------------------------------------
%              Smooth Vancouver's precipitation with a 
%                       positive function.
%  -------------------------------------------------------------------

%  select Vancouver's precipitation

VanPrec  = precav(:,30);

%  smooth the data using 65 basis functions

lambda    = 1e4;
fdParobj  = fdPar(smallbasis, harmaccelLfd, lambda);
VanPrecfd = smooth_basis(daytime, VanPrec, fdParobj);
                      
%  Plot temperature curves and values

plotfit_fd(VanPrec, daytime, VanPrecfd)

%  set up the functional parameter object for positive smoothing

dayfdPar = fdPar(smallbasis, harmaccelLfd, lambda);

%  smooth the data with a positive function

[Wfd1, Fstr, iternum, iterhist] = ...
   smooth_pos(daytime, VanPrec, dayfdPar);

%  plot both the original smooth and positive smooth

VanPrecvec    = eval_fd(daytime, VanPrecfd);
VanPrecposvec = eval_pos(daytime, Wfd1);

plot(daytime, VanPrec, '.', daytime, VanPrecposvec, 'b-', ...
     daytime, VanPrecvec, 'r--')
legend('Observed', 'Positive smooth', 'Unrestricted smooth')
 
%  plot the residuals

VanPrecres = VanPrec - VanPrecposvec;
plot(daytime, VanPrecres.^2, '.', [0,365], [0,0], ':')
title('Residuals from positive fit')

%  compute a postive smooth of the squared residuals

lambda = 1e3;
dayfdPar = fdPar(smallbasis, harmaccelLfd, lambda);
[Wfd, Fstr, iternum, iterhist] = ...
   smooth_pos(daytime, VanPrecres.^2, dayfdPar);

%  plot the square root of this smooth along with the residuals

VanPrecvarhat = eval_pos(daytime, Wfd);
VanPrecstdhat = sqrt(VanPrecvarhat);
plot(daytime, VanPrecres.^2, '.', daytime, VanPrecvarhat, 'b-', ...
     [0,365], [0,0], ':')

%  set up a weight function for revised smoothing

wtvec = 1./VanPrecvarhat;
 
lambda   = 1e3;
dayfdPar = fdPar(smallbasis, harmaccelLfd, lambda);
[Wfd2, Fstr, iternum, iterhist] = ...
         smooth_pos(daytime, VanPrec, dayfdPar, wtvec);

%  plot the two smooths, one with weighting, one without

VanPrecposvec2 = eval_pos(daytime, Wfd2);

plot(daytime, VanPrec, '.', daytime, VanPrecposvec2, 'b-', ...
     daytime, VanPrecposvec, 'r--')
legend('Observed', 'Weighted', 'Unweighted')

 