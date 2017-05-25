Windows:

addpath ('c:\matlab\fdaM')
addpath ('c:\matlab\fdaM\examples\melanoma')

Unix:

addpath('/export/home/steve/ramsay/u1/fdaM')
addpath('/export/home/steve/ramsay/u1/fdaM/examples/melanoma')

%  Last modified  9 January 2001

%  -----------------------------------------------------------------------
%                      Melanoma Incidence data
%  -----------------------------------------------------------------------

%  input the data  

tempmat = load('melanoma.dat');

year  = tempmat(:,2);
mela  = tempmat(:,3);
nyear = length(year);
rng   = [min(year),max(year)];

xmat1 = [ones(37,1),year];
coef1 = xmat1\mela;
melahat1 = xmat1*coef1;
sse1 = sum((mela-melahat1).^2);

plot(year, mela, 'o', year, melahat1, '--')

xmat2 = [ones(37,1),year,sin(0.65*year)];
coef2 = xmat2\mela;
melahat2 = xmat2*coef2;
sse2 = sum((mela-melahat2).^2);

plot(year, mela, 'ko', year, melahat1, 'k--', year, melahat2, 'k-')
xlabel('\fontsize{16} Year')
ylabel('\fontsize{16} Cases of Melanoma per 100,000')

print -dps2 'c:/MyFiles/P651/melanoma.ps'

lnmela = log(mela);

xmat1 = [ones(37,1),year];
coef1 = xmat1\lnmela;
lnmelahat1 = xmat1*coef1;
sse1 = sum((lnmela-lnmelahat1).^2);

plot(year, lnmela, 'o', year, lnmelahat1, '--')

xmat2 = [ones(37,1),year,sin(0.65*year)];
coef2 = xmat2\lnmela;
lnmelahat2 = xmat2*coef2;
sse2 = sum((lnmela-lnmelahat2).^2);

plot(year, lnmela, 'ko', year, lnmelahat1, 'k--', year, lnmelahat2, 'k-')
xlabel('\fontsize{16} Year')
ylabel('\fontsize{16} Log Cases of Melanoma per 100,000')

print -dps2 'c:/MyFiles/P651/melanoma.ps'




%  -----------------------------------------------------------------------
%             smooth data using B-splines
%  -----------------------------------------------------------------------

%  set up the basis with a knot at every year

knots  = year';
nbasis = nyear + 4;
norder = 6;
basis  = create_bspline_basis(rng, nbasis, norder, knots);

%  smooth the data by penalizing the second derivative

Lfd    = 2;
lambda = 1;
melafdPar = fdPar(basis, Lfd, lambda);

melafd = smooth_basis(mela, year, melafdPar);

%  plot the data and the smooth

plotfit_fd(mela, year, melafd)

%  plot the residuals

plotfit_fd(mela, year, melafd, [], [], 1, 1)

%  set up operator to remove sinusoid plus trend

omega = 0.65;
Lbasis = create_constant_basis(rng);
Lcoef  = [0, 0, omega^2, 0];

Lfd       = fd(Lcoef,Lbasis);
lambda    = 1e-2;
melafdPar = fdPar(basis, Lfd, lambda);

%  smooth the data

melafd = smooth_basis(mela, year, melafdPar);
melafd   = melalist.fdobj;

%  plot the results

plotfit_fd(mela, year, melafd)


%  -----------------------------------------------------------------------
%             smooth data using seasonfit function
%  -----------------------------------------------------------------------

basisA = create_power_basis(rng, 2);
basisB = basisA;
nphi   = 2;
period = 2*pi/0.65;

seasonfitstr = seasonfit(mela, year, nphi, period, basisA, basisB);

amat     = seasonfitstr.amat;
bmat     = seasonfitstr.bmat;
coeffn   = seasonfitstr.coeffn;
alpha    = seasonfitstr.alpha;
seasonal = seasonfitstr.seasonal;
melahat  = seasonfitstr.yhat;

plot(year, mela, 'o', year, melahat, '-', year, alpha, '--')

yearfine = (1936:0.1:1972);

seasonevalstr = seasoneval(seasonfitstr, yearfine);
melafine = seasonevalstr.yhat;

plot(year, mela, 'o', yearfine, melafine, '-', year, alpha, '--')






