addpath ('c:\matlab\fdaM')
addpath ('c:\matlab\fdaM\examples\pinch')

%  Last modified 15 January 2003

%  -----------------------------------------------------------------------
%                Pinch force data
%  -----------------------------------------------------------------------

%  ------------------  input the data  --------------------

fid = fopen('pinch.dat','rt');
pinchvec = fscanf(fid,'%f');
pinchmat = reshape(pinchvec, [20,151])';

pinchtime  = linspace(0,150,151);

pinchbasis = create_bspline_basis([0,150], 42, 6);

%  -----------  create fd object (no smoothing)  --------------------

pinchfd = data2fd(pinchmat, pinchtime, pinchbasis);
pinchfd_fdnames{1} = 'Point index';
pinchfd_fdnames{2} = 'Replications';
pinchfd_fdnames{3} = 'Force (N)';
pinchfd = putnames(pinchfd, pinchfd_fdnames);

%  plot all curves

subplot(1,1,1)
plot(pinchfd)
title('Pinch Force Curves')

%  plot each curve along with the data

plotfit_fd(pinchmat, pinchtime, pinchfd)

%  plot the residuals, with cases sorted by size of mean squared residuals

casenames = [];
varnames  = [];
residual  = 1;
sortwrd   = 1;

plotfit_fd(pinchmat, pinchtime, pinchfd, casenames, varnames, ...
           residual, sortwrd)

%  ---------------------  do a PCA (light smoothing)  --------------

nharm  = 2;
lambda = 1e-2;

pinchpcastr = pca(pinchfd, nharm, lambda, int2Lfd(2));

plot_pca(pinchpcastr)

pincheigvals = pinchpcastr.eigvals;
x = ones(17,2);
x(:,2) = reshape((3:19),[17,1]);
y = log10(pincheigvals(3:19));
c = x\y;
subplot(1,1,1)
plot(1:19,log10(pincheigvals(1:19)),'-o', ...
     1:19, c(1)+ c(2).*(1:19), ':')
xlabel('Eigenvalue Number')
ylabel('Log10 Eigenvalue')

