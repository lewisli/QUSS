function cycleplot_fd(fd, matplt, nx)
%CYCLEPLOT Performs a cycle plot of a functional data object FD,
%   assuming that FD is a bivariate function...the first component
%   of which is the x-coordinate and the second the y-coordinate
%
%  If MATPLT is 1, matplot is used to plot all curves in
%     a single plot.
%  Otherwise, each curve is plotted separately, and the
%     next curve is plotted when any key is pressed.
%  NX is the number of sampling points to use (default 128)

%  last modified 16 January 2001

if nargin < 3
    nx = 128;
end
if nargin < 2
    matplt = 1;
end

coef  = getcoef(fd);
coefd = size(coef);
ndim  = length(coefd);
if(ndim < 3)
    error('Univariate functions cannot be cycle plotted');
end
nbasis = coefd(1);
ncurve = coefd(2);
nvar   = coefd(3);
basisfd = getbasis(fd);
if(nvar > 2)
    warning('Only first two functions used');
end
rangex   = getbasisrange(basisfd);
x        = linspace(rangex(1), rangex(2), nx);
fdmat    = eval_fd(fd, x);
fdnames  = getnames(fd);

fdmean   = zeros(nx,2);
fdmean(:,1) = mean(squeeze(fdmat(:,:,1))')';
fdmean(:,2) = mean(squeeze(fdmat(:,:,2))')';

if matplt
    plot(fdmat(:,:,1), fdmat(:,:,2), '-')
else
    pltbot = min(min(min(fdmat))); 
    plttop = max(max(max(fdmat))); 
    for icurve = 1:ncurve
        plot(fdmat(:,icurve,1), fdmat(:,icurve,2), '-', ...
            fdmean(:,1), fdmean(:,2), '--');
        axis([pltbot, plttop, pltbot, plttop]);
        title(['Curve', num2str(icurve)]);
        pause;
    end
end

