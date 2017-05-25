addpath ('c:\matlab\fdaM')
addpath ('c:\matlab\fdaM\examples\handwrit')

%  Last modified 6 September 2004

%  -----------------------------------------------------------------------
%                 Registered Handwriting Data
%  -----------------------------------------------------------------------

%  Input the data.  These handwriting functions have already been
%  registered.

fid = fopen('fdareg.dat','rt');
fdarray = reshape(fscanf(fid,'%f'), [20,2,1401]);
fdarray = permute(fdarray,[3,1,2]);
fdarray = fdarray/1000;   %  convert unit to meters

%  set up time values and range

fdatime   = linspace(0, 2.3, 1401)';

fdarange  = [0, 2.3];

%  After some experimentation, 205 order 6 splines without smoothing
%    seemed to do an adequate job of representing the original data.
%  Order 6 was used to get a reasonable estimate of the third deriv.

fdabasisobj = create_bspline_basis([0, 2.3], 205, 6);

%  set up the functional data structure without smoothing

fdafd = data2fd(fdarray, fdatime, fdabasisobj);
fdafd_fdnames{1} = 'Seconds';
fdafd_fdnames{2} = 'Replications';
fdafd_fdnames{3} = 'mm';
fdafd = putnames(fdafd, fdafd_fdnames);

%  plot all curves

plot(fdafd)

%  compute values of curves and the values of the curve

fdamat     = eval_fd(fdafd, fdatime);
fdameanvec = squeeze(eval_fd(mean(fdafd), fdatime));

%  plot individual curves, including both sampling points and fit
%  also plot the mean curve in the background

subplot(1,1,1);
for i = 1:20
  plot(fdarray(:,i,1),  fdarray(:,i,2), 'go', ...
       fdamat(:,i,1),   fdamat(:,i,2), '-', ...
       fdameanvec(:,1), fdameanvec(:,2), 'r--');
  axis([-.040, .040,-.040, .040]);
  title(['Record ', num2str(i)]);
  pause;
end

%  apply a small amount of smoothing to get better acceleration
%  curves

%  set up the functional parameter object to define amount of
%  smoothing

Lfdobj   = int2Lfd(4);
lambda   = 1e-12;
fdParobj = fdPar(fdafd, Lfdobj, lambda);

%  smooth the functional data object

fdasmthfd = smooth_fd(fdafd, fdParobj);

%  plot the acceleration records

subplot(1,1,1);
for i=1:20
    accel = squeeze(eval_fd(fdatime, fdasmthfd(i,:), int2Lfd(2)));
    plot(fdatime, accel(:,1), '-', fdatime, accel(:,2), '-', ...
         [0,2.3],[0,0], 'r:');
    axis([0, 2.3, -8, 8]);
    xlabel('\fontsize{16} Seconds')
    ylabel('\fontsize{16} meters/sec/sec')
    title(['\fontsize{16} Curve ',num2str(i)])
    legend('X', 'Y')
    pause
end

%  compute the acceleration magnitudes

D2fdamat = eval_fd(fdasmthfd, fdatime, int2Lfd(2));
D2mag = sqrt(D2fdamat(:,:,1).^2 + D2fdamat(:,:,2).^2);

%  plot the acceleration magnitudes

subplot(1,1,1);
plot(fdatime, D2mag)
axis([0,2.3,0,15000]);
xlabel('Seconds')
ylabel('mm/sec/sec')
title('Acceleration Magnitude');
axis([0,2.3,0,10]);

% ---------------  do a PCA of handwriting data  -------------------------------

%  do the PCA with varimax rotation

nharm   = 4;
lambda  = 1e-9;
gaitpca = pca(fdasmthfd, nharm, lambda, Lfdobj);
gaitpca = varmx_pca(gaitpca);

%  plot harmonics using cycle plots

subplot(1,1,1)
plot_pca(gaitpca, 101, 1, 0, 0, 1);

%  plot eigenvalues

gaiteigvals = gaitpca.eigvals;
x = ones(16,2);
x(:,2) = reshape((5:20),[16,1]);
y = log10(gaiteigvals(5:20));
c = x\y;
subplot(1,1,1)
plot(1:20,log10(gaiteigvals(1:20)),'-o', ...
     1:20, c(1)+ c(2).*(1:20), ':')
xlabel('Eigenvalue Number')
ylabel('Log10 Eigenvalue')

gaitharmfd = gaitpca.harmfd;
gaitharmmat = eval_fd(gaitharmfd,gaittime);
gaitvarprop = gaitpca.varprop;
gaitmeanvec = squeeze(eval_fd(gaitmeanfd,gaittime));

con = 5.*ones(1,4);
for j=1:4
    subplot(2,2,j)
    yplus = gaitmeanvec + con(j).*squeeze(gaitharmmat(:,j,:));
    plot(gaitmeanvec(:,1),gaitmeanvec(:,2),'g.')
    hold on
    for i=1:20
        plot([gaitmeanvec(i,1),yplus(i,1)],...
             [gaitmeanvec(i,2),yplus(i,2)],'b-')
    end
    hold off
    xlabel('Hip Angle')
    ylabel('Knee Angle')
    title(['PC ',num2str(j),' (',num2str(round(gaitvarprop(j)*1000)/10),'%)'])
    axis([-20,60,0,80])
end
