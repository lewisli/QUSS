function plot_pca_mult(pcastr, nx, pointplot, harm, expand, cycle)
%  PLOT_PCA_MULT  Plots the harmonics for a principal components analysis.
%  It differs from PLOT_PCA in plotting these as multiple panels without
%    prompting.  It only accepts univariate functions.
%  Arguments:
%  PCASTR    ... Struct object returned by PCA_FD.
%  NX        ... Number of argument values for plotting. Default = 101.
%  POINTPLOT ... If pointplot=1, then the harmonics are plotted with
%                +'s and -'s on either side of the mean function.
%                Otherwise lines are used.  Default is 1.
%  HARM      ... If harm = 0 (the default) then all the computed harmonics
%                are plotted.   Otherwise those in HARM are plotted.
%  EXPAND    ... If expand =0 then effect of +/- 2 standard deviations of
%                each harmonic are given.
%                Otherwise the factor expand is used.
%  CYCLE     ... If cycle=T and there are 2 variables then a cycle plot
%                will be drawn.  If the number of variables is anything else,
%                CYCLE will be ignored.

%  last modified 26 November 2001

  %  set up default argument values

  if nargin < 6
    cycle = 0;
  end
  if nargin < 5
    expand = 0;
  end
  if nargin < 4
    harm = 0;
  end
  if nargin < 3
    pointplot = 1;
  end
  if nargin < 2
    nx = 51;
  end

  harmfd  = pcastr.harmfd;
  basis   = getbasis(harmfd);
  rangex  = getbasisrange(basis);
  fdnames = getnames(harmfd);
  x       = linspace(rangex(1), rangex(2), nx);
  fdmat   = eval_fd(harmfd, x);
  meanmat = squeeze(eval_fd(pcastr.meanfd, x));
  dimfd   = size(fdmat);
  nharm   = dimfd(2);
  if harm == 0
    harm = (1:nharm);
  end
  onesharm = ones(1,nharm);
  if expand == 0
      fac = ones(nx,1)*sqrt(pcastr.eigvals(1:nharm))';
  else
      fac = expand;
  end 
  meanplus  = meanmat*onesharm + fac.*fdmat;
  meanminus = meanmat*onesharm - fac.*fdmat;
  plottop = max(max([meanplus;meanminus]));
  plotbot = min(min([meanplus;meanminus]));
  if length(dimfd) == 2
    %  plotting for univariate functions
    ijtable = zeros(nharm,2);
    for iharm = harm
      percentvar = round(100 * pcastr.varprop(iharm));
      if nharm < 4
          subplot(1,nharm,iharm)
      else
          subplot(2,2,iharm)
      end
      if plottop*plotbot < 0
          plot(x, meanmat, '-', [min(x),max(x)], [0,0], '--')
      else
          plot(x, meanmat, '-')
      end
      text(x-.2, meanplus(:,iharm),  '+')
      text(x-.2, meanminus(:,iharm), '-')
      axis([rangex(1), rangex(2), plotbot, plottop])
      if ~(nharm == 4), axis('square'); end
      title(['\fontsize{12} Component ', num2str(iharm), '  ', ...
              num2str(percentvar), '%'])
    end
  else
    error('Multivariate functions cannot be plotted.');
  end

