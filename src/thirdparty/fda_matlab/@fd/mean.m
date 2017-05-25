function meanfd = mean(fd)
%  MEAN  Compute mean function for functional observations
%  Argument:
%  FD ... a functional data object for the observations
%  Return:
%  MEANFD ... a functional data object for the mean

%  last modified 1 July 1998

  if ~isa_fd(fd)
    error ('Argument FD is not a functional data object.');
  end

  coef     = getcoef(fd);
  coefd    = size(coef);
  ndim     = length(coefd);
  nbasis   = coefd(1);
  basisobj = getbasis(fd);
  if ndim == 2
    coefmean = mean(coef')';
  else
    nvar = coefd(3);
    coefmean = zeros(coefd(1),1,nvar);
    for j = 1:nvar
      coefmean(:,1,j) = mean(coef(:,:,j)')';
    end
  end

  fdnames = getnames(fd);
  fdnames{2} = 'Mean';
  fdnames{3} = ['Mean ', fdnames{3}];

  meanfd.coef     = coefmean;
  meanfd.basisobj = basisobj;
  meanfd.fdnames  = fdnames;

  meanfd = class(meanfd, 'fd');

