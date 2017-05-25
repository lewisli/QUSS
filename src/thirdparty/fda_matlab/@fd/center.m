function centerfd = center(fd)
%  CENTER Center functional observations by subtracting mean function
%  Returns CENTERFD, the centered functional data object

%  Last modified 1 July 1998

  if ~isa_fd(fd)
    error ('Argument is not a functional data object.');
  end

  coef   = getcoef(fd);
  coefd  = size(coef);
  ndim   = length(coefd);
  nbasis = coefd(1);
  nrep   = coefd(2);
  onebas = ones(1, nrep);
  basisobj = getbasis(fd);
  if ndim == 2
    coefmean = mean(coef')';
    coef = coef - coefmean * onebas;
  else
    nvar = coefd(3);
    for j = 1:nvar
      coefmean = mean(coef(:,:,j)')';
      coef(:,:,j) = coef(:,:,j) - coefmean * onebas;
    end
  end
  fdnames = getnames(fd);
  fdnames{3} = ['Centered ', fdnames{3}];

  centerfd.coef     = coef;
  centerfd.basisobj = basisobj;
  centerfd.fdnames  = fdnames;

  centerfd = class(centerfd, 'fd');

