function pcastr = pcaextend(fd, pcabasis, nharm, lambda, Lfd, centerfns)
%  PCA Functional principal components analysis with regularization
%
%  Arguments:
%  FD        ... Functional data object (a struct object)
%  PCABASIS  ... Basis for principal component weight functions
%  NHARM     ... Number of principal components to be kept. Default 2
%  LAMBDA    ... Smoothing or regularization parameter. Default 0
%  LFD       ... The order of the derivative, or a linear differential
%                operator structure, to be penalized.  Default 2.
%  CENTERFNS ... If 1, the mean function is first subtracted from each function
%                1 is the default.
%
%  Returns:  
%  A struct object PCASTR with the fields:
%  HARMFD  ... A functional data object for the harmonics or eigenfunctions
%  EIGVALS ... The complete set of eigenvalues
%  HARMSCR ... A matrix of scores on the principal components or harmonics
%  VARPROP ... A vector giving the proportion of variance explained
%                 by each eigenfunction
%  MEANFD  ... A functional data object giving the mean function
%

%  Last modified:  14 January 2003

  if ~isa_fd(fd)
    error ('First argument is not a functional data object.');
  end

%  set up default values

  if nargin < 6
    centerfns = 1;
  end
  if nargin < 5
    Lfd = int2Lfd(2);
  end
  if nargin < 4
    lambda = 0;
  end
  if nargin < 3
    nharm = 2;
  end
  if nargin < 2
      pcabasis = getbasis(fd);
  end
  
  %  check arguments

  if ~isa_basis(pcabasis)
    error ('Second argument is not a basis object.');
  end

  if ~isa_Lfd(Lfd)
    error (['LFD is not a linear differential operator object.');
  end

  %  compute mean function and center if required

  meanfd = mean(fd);
  if centerfns ~= 0
    fd = center(fd);
  end

  coef   = getcoef(fd);
  coefd  = size(coef);
  nbasis = coefd(1);
  nrep   = coefd(2);
  ndim   = length(coefd);

  if nrep < 2
    error('PCA not possible without replications');
  end

  type     = getbasistype(pcabasis);
  rangeval = getbasisrange(pcabasis);

  nrep   = coefd(2);
  if nrep < 2
    error ('PCA not possible without replications.');
  end

  if ndim == 3
    nvar  = coefd(3);
    ctemp = zeros(nvar*nbasis,nrep);
    for j = 1:nvar
      index = (1:nbasis) + (j-1)*nbasis;
      ctemp(index,:) = coef(:,:,j);
    end
  else
    nvar = 1;
    ctemp = coef;
  end

  %  set up cross product and penalty matrices

  fdbasis = getbasis(fd);
  Wmat = inprod(fdbasis,pcabasis);
  Jmat = eval_penalty(pcabasis, 0);
  Cmat = ctemp*ctemp'*Wmat./nrep;
  if lambda > 0
    Kmat = eval_penalty(pcabasis, Lfd);
    Wmat = Jmat + lambda .* Kmat;
  end

  %  compute the Choleski factor of Jmat

  Lmat    = chol(Jmat);
  Lmatinv = inv(Lmat);

  %  set up matrix for eigenanalysis

  if nvar == 1
    if lambda > 0
      Cmat = Lmatinv' * Jmat * Cmat * Jmat * Lmatinv;
    else
      Cmat = Lmat * Cmat * Lmat';
    end
  else
    for i = 1:nvar
      indexi =   (1:nbasis) + (i-1)*nbasis;
      for j = 1:nvar
        indexj = (1:nbasis) + (j-1)*nbasis;
        if lambda > 0
          Cmat(indexi,indexj) = ...
            Lmatinv' * Jmat * Cmat(indexi,indexj) * Jmat * Lmatinv;
        else
          Cmat(indexi,indexj) = Lmat * Cmat(indexi,indexj) * Lmat';
        end
      end
    end
  end

  % Eigenanalysis

  Cmat = (Cmat + Cmat')./2;
  [eigvecc, eigvalc]  = eig(Cmat);
  eigvalc = diag(eigvalc);
  [eigvals, indsrt] = sort(eigvalc);
  neig = length(eigvals);
  eigvals = eigvals(neig + 1 - (1:neig));
  eigvecs = eigvecc(:,indsrt);
  indx = nvar*nbasis + 1 - (1:nharm);
  eigvecs = eigvecs(:,indx);
  sumvecs = sum(eigvecs);
  eigvecs(:,sumvecs < 0) = -eigvecs(:,sumvecs < 0);

  varprop = eigvals(1:nharm)./sum(eigvals);

  if nvar == 1
    harmcoef = Lmatinv * eigvecs;
    harmscr  = ctemp' * Lmat' * eigvecs;
  else
    harmcoef = zeros(nbasis,nharm,nvar);
    harmscr  = zeros(nrep,nharm);
    for j = 1:nvar
      index = (1:nbasis) + (j-1)*nbasis;
      temp  = eigvecs(index,:);
      harmcoef(:,:,j) = Lmatinv * temp;
      harmscr = harmscr + ctemp(index,:)' * Lmat' * temp;
    end
  end

  harmnames = getnames(fd);
  harmnames{2} = 'Harmonics';
  harmnames{3} = ['Harmonics for',harmnames{3}];

  harmfd.coef     = harmcoef;
  harmfd.basisobj = pcabasis;
  harmfd.fdnames  = harmnames;

  harmfd = class(harmfd, 'fd');

  pcastr.harmfd  = harmfd;
  pcastr.eigvals = eigvals;
  pcastr.harmscr = harmscr;
  pcastr.varprop = varprop;
  pcastr.meanfd  = meanfd;
 
