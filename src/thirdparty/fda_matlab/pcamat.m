function pcastr = pcamat(fdmat, wt, nharm, centerfns)
%PCAMAT Functional principal components analysis of function values
%  This function resembles PCA for FD objects, but accepts as input
%  a matrix FDMAT of values of curves.  In this case no regularization
%  is possible, so it is assumed that the curve values themselves are
%  sufficiently smooth to ensure smooth eigenvectors.  
%  Only one variable is permitted.
%  PCAMAT also accepts a weight vector WT.
%
%  Arguments:
%  FDMAT     ... Values of a functional data object.  
%                Rows    correspond to a fine mesh of argument values,
%                Columns correspond to replications.
%  WT        ... A vector of weights.
%  NHARM     ... Number of principal components to be kept. Default 2
%  CENTERFNS ... If 1, the mean function is first subtracted from each function
%                1 is the default.
%
%  Returns:  
%  A struct object PCASTR with the fields:
%  HARMMAT ... A matrix of harmonic or eigenfunction values
%  EIGVALS ... The complete set of eigenvalues
%  SCORES  ... A matrix of scores on the principal components or harmonics
%  VARPROP ... A vector giving the proportion of variance explained
%                 by each eigenfunction

%  Last modified:  7 December 2000

[nfine, nrep] = size(fdmat);
if nrep < 2
    error('PCA not possible without replications');
end

if size(wt,1) ~= nfine
    error('WT does not have same number of rows as FDMAT.');
end
if size(wt,2) ~= 1
    error('WT does not have a single column.');
end
if any(wt <= 0)
    error('One or more values of WT are not positive.');
end

%  set up default values

if nargin < 4
    centerfns = 1;
end
if nargin < 3
    nharm = 2;
end
if nargin < 2
    wt = ones(nfine,1);
end

%  compute mean function and center if required

if centerfns ~= 0
    fdmean = mean(fdmat')';
    fdmat  = fdmat - fdmean*ones(1,nrep);
end

%  set up cross product and penalty matrices

wtmat = sqrt(wt)*ones(1,nrep);
ctemp = fdmat.*wtmat;
Cmat  = ctemp'*ctemp;

% Eigenanalysis

[eigvecc, eigvalc]  = eig(Cmat);
eigvalc = diag(eigvalc);
[eigvals, indsrt] = sort(eigvalc);
eigvals = eigvals(nrep + 1 - (1:nrep));
eigvecs = eigvecc(:,indsrt);
indx = nrep + 1 - (1:nharm);
eigvecs = eigvecs(:,indx);
sumvecs = sum(eigvecs);
eigvecs(:,sumvecs < 0) = -eigvecs(:,sumvecs < 0);

pcastr.eigvals = eigvals;
pcastr.varprop = eigvals(1:nharm)./sum(eigvals);
pcastr.harmmat = (fdmat * eigvecs)./sqrt(wt*ones(1,nharm));
pcastr.scores  = eigvecs * diag(eigvals(1:nharm));
