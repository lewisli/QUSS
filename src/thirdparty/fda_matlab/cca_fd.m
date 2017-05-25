function [ccawtfd1, ccawtfd2, ccavar1, ccavar2, corrs] = ...
              cca_fd(fdobj1, fdobj2, ncan, ccafdPar1, ccafdPar2)
%  CCA_FD   Functional canonical correlation analysis with regularization.
%
%  Arguments:
%  FDOBJ1    ... Functional data object for the first  set of functions.
%  FDOBJ2    ... Functional data object for the second set of functions.
%  NCAN      ... Number of pairs of canonical variates to be found. Default 2
%  CCAFDPAR1 ... A functional parameter object for the canonical weight
%                functions for the first  set of functions.
%  CCAFDPAR2 ... A functional parameter object for the canonical weight
%                functions for the second set of functions.
%  Returns:  
%  CCAWTFD1 ... A functional data object for the canonical weight
%                functions for the first  set of functions.
%  CCAWTFD2 ... A functional data object for the canonical weight
%                functions for the second set of functions.
%  CCAVAR1  ... Canonical variate scores for first  set of functions.
%  CCAVAR2  ... Canonical variate scores for second set of functions.
%  CORRS    ... The corresponding set of canonical correlations.

%  Last modified on:  19 March 2005

%  check to see that all five arguments are present

if nargin < 5
    error('There are less than five arguments.');
end

%  check that functions have the same number of replications

coef1  = getcoef(fdobj1);
coef2  = getcoef(fdobj2);
coefd1 = size(coef1);
coefd2 = size(coef2);
nrep1  = coefd1(2);
nrep2  = coefd2(2);

if (nrep1 ~= nrep2)
    error('The numbers of replications are not equal.');
end

%  check that there are more one replication

if nrep1 < 2
    error('There is only one replication.');
end

nrep = nrep1;

%  Center functions 

fdobj1 = center(fdobj1);
fdobj2 = center(fdobj2);

%  get basis information

fdbasis1  = getbasis(fdobj1);
fdbasis2  = getbasis(fdobj2);

%   Set up essential cross product matrices

Jmat1 = eval_penalty(fdbasis1, int2Lfd(0));
Jmat2 = eval_penalty(fdbasis2, int2Lfd(0));
Jx    = (Jmat1 * coef1)';
Jy    = (Jmat2 * coef2)';
PVxx  = Jx' * Jx./nrep; 
PVyy  = Jy' * Jy./nrep;

%  add penalty if either LAMBDA positive

lambda1 = getlambda(ccafdPar1);
lambda2 = getlambda(ccafdPar2);

if lambda1 > 0
    Lfdobj1 = getLfd(ccafdPar1);
    Kmat1   = eval_penalty(fdbasis1, Lfdobj1);
    PVxx    = PVxx + lambda1 * Kmat1;
end
if lambda2 > 0 
    Lfdobj2 = getLfd(ccafdPar2);
    Kmat2   = eval_penalty(fdbasis2, Lfdobj2);
    PVyy    = PVyy + lambda2 * Kmat2;
end

%  set up matrix to be analyzed

Vxy   = Jx' * Jy./nrep;

%  do eigenanalysis

geigstr = geigen(Vxy, PVxx, PVyy);

%  set up canonical correlations and coefficients for weight functions

canwtcoef1 = geigstr.Lmat(:,1:ncan);
canwtcoef2 = geigstr.Mmat(:,1:ncan);

corrs = diag(geigstr.values);
% corrs = corrs(1:ncan);

%   Normalize the weight functions

for j = 1:ncan
    temp = squeeze(canwtcoef1(:,j));
    temp = temp./sqrt(sum(temp.^2));
    canwtcoef1(:,j) = temp;
    temp = squeeze(canwtcoef2(:,j));
    temp = temp./sqrt(sum(temp.^2));
    canwtcoef2(:,j) = temp;
end

%  set up final results in struct object CCASTR

fdnames1  = getnames(fdobj1);
fdnames2  = getnames(fdobj2);
ccawtfd1  = fd(canwtcoef1, fdbasis1, fdnames1);
ccawtfd2  = fd(canwtcoef2, fdbasis2, fdnames2);

ccavar1 = Jx * canwtcoef1;
ccavar2 = Jy * canwtcoef2;

