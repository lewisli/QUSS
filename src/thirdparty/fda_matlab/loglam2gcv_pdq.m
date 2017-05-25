function [GCV, DGCV, df] = ...
           loglam2gcv_pdq(loglambda, Y, Z, V, D, factor)
% LOGLAM2GCV computes the GCV criterion and its derivative
%  as a function of log_{10} \lambda where \lambda is the
%  smoothing parameter in a call to SMOOTH_BASIS
%  Arguments:
%  LOGLAMBDA ... logarithm to base 10 of \lambda
%  Y         ... data that are smoothed, n by N matrix where
%                n is the number of values per record that are
%                smoothed and N is the number of records.
%  Z         ... n by K matrix of basis function values where K
%                is the number of basis functions
%  V and D   ... these result from the generalized 
%                eigenanalysis  [V,D] = EIG(R, Z'Z)
%                where R is the K by K penalty matrix.  
%                D is a vector of eigenvalues (not a diagonal
%                matrix as produced by function EIG).  
%                Note that the rank of R is n - m where m is
%                largest derivative in the differential operator
%                defining the roughness penalty. It is
%                assumed that these m zero eigenvalues and the 
%                corresponding columns of V have been removed.  
%  FACTOR    ... Chong Gu's correction factor, suggested as 
%                1.2 or 1.4.  

%                This analysis assumes that Z'Z is of full rank.
%                This will not be the case, for example, if
%                the number of basis functions K exceeds the
%                number of sampling points n.  In such cases,
%                the slower but more general function
%                LOGLAM2GCV should be used.

%  Last modified 3 May 2004

%  In the following computation, no attempt is made to use
%  matrix decompositions since normally both Z and R will be
%  band-structured.

if nargin < 6, factor = 1;  end
if nargin < 5
    error('Less than five arguments.');
end

lambda = 10^loglambda;   %  lambda

Dlambda = log(10)*lambda; %  Derivative of lambda wrt log10 lambda

[n, N] = size(Y);

[K, nminusm] = size(V);

ZV = Z*V;

matfac = ones(n,1)*(1./(1 + lambda.*D'));

A = (ZV.*matfac)*ZV';    %  the smoothing or hat matrix

Yhat = A*Y;              %  matrix of smoothed values

Res  = Y - Yhat;         %  matrix of residuals

SSE = sum(sum(Res.^2));  %  error sum of squares

df = n - factor*trace(A);  %  degrees of freedom in the smooth

numGCV = SSE/n;          %  numerator of GCV

denGCV = (df/n)^2;       %  denominator of GCV

GCV = numGCV/denGCV;      %  GCV

DA = -Dlambda.*(ZV.*matfac.^2)*ZV'; %  derivative of A

DnumGCV = -2*trace(DA*Y*Res')/n;  %  derivative of SSE/n wrt lambda

DdenGCV = -2*df*factor*trace(DA)/n^2;

DGCV = ((denGCV*DnumGCV - numGCV*DdenGCV)/denGCV^2);

