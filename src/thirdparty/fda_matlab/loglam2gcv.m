function [GCV, DGCV, df] = loglam2gcv(loglambda, Y, Z, R, factor)
% LOGLAM2GCV computes the GCV criterion and its derivative
%  as a function of log_{10} \lambda where \lambda is the
%  smoothing parameter in a call to SMOOTH_BASIS
%  Arguments:
%  LOGLAMBDA ... logarithm to base 10 of \lambda
%  Y         ... data that are smoothed, n by N matrix where
%                n is the number of values per record that are
%                smoothed and N is the number of records.
%  Z         ... n by K matrix of basis functions where K
%                is the number of basis functions
%  R         ... K by K penalty matrix
%  FACTOR    ... Chong Gu's correction factor, suggested as 
%                1.2 or 1.4.  

%  Although it is faster to find the minimizing GCV value of
%  lambda with this function than by simplying smoothing the 
%  data over and over again, this function does carry a matrix
%  inversion of order K each time it is called.  

%  A still faster version is available under the name 
%  LOGLAM2GCV_PDQ if (1) matrix Z'Z is of full rank and
%  (2) a preliminary generalized eigenequation
%        RV = Z'ZVD, V'Z'ZV = I
%  is solved using [V,D] = EIG(R,Z'Z).  That version
%  avoids matrix inversion and any other O(n^3) calculation.

%  Nevertheless,  when R and Z'Z are band-structured, as they are
%  with a B-spline basis, or diagonal, as they are with Fourier
%  bases and equally spaced arguments, the matrix inversion
%  using sparse computation methods is fast, and this version
%  will consequently serve most purposes.

%  Last modified 29 April 2004

%  In the following computation, no attempt is made to use
%  matrix decompositions since normally both Z and R will be
%  band-structured.

if nargin < 5, factor = 1;  end
if nargin < 4
    error('Less than four arguments.');
end

logtrR = log10(trace(R));

if logtrR + loglambda > 10
    warning('Condition number is too large for stable results.');
end

lambda = 10^loglambda;   %  lambda

Dlambda = log(10)*lambda; %  Derivative of lambda  wrt log10 lambda

[n, N] = size(Y);

M = Z'*Z + lambda.*R;    %  perturbed cross-product matrix

Minv = inv(M);           %  inverse of this

ZMinv = Z*Minv;

A = ZMinv*Z';            %  smoothing or hat matrix

Yhat = A*Y;              %  matrix of smoothed values

Res  = Y - Yhat;         %  matrix of residuals

SSE = sum(sum(Res.^2));  %  error sum of squares

df = n - factor*trace(A);  %  degrees of freedom in the smooth

numGCV = SSE/n;          %  numerator of GCV

denGCV = (df/n)^2;       %  denominator of GCV

GCV = numGCV/denGCV;      %  GCV

DA = -Dlambda.*ZMinv*R*ZMinv'; %  derivative of A

DnumGCV = -2*trace(DA*Y*Res')/n;  %  derivative of SSE/n wrt lambda

DdenGCV = -2*df*factor*trace(DA)/n^2;

DGCV = ((denGCV*DnumGCV - numGCV*DdenGCV)/denGCV^2);

