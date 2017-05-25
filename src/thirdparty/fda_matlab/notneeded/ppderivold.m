function [CoeffD] = ppderiv(Coeff, Deriv)
%PPDERIV computes the DERIV-th derivatives of the polynomials 
% with coefficients COEFF such that the i-th polynomial is
% COEFF(i,1)*x^(k-1) + COEFF(i,2)*x^(k-2) + ... + COEFF(i,k)
% It returns a matrix COEFFD with the same number of rows as COEFF, 
% but with k-DERIV columns such that the DERIV-th derivative 
% of the i-th polynomial is expressed as
% COEFFD(i,1)*x^(k-1-DERIV) + COEFFD(i,k-DERIV-1)*x + COEFFD(i,k-DERIV),
% Note that if k-DERIV < 1, then COEFFD is the zero vector,
% and if DERIV < 1 we are not differentiating.

if nargin < 2
    Deriv = 0;
end

[m,k] = size(Coeff); % k is the order of the polynomials.

% If DERIV is not a positive integer, we are not differentiating.
if Deriv < 1
    CoeffD = Coeff; 
    return
end

% Compute the coefficient of the DERIV-th derivative of the function

if k-Deriv < 1
    CoeffD = zeros(m,1); % The derivative is zero everywhere
    return;
else
    % initialize COEFFD with the coefficients from COEFF we will need
    CoeffD = Coeff(:,1:k-Deriv); 
    for j=1:k-2
        bound1 = max(1,j-Deriv+1);
        bound2 = min(j,k-Deriv);
        CoeffD(:,bound1:bound2) = (k-j)*CoeffD(:,bound1:bound2);
    end
end