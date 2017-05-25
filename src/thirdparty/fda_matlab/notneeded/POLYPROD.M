function convmat = polyprod(Coeff1, Coeff2)
% POLYCONV computes products of polynomials defined by columns of 
%   coefficient matrices Coeff1 and Coeff2

%  Last modified 30 October 2002

[polyorder1, norder] = size(Coeff1);
[polyorder2, norder] = size(Coeff2);
ndegree1 = polyorder1 - 1;
ndegree2 = polyorder2 - 1;

%  if the degrees are not equal, pad out the smaller matrix with 0s

if ndegree1 ~= ndegree2
    if ndegree1 > ndegree2
        Coeff2 = [Coeff2;zeros(ndegree1-ndegree2,norder)];
    else
        Coeff1 = [Coeff1;zeros(ndegree2-ndegree1,norder)];
    end
end

%  find order of the product

D = max([ndegree1,ndegree2]);  % maximum degree
N = 2*D+1;                     % order of product

%  compute the coefficients for the products

convmat = zeros(norder,norder,N);
for i=0:D-1
    ind = (0:i) + 1;
    convmat(:,:,i+1) = Coeff1(ind,    :)'*Coeff2(i-ind+2,:);
    convmat(:,:,N-i) = Coeff1(D-ind+2,:)'*Coeff2(D-i+ind,:);
end
ind = (0:D)+1;
convmat(:,:,D+1) = Coeff1(ind,:)'*Coeff2(D-ind+2,:);

if ndegree1 ~= ndegree2
    convmat = convmat(:,:,1:(ndegree1+ndegree2+1));
end
