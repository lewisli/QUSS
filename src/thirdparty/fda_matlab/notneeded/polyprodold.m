function convmat = polyprod(c1, c2)
% POLYCONV computes products of polynomials defined by columns of 
%   coefficient matrices C1 and C2

%  Last modified 3 April 2002

[ord1, n1] = size(c1);
[ord2, n2] = size(c2);
deg1 = ord1 - 1;
deg2 = ord2 - 1;

%  if the degrees are not equal, pad out the smaller matrix with 0s

if deg1 ~= deg2
    if deg1 > deg2
        c2 = [c2;zeros(deg1-deg2,n2)];
    else
        c1 = [c1;zeros(deg2-deg1,n1)];
    end
end
       
D = max([deg1,deg2]);
N = 2*D+1;

convmat = zeros(n1,n2,N);
for i=0:D-1
    ind = (0:i) + 1;
    convmat(:,:,i+1) = c1(ind,    :)'*c2(i-ind+2,:);
    convmat(:,:,N-i) = c1(D-ind+2,:)'*c2(D-i+ind,:);
end
ind = (0:D)+1;
convmat(:,:,D+1) = c1(ind,:)'*c2(D-ind+2,:);

if deg1 ~= deg2
    convmat = convmat(:,:,1:(deg1+deg2+1));
end
