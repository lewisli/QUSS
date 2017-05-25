function msplinemat = MsplineM(x, knots, norder, nderiv, sparsewrd)
%   Values of NDERIV-th derivative of M-splines corresponding to args X.
%
%   X         ... a nondecreasing set of argument values
%   KNOTS     ... a set of strictly increasing knot values containing all X
%                    within their range
%   NORDER    ... a order of splines (degree + 1)
%   NDERIV    ... order of derivative required (< norder)
%                 if -1, value of I-spline is returned
%                 if  0, value of M-spline is returned
%                 if  m, value of mth derivative of M-spline is returned
%   SPARSEWRD ... if 1, return in sparse form

if ~isempty(find(diff(knots) <= 0))
   error('The knot sequence KNOTS should be increasing.')
end
if ~isempty(find(diff(x)<0))
   error('The point sequence x should be nondecreasing.')
end

%  Compute the number  n  of B-splines of order K supported by the given
%  knot sequence and return empty matrix in case there aren't any.

% Settle the options:
if nargin < 5 sparsewrd = 0; end
if nargin < 4 nderiv    = 0; end
if nargin < 3 norder    = 4; end
if nderiv >= norder
   error('Order of derivative not less than order of spline');
end

%set some abbreviations

% order of splines
if nderiv == -1
  k = norder + 1;
else
  k   = norder;        
end
km1 = k-1;
nk  = length(knots); % number of break points
nx  = length(x);     % number of argument values
% nd is order of derivative plus one, or 1 if NDERIV = -1
if nderiv == -1
  nd  = 1; 
else
  nd = nderiv + 1;
end     
ns  = nk - 2 + k;    % number of splines to compute
if ns < 1 
   fprintf('There are no B-splines for the given input.\n')
   bsplinemat = []; 
   return
end
onenx = ones(nx,1);
onenk = ones(k, 1);
onens = ones(ns,1);

%  augment knot sequence to provide a K-fold knot at each end.

if size(knots,1) > 1 knots = knots'; end
augknot = [knots(1)*ones(1,km1), knots, knots(nk)*ones(1,km1)]';
naug = length(augknot) - k;

%  For each  i , determine  left(i)  so that  K <= left(i) < naug+1 , and,
% within that restriction,
%        augknot(left(i)) <= pts(i) < augknot(left(i)+1) .

if size(x,1) == 1 x = x'; end
[ignored,index] = sort([augknot(1:naug)', x']);
pointer = find(index>length(augknot(1:naug)'))-[1:length(x)];
left = max([pointer; k*onenx']);  

% compute bspline values and derivatives, if needed:

% initialize the  b  array.
temp = [1 zeros(1,km1)]; 
b   = temp(ones(nd*nx,1),:);
nxs  = nd*[1:nx];
% run the recurrence simultaneously for all  x(i) .
% First, bring it up to the intended level.  
%  This is order NORDER + 1 if NDERIV = -1.
if nderiv > -1
  kmnd = k - nd;
else
  kmnd = k - nd + 1;
end
for j=1:k-nd
   saved = zeros(nx,1);
   for r=1:j
      leftpr    = left + r;
      tr        = augknot(leftpr) - x;
      tl        = x - augknot(leftpr-j);
      term      = b(nxs,r)./(tr+tl);
      b(nxs,r) = saved + tr.*term;
      saved     = tl.*term;
   end
   if nderiv = 0
     %  convert B-spline to M-spline by re-normalizing
     b(nxs,j+1) = saved.*k./(augknot(leftpr)-augknot(leftpr-j-1);
   else
     b(nxs,j+1) = saved;
   end
end

% Save the B-spline values in successive blocks in  b if NDERIV > 0.

for jj=1:nd-1
   j = k - nd + jj; 
   saved = zeros(nx,1); 
   nxn = nxs - 1;
   for r=1:j
      leftpr    = left + r;
      tr        = augknot(leftpr) - x;
      tl        = x - augknot(leftpr-j);
      term      = b(nxs,r)./(tr+tl);
      b(nxn,r) = saved + tr.*term;
      saved     = tl.*term;
   end
   b(nxn,j+1) = saved;
   nxs = nxn;
end

% now use the fact that derivative values can be obtained by differencing:

for jj=nd-1:-1:1
   j = k - jj;
   temp = [jj:nd-1].'*onenx' + ones(nd-jj,1)*nxn; 
   nxs=temp(:);
   for r=j:-1:1
      leftpr = left + r;
      temp   = ones(nd-jj,1)*(augknot(leftpr) - augknot(leftpr-j)).'/j;
      b(nxs,r)   = -b(nxs,r)./temp(:);
      b(nxs,r+1) =  b(nxs,r+1) - b(nxs,r);
   end
end

%  If NDERIV > 0, renormalize to make into M-spline derivative

if nderiv > 0
  j = k - 1;
  leftpr = left + j;
  temp = [1:nd-1].'*onenx' + ones(nd-1,1)*nxn; 
  nxs=temp(:);
  b(nxs,1) = b(nxs).*norder./(augknot(leftpr)-augknot(leftpr-j-1);
end

% Finally, zero out all rows of  b  corresponding to x outside the basic
% interval,  [knots(1) .. knots(nk)] .

index = find(x<knots(1)|x>knots(nk));
if ~isempty(index)
   temp = [1-nd:0].'*ones(1,length(index))+nd*ones(nd,1)*index(:).';
   b(temp(:),:) = zeros(nd*length(index),k);
end

% set up output matrix bsplinemat

width = max([ns,naug]) + km1 + km1;
cc = zeros(nx*width,1);
index = [1-nx:0].'*onenk' + ...
        nx*(left.'*onenk' + onenx*[-km1:0]);
cc(index) = b(nd*[1:nx],:);
% (This uses the fact that, for a column vector  v  and a matrix  A ,
%  v(A)(i,j)=v(A(i,j)), all i,j.)
msplinemat = reshape(cc([1-nx:0].'*onens' + ...
                        nx*onenx*([1:ns])), nx, ns);
                        
%  If NDERIV is -1, compute I-spline values

if nderiv == -1
  sum = zeros(nxs,1);
  for i=1:norder
    j = norder - i + 1;
    sum = sum + msplinemat(:,j+1);
    msplinemat(:,j) = sum;
  end
end
  

if sparsewrd msplinemat = sparse(msplinemat); end
