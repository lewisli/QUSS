function bsplinemat = BsplineM(x, breaks, norder, nderiv, sparsewrd)
%   Matrix of values of Mth derivative of B-splines corresponding to args X.
%
%   X         ... a nondecreasing set of argument values
%   BREAKS    ... a set of strictly increasing break values containing all X
%                    within their range
%   NORDER    ... a order of splines (degree + 1)
%   NDERIV    ... order of derivative required (< norder)
%   SPARSEWRD ... if 1, return in sparse form

%  Last modified 30 May 2002

if ~isempty(find(diff(breaks) <= 0))
   error('The knot sequence breaks should be increasing.')
end
if ~isempty(find(diff(x)<0))
   error('The point sequence x should be nondecreasing.')
end

%  Compute the number  n  of B-splines of order K supported by the given
%  knot sequence and return empty matrix in case there aren't any.

% Set default values

if nargin < 5 sparsewrd = 0; end
if nargin < 4 nderiv    = 0; end
if nargin < 3 norder    = 4; end
if nderiv >= norder
   error('Order of derivative not less than order of spline');
end

% set some abbreviations

k   = norder;        % order of splines
km1 = k-1;
nb  = length(breaks); % number of break points
nx  = length(x);     % number of argument values
nd  = nderiv+1;      % nd is order of derivative plus one
ns  = nb - 2 + k;    % number of splines to compute
if ns < 1 
   fprintf('There are no B-splines for the given input.\n')
   bsplinemat = []; 
   return
end
onenx = ones(nx,1);
onenb = ones(k, 1);
onens = ones(ns,1);

%  augment break sequence to get knots by adding a K-fold knot at each end.

if size(breaks,1) > 1 breaks = breaks'; end
knots  = [breaks(1)*ones(1,km1), breaks, breaks(nb)*ones(1,km1)]';
nbasis = length(knots) - k;

%  For each  i , determine  left(i)  so that  K <= left(i) < nbasis+1 , and,
% within that restriction,
%        knots(left(i)) <= pts(i) < knots(left(i)+1) .

if size(x,1) == 1 x = x'; end
[ignored,index] = sort([knots(1:nbasis)', x']);
pointer = find(index>length(knots(1:nbasis)'))-[1:length(x)];
left = max([pointer; k*onenx']);  

% compute bspline values and derivatives, if needed:

% initialize the  b  array.
temp = [1, zeros(1,km1)]; 
b   = temp(ones(nd*nx,1),:);
nxs  = nd*[1:nx];
% run the recurrence simultaneously for all  x(i) .
% First, bring it up to the intended level:
for j=1:k-nd
   saved = zeros(nx,1);
   for r=1:j
      leftpr    = left + r;
      tr        = knots(leftpr) - x;
      tl        = x - knots(leftpr-j);
      term      = b(nxs,r)./(tr+tl);
      b(nxs,r) = saved + tr.*term;
      saved     = tl.*term;
   end
   b(nxs,j+1) = saved;
end

% save the B-spline values in successive blocks in  b .

for jj=1:nd-1
   j = k - nd + jj; 
   saved = zeros(nx,1); 
   nxn = nxs - 1;
   for r=1:j
      leftpr    = left + r;
      tr        = knots(leftpr) - x;
      tl        = x - knots(leftpr-j);
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
   nxs = reshape(temp,(nd-1-jj+1)*nx,1);
   for r=j:-1:1
      leftpr = left + r;
      temp   = ones(nd-jj,1)*(knots(leftpr) - knots(leftpr-j)).'/j;
      b(nxs,r)   = -b(nxs,r)./temp(:);
      b(nxs,r+1) =  b(nxs,r+1) - b(nxs,r);
   end
end

% Finally, zero out all rows of  b  corresponding to x outside the basic
% interval,  [breaks(1) .. breaks(nb)] .

index = find(x<breaks(1)|x>breaks(nb));
if ~isempty(index)
   temp = [1-nd:0].'*ones(1,length(index))+nd*ones(nd,1)*index(:).';
   b(temp(:),:) = zeros(nd*length(index),k);
end

% set up output matrix bsplinemat

width = max([ns,nbasis]) + km1 + km1;
cc = zeros(nx*width,1);
index = [1-nx:0].'*onenb' + ...
        nx*(left.'*onenb' + onenx*[-km1:0]);
cc(index) = b(nd*[1:nx],:);
% (This uses the fact that, for a column vector  v  and a matrix  A ,
%  v(A)(i,j)=v(A(i,j)), all i,j.)
bsplinemat = reshape(cc([1-nx:0].'*onens' + ...
                        nx*onenx*([1:ns])), nx, ns);
                    
%  store in sparse mode if required

if sparsewrd, bsplinemat = sparse(bsplinemat); end
