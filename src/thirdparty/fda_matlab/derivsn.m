function dy = derivsn(tnow, y) 
% DERIVS sets up the 1st order system corresponding to a
%   coupled linear differential operator defined by wfd for n variables.
  n = 3; 
  global wfd;   
  w  = squeeze(eval_fd(wfd, tnow));
  nv = length(y);
  m  = nv/n;
  wmat = zeros(nv, nv);
  wmat(1:(nv-m),(m+1):nv) = eye(nv-m);
  wmat((nv-m+1):nv,:) = -w';
  dy = wmat * y;

