function fx = Trap(fx, delta)
%  TRAP ... values of integrated function by trapezoidal rule
  %  FX is a set of integrand values computed on 
  %    an equally spaced grid.
  %  DELTA is the width of the intervals.
  n = length(fx);
  if nargin<2, delta=1/(n-1); end
  fx = delta.*(cumsum(fx) - 0.5.*(fx(1)+fx));
