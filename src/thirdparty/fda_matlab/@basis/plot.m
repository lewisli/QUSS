function h=plot(basisobj, nx)
%  Plot a basis object.

%  last modified 10 May 2004

if nargin < 2, nx = 101;  end
 
rangex   = getbasisrange(basisobj);
x        = linspace(rangex(1),rangex(2),nx)';
basismat = full(eval_basis(x, basisobj));
h=plot (x, basismat, '-');
minval = min(min(basismat));
maxval = max(max(basismat));
if minval == maxval
    if abs(minval) < 1e-1
        minval = minval - 0.05;
        maxval = maxval + 0.05;
    else
        minval = minval - 0.05*minval;
        maxval = maxval + 0.05*minval;
    end
end
axis([x(1), x(nx), minval, maxval])
