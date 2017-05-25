%  test function posfd

N = 101;

x = linspace(0,1,N)';

y0 = exp(sin(4*pi*x));

sigma = 0.2;

y = y0 + randn(N,1).*sigma;

plot(x, y, 'o', x, y0, '--')

nbasis = 13;
norder = 4;
basis  = create_bspline_basis([0,1],nbasis,norder);

basismat = getbasismatrix(x, basis);
cvec = basismat\log(y);

Wfd0 = fd(cvec, basis);

Lfd     = 2;
lambda  = 1e-5;
iterlim = 20;
conv    = 1e-4;
dbglev  = 1;

[Wfd, Fstr, iternum, iterhist] = ...
   posfd(x, y, Wfd0, Lfd, lambda, iterlim, conv, dbglev);

plot(Wfd)

Wvec = eval_fd(Wfd, x);
fdvec = eval_pos(x, Wfd);
plot(x, y, 'o', x, fdvec, 'b-', x, y0, 'g--')

Dyhat = eval_pos(x, Wfd, 1);

plot(x, Dyhat, '-')

