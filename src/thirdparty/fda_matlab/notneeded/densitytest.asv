N = 51;
x = randn(N,1);
x = [  0.2402, 0.8777,  0.2234, -1.8073, -2.6152, -0.1879, -0.1391, -0.0609,  0.3443, ...
      -2.0978,-1.1409,  1.9057,  0.2563,  1.5610,  0.3364,  1.3327, -0.4593, -0.1835, ...
       1.8356, 0.6365, -1.3793, -0.2455,  0.3575, -0.1923, -1.0962,  0.3221,  1.1339, ...
       0.6079,-0.5165, -3.2414,  1.6928,  1.3289, -0.5835,  0.1294,  0.9806,  1.5827, ...
       0.7757,-0.1011, -0.8097,  0.5285, -0.5784, -0.7105,  0.1749,  0.8444,  0.1769, ...
       0.0388,-0.4058, -0.8430, -1.5724,  0.8137, -1.3632];
 
nbasis = 13;
norder = 6;
range  = [min([-3,min(x)]),max([3,max(x)])];
basis  = create_bspline_basis(range,nbasis,norder);

Wfd0 = fd(zeros(nbasis,1), basis);

dbglev  = 1; 
active  = 2:nbasis;
conv    = 1e-2;
iterlim = 20; 
lambda  = 1e-2;
Lfd     = 3;

[Wfd, C, Fstr, iternum, iterhist] = densityfd(x, Wfd0, Lfd, lambda, ...
                                         iterlim, conv, active, dbglev);

p = exp(eval_fd(Wfd, xfine))/C;

subplot(1,3,1)
plot(Wfd)
xlabel('x')
ylabel('')
title('Function W(x)')
axis('square')

subplot(1,3,2)
plot(deriv(Wfd))
title('Function w(x) = DW(x)')
xlabel('x')
ylabel('')
axis('square')

xfine = linspace(-3,3,101);
subplot(1,3,3)
plot(xfine,p,'-',x,zeros(N,1),'b.')
title('Density Function p(x)')
xlabel('x')
ylabel('')
axis([range(1),range(2),0,max(p)])
axis('square')


x = linspace(-3,3,51)';
f = exp(-x.^2/2)/sqrt(2*pi);
plot(x,f)
x = [x,f];

[Wfd, C, Fstr, iternum, iterhist] = densityfd(x, Wfd0, Lfd, lambda);

