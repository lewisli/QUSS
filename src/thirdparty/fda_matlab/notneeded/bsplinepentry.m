nbasis = 3;
norder = 2;

basis = create_bspline_basis([0,1],nbasis,norder);

nderiv = 0;

penmat = bsplinepen(basis,nderiv);

nbasis = 4;
norder = 3;

basis = create_bspline_basis([0,1],nbasis,norder);

nderiv = 0;

penmat2 = bsplinepen(basis,nderiv)

penmat = inprod(basis, basis, nderiv, nderiv)

coef = randn(nbasis,1);
fd1 = fd(coef,basis);


x = 0:0.1:1;

fx = eval_fd(x, fd1, nderiv);

fx
