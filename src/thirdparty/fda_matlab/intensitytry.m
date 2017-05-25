N = 100;
t = sort([  N*sort(rand(N,1)); ...
          2*N*sort(rand(N,1))+N; ...
            N*sort(rand(N,1))+3*N;]);

nbasis = 23;
norder = 4;
wbasis = create_bspline_basis([0,max(t)], nbasis, norder);

Wfd0    = fd(zeros(nbasis,1), wbasis);
lambda  = 1e5;
Lfdobj  = 2;
Wfd0Par = fdPar(Wfd0, Lfdobj, lambda);

Wfdobj = intensity_fd(t, Wfd0Par);

subplot(2,1,1)
plot(Wfdobj)
subplot(2,1,2)
plot(t, exp(eval_fd(t, Wfdobj)))


