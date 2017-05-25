%  test LOGLAM2GCV function

%  set up a test problem

t = (0:0.01:1)';
y = sin(4*pi*t) + randn(101,1).*0.2;
t = t.*100;

nbasis   = 103;
basisobj = create_bspline_basis([0,100], nbasis);

Lfdobj   = 2;
lambda   = 1e-4;
fdParobj = fdPar(basisobj, Lfdobj, lambda);

Z = eval_basis(t, basisobj);
R = sparse(eval_penalty(basisobj, Lfdobj)); 

loglam = 0:0.5:3;

for ilam = 1:length(loglam)
    [GCV, DGCV, df] = loglam2gcv(loglam(ilam), y, Z, R);
    display([loglam(ilam), GCV, df])
end

fdParobj = putlambda(fdParobj, 10^1.5);
fdobj = smooth_basis(t, y, fdParobj);
plotfit_fd(y, t, fdobj)

%  This doesn't work because Z'Z is not of full rank.
%  Two eigenvalues are Inf

H = Z'*Z;
[V, D] = eig(full(R), full(Z'*Z));
[D, index] = sort(diag(D));
log10(D)

%  now make Z'Z of full rank by dropping two knots

nbasis   = 101;
breaks   = [0, 3:98, 100];
basisobj = create_bspline_basis([0,100], nbasis, 4, breaks);

Z = eval_basis(t, basisobj);
H = Z'*Z;

R = eval_penalty(basisobj, Lfdobj); 
R = full(R);
R = (R' + R)./2;
H = full(Z'*Z);

[V, D] = eig(R, H);
D = diag(D);
max(max(abs(V'*R*V - diag(D))))
max(max(abs(V'*H*V - eye(101))))

[Dsrt,srtind] = sort(D);
[log10([D, Dsrt]),srtind]

index = 1:2;
Dsrt(index) = 0;

Vsrt = V(:,srtind);
max(max(abs(Vsrt'*R*Vsrt - diag(Dsrt))))
max(max(abs(Vsrt'*H*Vsrt - eye(101))))

lambda = 1;

M = H + lambda*R;

Minv = inv(M);

Meigvals = ones(nbasis,1) + lambda.*Dsrt;

Minveigvals = 1./Meigvals;

Minvhat = Vsrt*diag(Minveigvals)*Vsrt';

max(max(abs(Minv - Minvhat)))

A = Z*Minv*Z';

ZV = Z*V;

Ahat = ZV*diag(Minveigvals)*ZV';

max(max(abs(A - Ahat)))

loglam = 0:0.5:3;

for ilam = 1:length(loglam)
    [GCV, DGCV, df] = loglam2gcv(loglam(ilam), y, Z, R);
    display([loglam(ilam), GCV, df])
    [GCV, DGCV, df] = loglam2gcv_pdq(loglam(ilam), y, Z, V, D);
    display([loglam(ilam), GCV, df])
end



