%  test llbda2gcv function

Y = y;
Z = eval_basis(tval, basis);
R = sparse(eval_penalty(basis, 4)); 

[GCV1, DGCV1, df1, trA1, DtrA1, trMinv1, DtrMinv1] = llbda2gcv(-9.000, Y, Z, R, factor);
[GCV2, DGCV2, df2, trA2, DtrA2, trMinv2, DtrMinv2] = llbda2gcv(-3.999, Y, Z, R, factor);
[DGCV1,(GCV2-GCV1)/0.001]
full([DtrA1,(trA2-trA1)/0.001])
full([DtrMinv1,(trMinv2-trMinv1)/0.001])

delta = 0.01;
log10lam = -6:delta:-5;
m = length(log10lam);
results = zeros(m,6);
for i=1:m
    [GCV, DGCV, df, trA, DtrA, trMinv, DtrMinv] = ...
        llbda2gcv(log10lam(i), Y, Z, R, factor);
    results(i,1) = GCV;
    results(i,2) = DGCV;
    results(i,3) = trA;
    results(i,4) = DtrA;
    results(i,5) = trMinv;
    results(i,6) = DtrMinv;
end
for j=1:6
    subplot(3,2,j)
    plot(log10lam, results(:,j), 'r-')
    if floor(j/2)*2 == j
        line(log10lam(2:m), diff(results(:,j-1)')./delta)
    end
end

