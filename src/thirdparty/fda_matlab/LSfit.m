function beta = LSfit(yroot, zmat, f, wtrtmt)
xmat  = [zmat,f];
ncol  = size(xmat,2);
xroot = xmat.*wtrtmt;
[Q, R, E] = qr(xroot);
tol       = size(xmat,1)*3e-16*abs(R(1,1));
y         = Q'*yroot;
for j=1:ncol
  if abs(R(j,j)) < tol
    if R(j,j) < 0, R(j,j) = -tol;  else R(j,j) = tol; end
  end
end
beta      = R\y;
beta      = E*beta;

