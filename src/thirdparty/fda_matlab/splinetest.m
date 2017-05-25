
x = linspace(0,1,101)';


y0 = sin(4*pi*x);
%y0 = ones(length(x),1);

n = 5;

onen = ones(1,n);

sigma = 1;

y = y0*onen + sigma.*randn(101,n);

basis = create_bspline_basis([0,1],5,4);
yfd = data2fd(y,x,basis);

plot(yfd)

fd1 = yfd;
fd2 = yfd;
deriv1 = 2;
deriv2 = 1;

tic;
prodmat1 = inprod_bspline(yfd,yfd,2,1)
toc;
tic;
prodmat2 = inprod(yfd,yfd,2,1)
toc;

max(max(abs(prodmat1-prodmat2)))/max(max(abs(prodmat1)))

yfine = eval(yfd,x);
ysqr = yfine.^2;
inty  = 0.01.*(sum(ysqr)        -0.5.*(ysqr( 1)+ysqr(101)))
inty1 = 0.01.*(sum(ysqr( 1: 51))-0.5.*(ysqr( 1)+ysqr( 51)))
inty2 = 0.01.*(sum(ysqr(51:101))-0.5.*(ysqr(51)+ysqr(101)))


penaltymat1 = bsplinepen(basis,2);
penaltymat2 = sparse(inprod(basis,basis,Lfd,Lfd));

penaltymat1-penaltymat2

