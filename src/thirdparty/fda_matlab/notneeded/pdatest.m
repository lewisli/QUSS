%  test pda and pdan using solutions to some random first order linear equations
%            f + wx + Dx = 0
%  having the solution

%  x(t) = C exp(-int w)(t) + exp(-int w)[int [exp(int w) f]](t)

addpath('../fdaM')

%  define sample size

N = 10;

%  define number of forcing functions

nforce = 1;

%  set up forcing functions to vary randomly around 0

fnbasis = 5;
fbasis = create_bspline_basis([0,1],fnbasis);

fcoef0 = zeros(fnbasis,N,nforce);
fsigma = 1.0;
fcoef  = fcoef0 + fsigma.*randn(fnbasis,N,nforce);

ffd = fd(fcoef,fbasis);

plot(ffd);

%  set up x coefficient functions to vary randomly around 1

wnbasis = 5;
wbasis = create_bspline_basis([0,1],wnbasis);

wcoef0 = ones(wnbasis,N);
wsigma = 1.0;
wcoef  = wcoef0 + wsigma.*randn(wnbasis,N);

wfd = fd(wcoef,wbasis);

plot(wfd);

%  compute solutions to the N differential equations

tfine = (0:0.01:1)';
onesn = ones(length(tfine),1);

fmat   = eval_fd(ffd, tfine);
wmat   = eval_fd(wfd, tfine);

iwmat  = 0.01.*(cumsum(wmat) - ...
                0.5.*(onesn*wmat(1,:) + wmat));
eiwmat = exp(-iwmat);

xmat = eiwmat;  %  homogeneous solution

for iforce=1:nforce
    integrand = fmat(:,:,iforce)./eiwmat;
    intmat = 0.01.*(cumsum(integrand) - ...
                0.5.*(onesn*integrand(1,:) + integrand));
    xmat = xmat + eiwmat.*intmat;
end

subplot(1,1,1)
plot(tfine, xmat, 'b-', tfine, eiwmat, 'r-')

xnbasis = 23;
xbasis = create_bspline_basis([0,1],xnbasis);
xfd = data2fd(xmat, tfine, xbasis);

plot(xfd)

%  ---------------------------------------------------------------------------
%                        Test PDASCALAR
%  ---------------------------------------------------------------------------

norder  = 1;
flambda = 0;
wlambda = 0;
ffd0 = fd(zeros(fnbasis,1),fbasis);
wfd0 = fd( ones(wnbasis,1),wbasis);
festimate = 1;
westimate = 1;

%  PDA

[ffd, wfd, resfd] = pdascalar(xfd, norder, ...
                 fbasis, flambda, ffd0, festimate, ...
                 wbasis, wlambda, wfd0, westimate);

%  plot estimated forcing function and estimated weight function
subplot(2,1,1)
plot(ffd)
subplot(2,1,2)
plot(wfd)

%  compute derivatives
Dxfd = deriv(xfd,1);

%  plot derivatives and residuals
subplot(211)
plot(Dxfd)
subplot(212)
plot(resfd)

%  compute sums of squares of derivatives and of residuals
SSE0fd = sum(Dxfd.*Dxfd);
SSE1fd = sum(resfd.*resfd);

%  plot sums of squares
subplot(211)
plot(SSE0fd)
subplot(212)
plot(SSE1fd)

%  compute squared multiple correlation function

RSQfd = (SSE0fd - SSE1fd)./SSE0fd;
subplot(1,1,1)
plot(RSQfd)

%  ---------------------------------------------------------------------------
%                        Test PDACELL
%  ---------------------------------------------------------------------------
% 
%  define sample size

N = 10;

%  define number of u-variabless

nu = 2;

%  set up u-variables to vary randomly around 0

unbasis = 5;
ubasis  = create_bspline_basis([0,1],unbasis);

ucoef0 = zeros(unbasis,N,nu);
usigma = 1.0;
ucoef  = ucoef0 + usigma.*randn(unbasis,N,nu);

ufd = fd(ucoef,ubasis);

plot(ufd);

%  set up x coefficient functions w(t) to vary randomly around 1

wnbasis = 5;
wbasis  = create_bspline_basis([0,1],wnbasis);
wnbasis = 1;
wbasis  = create_constant_basis([0,1]);

wcoef0 = ones(wnbasis,1);
wsigma = 1.0;
wcoef  = wcoef0 + wsigma.*randn(wnbasis,1);

wfd = fd(wcoef,wbasis);

plot(wfd);

%  compute solutions to the N differential equations

delta = 0.001;
tfine = (0:delta:1)';
onesn = ones(length(tfine),1);

wmat   = eval_fd(wfd, tfine);
iwmat  = delta.*(cumsum(wmat) - ...
                0.5.*(onesn*wmat(1,:) + wmat));
eiwmat = exp(-iwmat);

xmat = eiwmat*ones(1,N);  %  homogeneous solution

umat = eval_fd(ufd, tfine);
if nu == 1
    integrand = umat./xmat;
    intmat = delta.*(cumsum(integrand) - ...
                0.5.*(onesn*integrand(1,:) + integrand));
    xmat = xmat - (eiwmat*ones(1,N)).*intmat;
else    
    for i=1:nu
        integrand = umat(:,:,i)./xmat;
        intmat = delta.*(cumsum(integrand) - ...
                    0.5.*(onesn*integrand(1,:) + integrand));
        xmat = xmat - (eiwmat*ones(1,N)).*intmat;
    end
end

subplot(1,1,1)
plot(tfine, xmat, 'b-', tfine, eiwmat, 'r-')

xnbasis = 23;
xbasis = create_bspline_basis([0,1],xnbasis);
xfd = data2fd(xmat, tfine, xbasis);

subplot(1,1,1)
plot(xfd)

%  ---------------------------------------------------------------------------
%  set up the PDACELL

clear awtcell
clear bwtcell

norder  = 1;

anbasis = 1;
abasis  = create_constant_basis([0,1]);

awtstruct.fd = fd(zeros(anbasis,1),abasis);
awtstruct.estimate = 1;
for i=1:nu
    awtcell{1,i} = awtstruct;
end;

bnbasis = 1;
bbasis  = create_constant_basis([0,1]);

bwtstruct.fd = fd(zeros(bnbasis,1),bbasis);
bwtstruct.estimate = 1;
bwtcell{1,1,1} = bwtstruct;

clear xfdcell

xfdcell{1} = xfd;

clear ufdcell

if nu == 1
    ufdcell{1} = ufd;
else
    for i=1:nu
        ufdcell{i} = ufd(:,i);
    end
end

%  PDA

[afdcell, bfdcell, resfdcell] = ...
    pdacell(xfdcell, ufdcell, awtcell, bwtcell, norder);

%  plot estimated u-variables weight

subplot(211)
plot(afdcell{1})
subplot(212)
plot(afdcell{2})
%plot(bfdcell{1,1})

subplot(111)
plot(bfdcell{1,1})

%  compute derivatives

Dxfd = deriv(xfd,1);

%  plot derivatives and residuals
subplot(211)
plot(Dxfd)
subplot(212)
plot(resfdcell{1})

%  compute sums of squares of derivatives and of residuals
SSE0fd = sum(Dxfd.*Dxfd);
SSE1fd = sum(resfdcell{1}.*resfdcell{1});

%  plot sums of squares
subplot(211)
plot(SSE0fd)
subplot(212)
plot(SSE1fd)

%  compute squared multiple correlation function

RSQfd = (SSE0fd - SSE1fd)./SSE0fd;
subplot(1,1,1)
plot(RSQfd)

