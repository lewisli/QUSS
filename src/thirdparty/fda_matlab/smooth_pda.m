%   use profPDA to analyze some data

%  Note:  to run this code, you must be in Matlab/fdaM

%  ----------------------------------------------------------
%                    try some simulations
%  ----------------------------------------------------------

%  set options for fminunc

options = optimset('LargeScale', 'off', 'Display', ...
                   'iter', 'MaxIter', 40);

options = optimset('LargeScale', 'off', 'Display', ...
                   'off', 'MaxIter', 20);
               
N     = 1
nsam  = 100;
% sigma = 0.2;
sigma = 0.088 %  error level for Corpus Christi data

%  set up range of log10 lambda's

% loglam = (-3:1:2)';

loglam = 3;

nlam   = length(loglam);

df     = zeros(nlam,nsam);
gcv    = zeros(nlam,nsam);
SSE    = zeros(nlam,nsam);
PENSSE = zeros(nlam,nsam);
bmat   = zeros(nlam,nsam,length(bvec0));
hessdg = zeros(nlam,nsam,length(bvec0));

for isam=1:nsam
    y = y0*ones(1,N) + randn(n,N)*sigma;
    Dmat = (basismat .* (wtvec * ones(1,nbasis)))' * y;
    bvec = bvec0;
%     tic;
    for ilam=1:nlam
%         disp(['log10 lambda = ',num2str(loglam(ilam))])
        lambda = 10^loglam(ilam);
        [bvec, fval, exitflag, output, grad, hessian] = ...
            fminunc(@profPDA, bvec, options, ...
                y, basisobj, basismat, Bmat, Dmat, ...
                bwtcell, awtcell, ufdcell, lambda);
        bmat(ilam,isam,:) = bvec';
        [SSE(ilam,isam), PENSSE(ilam,isam), penmat, fdobj, ...
                df(ilam,isam), gcv(ilam,isam)] = ...
            profPDA(bvec, y, basisobj, basismat, Bmat, Dmat, ...
            bwtcell, awtcell, ufdcell, lambda);
        hessdg(ilam,isam,:) = diag(inv(hessian));
        [isam, loglam(ilam), df(ilam,isam), gcv(ilam,isam)]
    end
%     toc;
end

%  display RMS and penalized RMS

sqrt(SSE./((N*(n-df))))

%  display DF and GCV values

[df, gcv]

%  display sqrt of diagonal of inverse of hessian

stderr = squeeze(sqrt(mean(hessdg)./N).*sigma);

%  display coefficient vectors

[bmat]

bindex = 1;

bmat(bindex,:,:)

% output some results

meanb = squeeze(mean(bmat(bindex,:,:)));
stdvb = squeeze(sqrt(var(squeeze(bmat(bindex,:,:)))))';

[bvec0, meanb, meanb-bvec0, stdvb, stderr]

% for ilam=1:nlam
for ilam=3    
    bvec = bmat(ilam,:)';
    
    %  set up linear differential operator
    
    [Lfdobj, bfdcell, afdcell] = ...
        bvec2Lfd(bvec, bwtcell, awtcell, ufdcell);
    
    %  plot the operator
    
    plot(Lfdobj)
    hold on
    plot(tval, -2*4*pi*cos(4*pi*tval), 'r--')
    hold off
    pause;
    
    %  smooth the data with the operator
    
    lambda = 10^loglam(ilam);
    fdParobj = fdPar(basisobj, Lfdobj, lambda);
    [fdobj, df, gcv, coef, SSE, penmat] = ...
        smooth_basis(y, tval, fdParobj);
    
    %  plot the data and fit
    
    subplot(1,1,1)
    plot(tval, y, '.', tval, eval_fd(tval, fdobj), 'b-', ...
         tval, y0, 'r--')
    xlabel('\fontsize{16} t')
    ylabel('\fontsize{16} x(t)')
    legend('\fontsize{16} Data', 'Estimate', 'True')
    pause;
    
    plot(tval, eval_fd(tval, fdobj, 1),    'b-', ...
         tval, 1.5.*exp(2*sin(4*pi*tval)), 'r--')
    xlabel('\fontsize{16} t')
    ylabel('\fontsize{16} Dx(t)')
    legend('\fontsize{16} Estimate', 'True')
    pause;
end

%  ----------------------------------------------------------
%                    analysis of actual data
%  ----------------------------------------------------------

%  set options for fminunc

options = optimset('LargeScale', 'off', 'Display', ...
                   'iter', 'MaxIter', 40);

options = optimset('LargeScale', 'off', 'Display', ...
                   'off', 'MaxIter', 20);
               
%  set up range of log10 lambda's

loglam = (0:1:3)';
nlam   = length(loglam);

df     = zeros(nlam,1);
gcv    = zeros(nlam,1);
SSE    = zeros(nlam,1);
PENSSE = zeros(nlam,1);
bmat   = zeros(nlam,length(bvec0));
hessdg = zeros(nlam,length(bvec0));

for ilam=1:nlam
    lambda = 10^loglam(ilam);
    [bvec, fval, exitflag, output, grad, hessian] = ...
        fminunc(@profPDA, bvec0, options, ...
        y, basisobj, basismat, Bmat, Dmat, ...
        bwtcell, awtcell, ufdcell, lambda);
    bmat(ilam,:) = bvec';
    [SSE(ilam), PENSSE(ilam), penmat, fdobj, ...
            df(ilam), gcv(ilam)] = ...
        profPDA(bvec, y, basisobj, basismat, Bmat, Dmat, ...
        bwtcell, awtcell, ufdcell, lambda);
    hessdg(ilam,:) = diag(inv(hessian))';
    [loglam(ilam), df(ilam), gcv(ilam)]
end

%  display RMS and penalized RMS

sqrt(SSE./((N*(n-df))))

%  display DF and GCV values

[df, gcv]

%  display coefficient vectors

[bmat]

bindex = 1:nlam;

bmat(bindex,:,:)

% output some results

for ilam=bindex    
    bvec = bmat(ilam,:)';
    
    %  set up linear differential operator
    
    [Lfdobj, bfdcell, afdcell] = ...
        bvec2Lfd(bvec, bwtcell, awtcell, ufdcell);
    
    lambda   = 10^loglam(ilam);
    fdParobj = fdPar(basisobj, Lfdobj, lambda);

     %  plot the operator
    
    plot(Lfdobj)
    title(['Lambda = ',num2str(lambda)]);
    pause;
    
    %  smooth the data with the operator
    
    [fdobj, dfi, gcvi, coef, SSEi, penmati] = ...
                                smooth_basis(y, tval, fdParobj);
    
    %  plot the data and fit
    
    subplot(1,1,1)
    yhat = eval_fd(tval, fdobj);
    yres = y - yhat;
    plot(tval, y, '.', tval, yhat, 'b-')
    xlabel('\fontsize{16} t')
    ylabel('\fontsize{16} x(t)')
    title(['Lambda = ',num2str(lambda)]);
    legend('\fontsize{16} Data', 'Estimate')
    pause;
    
    plot(tval, eval_fd(tval, fdobj, 1), 'b-')
    xlabel('\fontsize{16} t')
    ylabel('\fontsize{16} Dx(t)')
    title(['Lambda = ',num2str(lambda)]);
    pause;
end



