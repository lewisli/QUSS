range    = [0,100];
nbasis   = 103;
basisobj = create_bspline_basis(range, nbasis);

nquad = 2001;
delta = 100/(nquad-1);
quadvals = zeros(nquad,2);
quadvals(:,1) = linspace(0,100,nquad)';
quadvals(:,2) = ones(nquad,1);
quadvals(2:2:nquad-1,2) = 4;
quadvals(3:2:nquad-2,2) = 2;
quadvals(:,2) = delta.*quadvals(:,2)/3;

basisobj = putquadvals(basisobj,quadvals);

newquadvals = getquadvals(basisobj);

for ivalue=1:3
    basismat = eval_basis(quadvals(:,1), basisobj, ivalue-1);
    values{ivalue} = basismat.*(sqrt(quadvals(:,2))*ones(1,nbasis));
end

basisobj = putvalues(basisobj, values);

values = getvalues(basisobj);

basismat = getvalues(basisobj, 0);

Dbasismat = getvalues(basisobj, 1);

D2basismat = getvalues(basisobj,2);

penmat1 = D2basismat'*D2basismat;

penmat2 = eval_penalty(basisobj,2);

max(max(abs(full(penmat1)-full(penmat2))))/ ...
max(max(abs(full(penmat2))))

%  basis for Lfdobj

nbasisw  = 5;
basiswobj = create_bspline_basis(range, nbasisw);

%  create Lfdobj

wfd1 = fd(rand(nbasisw,1),basiswobj);
wfdcell{1} = fdPar(wfd1);
Lfdobj = Lfd(1, wfdcell);

wfd2 = fd(rand(nbasisw,1),basiswobj);
wfdcell{2} = fdPar(wfd2);
Lfdobj = Lfd(2, wfdcell);

onesbas = ones(1,nbasis);
tic
basismat   = getvalues(basisobj);
Dbasismat  = getvalues(basisobj,1);
% D2basismat = getvalues(basisobj,2);
toc

tic
basismat   =  basismat.*(eval_fd(quadvals(:,1), wfd1)*onesbas);
Lbasismat  = basismat + Dbasismat;
% Dbasismat  = Dbasismat.*(eval_fd(quadvals(:,1), wfd2)*onesbas);
% Lbasismat  = basismat + Dbasismat + D2basismat;
penmat1    = Lbasismat'*Lbasismat;
toc

tic
[penmat2, iter] = eval_penalty(basisobj,Lfdobj);
toc
iter

mean(mean(abs(full(penmat1)-full(penmat2))))/ ...
mean(mean(abs(full(penmat2))))



