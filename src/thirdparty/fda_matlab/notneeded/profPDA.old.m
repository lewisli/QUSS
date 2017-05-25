function [PENSSE, SSE, penalty, cmat, rmat] = ...
    profPDA(bvec, y, t, phi_basis, psi_basis, norder, lambda)

npsi_basis = getnbasis(psi_basis);
%  set up wfd
wfd     = fd(reshape(bvec,npsi_basis,norder), psi_basis);
%  compute penalty matrix
Lfdobj  = Lfd(norder,fd2cell(wfd));
Kmat    = eval_penalty(phi_basis, Lfdobj);
%  compute profile SSE
phimat  = eval_basis(t, phi_basis);
Mmat    = phimat'*phimat + lambda.*Kmat;
cmat    = Mmat \ (phimat'*y);
rmat    = y - phimat*cmat;
SSE     = sum(sum(rmat.^2));
penalty = sum(diag(cmat'*Kmat*cmat));
PENSSE  = SSE + lambda.*penalty;

%  --------------------------------------------------------
%                all this code has been replaced
%  --------------------------------------------------------

% range   = getbasisrange(phi_basis);
% delta   = (range(2)-range(1))/(nfine-1);
% tfine   = linspace(range(1), range(2), nfine)';
% %  compute basis matrices
% phimat  = eval_basis(tfine, phi_basis);
% psimat  = eval_basis(tfine, psi_basis);
% Dphimat = eval_basis(tfine, phi_basis, 1);

% %  compute product of basis matrices
% phipsimat  = zeros([nfine, nbasis_phi, nbasis_psi]);
% for i=1:nbasis_phi
%     for j=1:nbasis_psi
%         phipsimat(:,i,j) = phimat(:,i).*psimat(:,j);
%     end
% end
% %  compute array of integral values for cross product of products
% Rarray = zeros([nbasis_phi, nbasis_psi, nbasis_phi, nbasis_psi]);
% for i1=1:nbasis_phi
%     for j1=1:nbasis_psi
%         for i2=1:nbasis_phi
%             for j2=1:nbasis_psi
%                 Rarray(i1,j1,i2,j2) = ...
%                     delta.*trapz(phipsimat(:,i1,j1).* ...
%                                  phipsimat(:,i2,j2));
%             end
%         end
%     end
% end
% %  compute array of integral values for product of
% %    product bases and first derivative of phi basis
% Sarray = zeros([nbasis_phi, nbasis_psi, nbasis_phi]);
% for i1=1:nbasis_phi
%     for j1=1:nbasis_psi
%         for i2=1:nbasis_phi
%             Sarray(i1,j1,i2) = ...
%                 delta.*trapz(phipsimat(:,i1,j1).*Dphimat(:,i2));
%         end
%     end
% end
% %  compute matrix of integrals of products of 1st derivs
% Tmat = zeros([nbasis_phi, nbasis_phi]);
% for i1=1:nbasis_phi
%     for i2=i1:nbasis_phi
%         Tmat(i1,i2) = ...
%             delta.*trapz(Dphimat(:,i1).*Dphimat(:,i2));
%         Tmat(i2,i1) = Tmat(i1,i2);
%     end
% end

% Kmat = zeros([nbasis_phi, nbasis_phi]);
% for i1=1:nbasis_phi
%     for i2=i1:nbasis_phi
%         Kmat(i1,i2) = bvec'*squeeze(Rarray(i1,:,i2,:))*bvec + ...
%                       2.0.*bvec'*squeeze(Sarray(i1,:,i2)) + ...
%                       Tmat(i1,i2);
%         Kmat(i2,i1) = Kmat(i1,i2);
%     end
% end
