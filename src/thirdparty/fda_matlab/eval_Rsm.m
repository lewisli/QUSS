function [penmatcell, Dpenmatcell] = ...
              eval_Rsm(npar, fitcell, gradwrd)
% EVAL_RS computes the matrices and vectors required to define the
%  penalty for an order one system of equations with forcing functions.
%  these are stored in cell arrays PENMATCELL and DPENMATCELL.
%  NCOEFS is the total number of coefficients defining the fits
%  for the M variables, that is NBASIS1 + ... + NBASISM.
%  If DERIVS is positive, derivatives of these vectors are also
%  computed

%  Last modified 10 August 2004

if nargin < 4
    gradwrd = 1;
end

nvar   = length(fitcell);
nbasisvec = zeros(nvar,1);
for ivar=1:nvar
    fitstruct = fitcell{ivar};
    basisobji = fitstruct.basisobj;
    nbasis    = getnbasis(basisobji);
    nbasisvec(ivar) = nbasis;
end
ncoefs = sum(nbasisvec);

%  loop through variables

mi2  = 0;
for ivar = 1:nvar
    fitstruct = fitcell{ivar};
    basisobji = fitstruct.basisobj;
    nbasis    = getnbasis(basisobji);
    mi1       = mi2 + 1;
    mi2       = mi2 + nbasis;
    indi      = mi1:mi2;
    
    %  set up the weight and forcing function cells
    
    bwtcell = fitstruct.bwtcell;
    awtcell = fitstruct.awtcell;
    ufdcell = fitstruct.ufdcell;
    if isempty(awtcell) | isempty(ufdcell)
        nforce = 0;
    else
        nforce = length(awtcell);
    end
    
    %  set up arrays to be filled
    
    Rmat = zeros(nbasis,nbasis);
    Smat = zeros(ncoefs,ncoefs);
    Tmat = zeros(nbasis,ncoefs);
    if nforce > 0
        Umat = zeros(nbasis,1);
        Vmat = zeros(ncoefs,1);
        Wmat = 0;
    else
        Umat = [];
        Vmat = [];
        Wmat = [];
    end
    
    %  retrieve quadrature points and weights
    
    quadvals = getquadvals(basisobji);
    tvalquad = quadvals(:,1);
    
    %  retrieve weighted basis function values
    
    basismati = getvalues(basisobji);
    if ivar == 1
        Dorder = fitstruct.Dorder;
    end

    %  Compute penalty matrix Rmat, penalty vector Umat, and
    %  penalty scalar Wmat that depend only on variable IVAR
    %  compute L_phi  values for this variable
    
    Dmbasismati = getvalues(basisobji, Dorder);
    Rmat = Dmbasismati'*Dmbasismati;
    
    for k=1:nforce
        auveck = eval_au(tvalquad, awtcell{k}, ufdcell{k});
        Umat   = Umat + Dbasismati'*auveck;
        Wmat   = Wmat + sum(auveck.^2);
    end
    
    %  loop through other variables in this equation to compute
    %  penalty matrices Smat and Tmat
    
    mj2 = 0;
    for jvar=1:nvar
        fitstructj = fitcell{jvar};
        basisobjj  = fitstructj.basisobj;
        nbasisj    = getnbasis(basisobjj);
        mj1  = mj2 + 1;
        mj2  = mj2 + nbasisj;
        indj = mj1:mj2;
        
        %  compute L_phi  values for this variable
    
        tempj = eval_Lphi(jvar, bwtcell, basisobjj);
       
        %  Compute penalty matrix Tmat entries for this J
        
        Tmat(:,indj) = Dmbasismati'*tempj;
        
        %  loop through all over variables in this equation 
        %  to compute penalty matrix Smat entries for fixed J

        mk2 = 0;
        for kvar=1:jvar
            fitstructk = fitcell{kvar};
            basisobjk  = fitstructk.basisobj;
            nbasisk    = getnbasis(basisobjk);
            mk1  = mk2 + 1;
            mk2  = mk2 + nbasisk;
            indk = mk1:mk2;
            
            %  compute L_phi  values for this variable
            
            tempk = eval_Lphi(kvar, bwtcell, basisobjk);
            
            Smatkj = tempk'*tempj;
            Smat(indk,indj) = Smatkj;
            if kvar ~= jvar
                Smat(indj,indk) = Smatkj';
            end
        end
        
        %  Loop through forcing functions to compute 
        %  penalty matrix Vmat for this J that sums over
        %  forcing functions
        
        for k=1:nforce
            auveck    = eval_au(tvalquad, awtcell{k}, ufdcell{k});
            basismatj = getvalues(basisobjj, 1);
            tempj     = basismatj'*(auveck.*bvecj);
            for jderiv = 1:Dorderj-1
                Dbfdj       = getfd(bwtcell{jvar},jderiv+1);
                Dbbasisj    = getbasis(bfdj);
                Dbcoefj     = getcoef(bfdj);
                Dbbasismatj = getvalues(bbasisj,jderiv+1);
                Dbvecj      = bbasismatj*Dbcoefj;
                Dbasismatj  = getvalues(basisobjj,jderiv+1);
                tempj       = tempj + Dbasismatj'*(auveck.*Dbvecj);
            end
            Vmat(indj) = Vmat(indj) + tempj;
        end
    end
    
    %  place these matrices in struct PENSTRUCT
    
    penstruct.Rmat = Rmat;
    penstruct.Smat = Smat;
    penstruct.Tmat = Tmat;
    penstruct.Umat = Umat;
    penstruct.Vmat = Vmat;
    penstruct.Wmat = Wmat;
    
    penmatcell{ivar} = penstruct;
    
end

%  ------------------------------------------------------------
%                 evaluate the derivatives
%  ------------------------------------------------------------

% if gradwrd
%     
%     %  loop through equations
%     
%     m2  = 0;
%     mi2 = 0;
%     for ivar=1:nvar
%         fitstruct = fitcell{ivar};
%         penstruct = penmatcell{ivar};
%         basisobj  = fitstruct.basisobj;
%         
%         %  set up the weight and forcing function cells
%         
%         bwtcell = fitstruct.bwtcell;
%         awtcell = fitstruct.awtcell;
%         
%         %  retrieve quadrature points and weights
%         
%         quadvals = getquadvals(basisobj);
%         tvalquad = quadvals(:,1);
%         
%         %  retrieve weighted basis function values
%         
%         D0basismat = getvalues(basisobj);
%         D1basismat = getvalues(basisobj, 1);
%         nbasis     = nbasisvec(ivar);
%         onebas     = ones(1,nbasis);
%         
%         %  set up arrays to be filled
%         
%         nbcoefs = sum(npar(ivar,1:Dorder));
%         nacoefs = npar(ivar,Dorder+1);
%         DSmat      = zeros(ncoefs,ncoefs,nbcoefs);
%         DTmat      = zeros(nbasis,ncoefs,nbcoefs);
%         DUmat      = zeros(nbasis,       nacoefs);
%         DaVmat     = zeros(ncoefs,       nacoefs);
%         DbVmat     = zeros(ncoefs,       nbcoefs);
%         nforce     = length(awtcell);
%         
%         mi1 = mi2 + 1;
%         mi2 = mi2 + nbasis;
%         
%         %  loop through variables and derivatives in this equation
%         
%         mj2   = 0;
%         nbpar = 0;
%         for jvar=1:nvar
%             fitstructj = fitcell{jvar};
%             basisobjj  = fitstructj.basisobj;
%             basismatj  = getvalues(basisobjj);
%             onebasj    = ones(1,nbasisvec(jvar));
%             mj1        = mj2 + 1;
%             mj2        = mj2 + nbasisvec(jvar);
%             indj       = mj1:mj2;
%             for ideriv=1:Dorder
%                 bfdParj = bwtcell{jvar,ideriv};
%                 if getestimate(bfdParj) 
%                     bfdj       = getfd(bwtcell{jvar,ideriv});
%                     bbasisj    = getbasis(bfdj);
%                     bbasismatj = getvalues(bbasisj);
%                     nbbasisj   = getnbasis(bbasisj);
%                     
%                     %  Compute DSmat by looping through all variables and
%                     %    derivatives for fixed variable J
%                     
%                     m1  = m2 + 1;
%                     m2  = m2 + nbbasisj;
%                     mk2 = 0;
%                     for kvar=1:nvar
%                         fitstructk = fitcell{kvar};
%                         basisobjk  = fitstructk.basisobj;
%                         basismatk  = getvalues(basisobjk);
%                         onebask    = ones(1,nbasisvec(kvar));
%                         mk1   = mk2 + 1;
%                         mk2   = mk2 + nbasisvec(kvar);
%                         indk  = mk1:mk2;
%                         tempk = eval_Lphi(kvar, bwtcell, basisobjk);
%                         for kderiv=1:Dorder
%                             for m=m1:m2
%                                 bfdvecm      = bbasismatj(:,m-m1+1);
%                                 bDjbasismatm = basismatj.*(bfdvecm*onebasj);
%                                 prodjkm      = bDjbasismatm'*tempk;
%                                 DSmat(indj,indk,indm) = DSmat(indj,indk,nbpar+m-m1+1) + ...
%                                     prodjkm;
%                                 DSmat(indk,indj,indm) = DSmat(indk,indj,nbpar+m-m1+1) + ...
%                                     prodjkm';
%                             end
%                         end
%                     end
%                     
%                     %  compute DTmat and DbVmat for variable J
%                     
%                     for m=m1:m2
%                         bfdvecm      = bbasismatj(:,m-m1+1);
%                         bDjbasismatm = basismatj.*(bfdvecm*onebasj);
%                         prod         = D1basismat'*bDjbasismatm;
%                         DTmat(:,indj,nbpar+m-m1+1) = prod;
%                         %  set up the part of the roughness penalty 
%                         %  affected by the presence of forcing function(s)
%                         Dpenvecm = zeros(nbasisvec(jvar),1);
%                         if nforce > 0
%                             for k=1:nforce
%                                 auveck    = eval_au(tvalquad, awtcell{k}, ufdcell{k});
%                                 Dpenvecmk = bDjbasismatm'*auveck;
%                                 Dpenvecm  = Dpenvecm + Dpenvecmk;
%                             end
%                             DbVmat(indj,nbpar+m-m1+1) = Dpenvecm;
%                         end
%                     end
%                     nbpar = nbpar + nbbasisj;
%                 end
%             end
%         end
%         
%         %  update penalized least squares for terms for 
%         %  smoothing forcing function coefficients
%         
%         napar = 0;
%         for k=1:nforce
%             afdPark = awtcell{k};
%             if getestimate(afdPark) 
%                 abasisk    = getbasis(getfd(afdPark));
%                 abasismatk = getvalues(abasisk);
%                 anbasis    = getnbasis(abasisk);
%                 ufdk       = ufdcell{k};
%                 uveck      = eval_fd(tvalquad, ufdk);
%                 m1 = m2 + 1;
%                 m2 = m2 + anbasis;
%                 for m=m1:m2
%                     napar    = napar + 1;
%                     Daveck   = abasismatk(:,m-m1+1);
%                     Dauveck  = Daveck.*uveck.*sqrt(quadvals(:,2));
%                     Dpenvecm = D1basismat'*Dauveck;
%                     DUmat(:,napar) = Dpenvecm;
%                     %  compute DaVmat by looping through variables
%                     %  for fixed forcing function K
%                     mk2 = 0;
%                     for kvar=1:nvar
%                         fitstructk = fitcell{kvar};
%                         basisobjk  = fitstructk.basisobj;
%                         mk1   = mk2 + 1;
%                         mk2   = mk2 + nbasisvec(kvar);
%                         indk  = mk1:mk2;
%                         tempk = eval_Lphi(kvar, bwtcell, basisobjk);
%                         tempk = tempk'*Dauveck;
%                         DaVmat(indk,napar) = tempk;
%                     end
%                 end
%             end
%         end  
%         Dpenstruct.DSmat  = DSmat;
%         Dpenstruct.DTmat  = DTmat;
%         Dpenstruct.DaVmat = DaVmat;
%         Dpenstruct.DbVmat = DbVmat;
%         Dpenstruct.DUmat  = DUmat;
%         Dpenmatcell{ivar} = Dpenstruct;
%         
%     end
%     
% else
    Dpenmatcell = {};
% end

%  ----------------------------------------------------------------

%  --------------------------------------------------------------

