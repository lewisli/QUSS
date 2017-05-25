function [penmat, penvec, DR, Ds] = ...
               eval_Rs(bwtcell, awtcell, ufdcell, basisobj, gradwrd)
% EVAL_RS computes the penalty matrix R(bvec) 
%  and penalty vector s(bvec) where BVEC is a vector of parameters
%  defining the weight coefficients in WFDCELL and AFDCELL
%  If DERVIVS is positive, derivatives of these vectors are also
%  computed

%  Last modified 10 August 2004

if nargin < 5
    gradwrd = 1;
end

nderiv = length(bwtcell);
nbasis = getnbasis(basisobj);

%  set up the homogeneous operator

Lfdobj = Lfd(nderiv, bwtcell);

%  compute the penalty matrix

penmat = eval_penalty(basisobj, Lfdobj);

%  evaluate the penalty vector if required

forcewrd = ~isempty(awtcell) & ~isempty(ufdcell);
if forcewrd
    %  set up the part of the roughness penalty affected by the
    %  presence of forcing function(s)
    nforce = length(ufdcell);
    penvec = zeros(nbasis,1);
    for k=1:nforce
        afdk    = getfd(awtcell{k});
        acoefk  = getcoef(afdk);
        abasisk = getbasis(afdk);
        ufdk    = ufdcell{k};
        ucoefk  = getcoef(ufdk);
        ubasisk = getbasis(ufdk);
        coefmat = acoefk*ucoefk';
        aubifdk = bifd(coefmat, abasisk, ubasisk);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        penveck = inprod(basisobj, aubifdk, Lfdobj);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        penvec  = penvec - penveck; 
        if gradwrd
            ffdcell{k} = aubifdk;
        end
    end
else
    nforce = 0;
    penvec = [];
    Ds     = [];
end

if gradwrd
    
    %  evaluate the derivatives

    m2 = 0;
    fdbasis = fd(eye(nbasis),basisobj);
    for j=1:nderiv
        Lfdobjj = int2Lfd(j-1);
        bfdParj = bwtcell{j};
        if getestimate(bfdParj) 
            bfdj    = getfd(bfdParj);
            basisj  = getbasis(bfdj);
            nbasisj = getnbasis(basisj);
            coefj   = getcoef(bfdj);
            m1 = m2 + 1;
            m2 = m2 + nbasisj;
            for m=m1:m2
                coefm = zeros(nbasisj,1);
                coefm(m-m1+1) = 1;
                coefmat = ...
                    reshape(kron(eye(nbasis),coefm),nbasisj,nbasis,nbasis);
                wDbasisbifd = bifd(coefmat, basisj, basisobj);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Dpenmatm = inprod(wDbasisbifd, basisobj, Lfdobjj, Lfdobj);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                DR(1:nbasis,1:nbasis,m) = Dpenmatm + Dpenmatm';
                if forcewrd
                    %  set up the part of the roughness penalty a
                    %  ffected by the presence of forcing function(s)
                    Dpenvecm = zeros(nbasis,1);
                    for k=1:nforce
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        Dpenvecmk = inprod(wDbasisbifd, ffdcell{k}, 0, 0);
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        Dpenvecm = Dpenvecm - Dpenvecmk;
                    end
                    Ds(1:nbasis,m) = Dpenvecm;
                end
                
                %                 lambdaj = getlambda(bfdParj);
                %                 if lambdaj > 0
                %                     Lfdobjj = getLfd(bfdParj);
                %                     penmatj = eval_penalty(basisj, Lfdobjj);
                %                     bcoefj  = bvec(m1:m2);
                %                     termj   = lambdaj.*bcoefj'*penmatj*bcoefj;
                %                     SSE     = SSE + termj;
                %                 end
            end
        end
    end
    
    %  update penalized least squares for terms for 
    %  smoothing forcing function coefficients
    
    for k=1:nforce
        afdPark = awtcell{k};
        if getestimate(afdPark) 
            afdk    = getfd(afdPark);
            abasisk = getbasis(afdk);
            nbasisk = getnbasis(abasisk);
            ufdk    = ufdcell{k};
            ucoefk  = getcoef(ufdk);
            ubasisk = getbasis(ufdk);
            m1 = m2 + 1;
            m2 = m2 + nbasisk;
            for m=m1:m2
                acoefk = zeros(nbasisk,1);
                acoefk(m-m1+1) = 1;
                coefmat = acoefk*ucoefk';
                Daubifdk = bifd(coefmat, abasisk, ubasisk);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Dpenvecm = -inprod(basisobj, Daubifdk, Lfdobj, 0);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Ds(1:nbasis,m) = Dpenvecm;
                
                %             lambdak = getlambda(afdPark);
                %             if lambdak > 0
                %                 Lfdobjk = getLfd(afdPark);
                %                 penmatk = eval_penalty(basisk, Lfdobjk);
                %                 bcoefk  = bvec(m1:m2);
                %                 SSE     = SSE + lambdak.*bcoefk'*penmatk*bcoefk;
                %             end
                
            end
        end
    end   
else
    DR = [];
    Ds = [];
end