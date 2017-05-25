function auvec = eval_au(tvalquad, awtcell, ufd)
%  evaluate weight function times forcing function
afd       = getfd(awtcell);
abasis    = getbasis(afd);
acoef     = getcoef(afd);
abasismat = getvalues(abasis);
avec      = abasismat*acoef;
uvec      = eval_fd(tvalquad, ufd);
auvec     = avec.*uvec.*sqrt(quadvals(:,2));
