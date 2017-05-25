function dy = derivcell(tnow, y, bfdcell) 
% DERIVCELL sets up the 1st order system corresponding to   
%   linear differential operator defined by BFD, UFD and AFD
%   returned by PDACELL.

% In this system, order varies inside variable inside equation.
% Example:  for two equations of the second order,
% y(1) is value of                     variable 1
% y(2) is value of first derivative of variable 1
% y(3) is value of                     variable 2
% y(4) is value of first derivative of variable 2

%  last modified 14 December 2004

bfddims = size(bfdcell);
nvar    = bfddims(1);
norder  = bfddims(length(bfddims));

%  for efficiency, there are two versions, one for a single
%  equation, and one for multiple equations

wmat = zeros(nvar*norder);

if nvar == 1

    %   -------------------------------------------------------
    %               Single equation
    %   -------------------------------------------------------

    for j=2:norder
        wmat(j-1,j) = 1;
    end
    m2 = 0;
    for i=1:nvar
        for j=1:norder
            m2 = m2 + 1;
            wmat(norder,m2) = -eval_fd(tnow, getfd(bfdcell{i,j}));
        end
    end
else

    %   -------------------------------------------------------
    %               Multiple equations
    %   -------------------------------------------------------

    m1 = 0;
    for i1=1:nvar
        for j=2:norder
            m1 = m1 + 1;
            wmat(m1,(i1-1)*norder + j) = 1;
        end
        m1 = m1 + 1;
        m2 = 0;
        for i2=1:nvar
            for j=1:norder
                m2 = m2 + 1;
                wmat(m1,m2) = -eval_fd(tnow, getfd(bfdcell{i1,i2,j}));
            end
        end
    end
end

dy = wmat*y;


