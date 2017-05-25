
function fdtype = use_proper_basis(fdtype)
%  USE_PROPER_BASIS recognizes type of basis by use of several variant spellings

%  Last modified 24 October 2003

switch fdtype
    
    case 'Fourier'
        fdtype = 'fourier';
    case 'fourier'
        fdtype = 'fourier';
    case 'Fou'
        fdtype = 'fourier';
    case 'fou'
        fdtype = 'fourier';
        
    case 'bspline'
        fdtype = 'bspline';
    case 'Bspline'
        fdtype = 'bspline';
    case 'Bsp'
        fdtype = 'bspline';
    case 'bsp'
        fdtype = 'bspline';
        
    case 'power'
        fdtype = 'power';
    case 'pow'
        fdtype = 'power';
        
    case 'polyg'
        fdtype = 'polyg';
    case 'polygon'
        fdtype = 'polyg';
    case 'polygonal'
        fdtype = 'polyg';
        
    case 'exp'
        fdtype = 'expon';
    case 'expon'
        fdtype = 'expon';
    case 'exponential'
        fdtype = 'expon';
        
    case 'con'
        fdtype = 'const';
    case 'const'
        fdtype = 'const';
    case 'const'
        fdtype = 'const';
        
    case 'mon'
        fdtype = 'monom';
    case 'monom'
        fdtype = 'monom';
    case 'monomial'
        fdtype = 'monom';
        
    case 'poly'
        fdtype = 'polynom';
    case 'polynom'
        fdtype = 'polynom';
    case 'polynomial'
        fdtype = 'polynom';
        
    otherwise
        fdtype = 'unknown';
        
end
