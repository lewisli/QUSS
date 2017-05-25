function plusfd = uplus(fdobj)
% Unary plus of functional data object.

%  last modified 24 April 2003

if ~(isa_fd(fdobj)
    error('Argument is not a functional data object.');
end

plusfd  = fdobj;


