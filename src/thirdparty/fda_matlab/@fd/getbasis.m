function basisobj = getbasis(fd)
%  GETBASIS    Extracts basis.fd object from functional data object.
%     FD, or, if FD is already a basis object, just returns it.

%  last modified 1 July 1998

  if isa_fd(fd) | isa_basis(fd)
    if isa_fd(fd)
      basisobj = fd.basisobj;
    else
      basisobj = fd;
    end
  else
    error('Structure is neither of fd type or of basis type');
  end

