function fd = putnames(fd, fdnames)
% PUTNAMES  Assigns fdnames to a functional data object FD

  if isa_fd(fd)
    if strcmp(class(fdnames), 'cell')
      fd.fdnames = fdnames;
    else
      error('Argument FDNAMES is not a cell object.');
    end
  else
    error('Argument FD is not a functional data object.');
  end

