function ufdcell = getufd(Lfdobj)
%  GETUFD   Extracts the forcing function from LFDOBJ.

%  last modified 12 January 2003

if ~isa_Lfd(Lfdobj)
    error('Argument is not a linear differential operator object');
end

ufdcell = Lfdobj.ufdcell;


