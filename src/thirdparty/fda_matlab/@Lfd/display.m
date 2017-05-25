function display(Lfd)
nderiv = Lfd.nderiv;
fprintf(['NDERIV = ', num2str(nderiv),'\n']);
if nderiv > 0
    fprintf('\nWFD:');
    bwtcell = Lfd.bwtcell;
    fprintf('\n\n-------------------');
    for ideriv=1:nderiv
        fprintf(['\n\nWFD(',num2str(ideriv-1),') fdPar object:\n'])
        display(Lfd.bwtcell{ideriv});
        fprintf('\n\n-------------------');
    end
end
% these statements removed because Lfd class is not
% going to containing forcing functions and their coefficients.
% fprintf('\nAFD:\n');
% display(Lfd.awtcell);
% fprintf('\nUFD:\n');
% display(Lfd.ufdcell);

