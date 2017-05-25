function display(fdParobj)
fprintf('\nFD:\n');
display(fdParobj.fd);
nderiv = getnderiv(fdParobj.Lfd);
if nderiv > 0
    fprintf('\nLFD:\n\n');
    display(fdParobj.Lfd);
else
    fprintf('\nLFD      = 0');
end
fprintf('\n\nLAMBDA   = %.6g', fdParobj.lambda);
fprintf('\nESTIMATE = %d', fdParobj.estimate);
fprintf('\n\nPENALTY MATRIX\n');
disp(fdParobj.penmat)

