function [y,dy] = polintarray(xa, ya, x) 
  %  YA is an array with up to 4 dimensions
  %     with 1st dim the same length same as the vector XA
  n     = length(xa);
  yadim = size(ya);
  nydim = length(yadim);
  if yadim(2) == 1, nydim = 1; end
  if yadim(1) ~= n, error('First dimension of YA must match XA'); end
  difx = xa - x;
  absxmxa = abs(difx);
  tmp = 1:n;
  ns = min(tmp(absxmxa == min(absxmxa)));
  cs = ya;
  ds = ya;
  if nydim == 1, y = ya(ns);  end
  if nydim == 2, y = ya(ns,:);  end
  if nydim == 3, y = ya(ns,:,:);  end
  if nydim == 4, y = ya(ns,:,:,:);  end
  ns = ns - 1;
  for m = 1:(n-1)
    if nydim == 1
      for i = 1:(n-m)
        ho      = difx(i);
        hp      = difx(i+m);
        w       = (cs(i+1) - ds(i))./(ho - hp);
        ds(i) = hp.*w;
        cs(i) = ho.*w;
      end
      if 2*ns < n-m
        dy = cs(ns+1);
      else 
        dy = ds(ns);
        ns = ns - 1;
      end
    end
    if nydim == 2
      for i = 1:(n-m)
        ho      = difx(i);
        hp      = difx(i+m);
        w       = (cs(i+1,:) - ds(i,:))./(ho - hp);
        ds(i,:) = hp.*w;
        cs(i,:) = ho.*w;
      end
      if 2*ns < n-m
        dy = cs(ns+1,:);
      else  
        dy = ds(ns,:);
        ns = ns - 1;
      end
    end
    if nydim == 3
      for i = 1:(n-m)
        ho      = difx(i);
        hp      = difx(i+m);
        w       = (cs(i+1,:,:) - ds(i,:,:))./(ho - hp);
        ds(i,:,:) = hp.*w;
        cs(i,:,:) = ho.*w;
      end
      if 2*ns < n-m
        dy = cs(ns+1,:,:);
      else
        dy = ds(ns,:,:);
        ns = ns - 1;
      end
    end
    if nydim == 4
      for i = 1:(n-m)
        ho      = difx(i);
        hp      = difx(i+m);
        w       = (cs(i+1,:,:,:) - ds(i,:,:,:))./(ho - hp);
        ds(i,:,:,:) = hp.*w;
        cs(i,:,:,:) = ho.*w;
      end
      if 2*ns < n-m
        dy = cs(ns+1,:,:,:);
      else 
        dy = ds(ns,:,:,:);
        ns = ns - 1;
      end
    end
    y = y + dy;
  end
