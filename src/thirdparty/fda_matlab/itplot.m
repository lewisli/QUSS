load itdata.txt

s = itdata(:,2);
b = itdata(:,3);
f = itdata(:,4);
span = max(s);
n    = length(s);
d    = span/10;

plot(s, f, 'o')
for i=1:n
   si = s(i);
   fi = f(i);
   line([si,si+d],[fi,fi+b(i)*d])
end
