function[s] = calkowanie_trapezy(n, m, f)
%Dane wejsciowe:
%   n, m - ilosc podprzedzialow
%   f - funkcja dwóch zniennych (x, y)
%Dane wyjsciowe
% s - obliczona wartosc calki

C = zeros(n+1, m+1);
c1(1)=1;
c2(1)=1;
c1(m+1)=1;
c2(n+1)=1;

for i=2:m
    c1(i)=2;
end

for i=2:n
   c2(i)=2;
end

C=c2'*c1;

%granice ca³kowanie we wspolrzednych biegunkowych:
a=0;
b=1;
c=0;
d=1;

h1=(b-a)/n;
h2=(d-c)/m;

x=linspace(a, b, n+1);
y=linspace(c, d, m+1);
s=0;

for i=1:n+1
   for j=1:m+1
     s = s + (C(i, j).*f(x(i),y(j))); 
   end
end
s=s*((h1.*h2)./4);
end
