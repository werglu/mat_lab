function[s] = calkowanie(n, m, f)
%Treœæ zadania: Obliczanie ca³ek ?f(?,?)dxdy po D = { (x, y): ?2 + ?2 ? 1 } przez transformacjê na prostok¹t [0,1]×[0, 2?] (wspó³rzêdne biegunowe) i zastosowanie z³o¿onych kwadratur Simpsona ze wzglêdu na ka¿d¹ zmienn¹.
%Dane wejsciowe:
%   n, m - ilosc podprzedzialow
%   f - funkcja dwóch zniennych (x, y)
%Dane wyjsciowe
% s - obliczona wartosc calki

C = zeros(2*n+1, 2*m+1);
c1(1)=1;
c2(1)=1;
c1(2*m+1)=1;
c2(2*n+1)=1;
pom=1;
for i=2:2*m
    if pom==1
        c1(i)=4;
        pom=0;
    else
        c1(i)=2;
        pom=1;
    end
end
pom=1;
for i=2:2*n
    if pom==1
        c2(i)=4;
        pom=0;
    else
        c2(i)=2;
        pom=1;
    end
end

C=c2'*c1;

%granice ca³kowanie we wspolrzednych biegunkowych:
a=0;
b=1;
c=0;
d=2*pi;

h1=(b-a)/(n);
h2=(d-c)/(m);

g=@(r, q)r*f(r*cos(q), r*sin(q)) %przeksztalcamy f(x, y) na wspol. biegunowe

x=linspace(a, b, 2*n+1);
y=linspace(c, d, 2*m+1);
s=0;

for i=1:2*n+1
   for j=1:2*m+1
     s = s + (((h1.*h2)./36).*C(i, j).*g(x(i),y(j))); 
   end
end

end

